"""Analyze actual Exodus fields from the controlled lookup/quadratic runs."""
from pathlib import Path
import argparse, json, html, hashlib, base64, re
import numpy as np
import pandas as pd
import netCDF4
from scipy.interpolate import CubicHermiteSpline
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/ni_al_lookup'
parser=argparse.ArgumentParser()
parser.add_argument('--partial',action='store_true')
parser.add_argument('--no-gif',action='store_true')
parser.add_argument('--output',default='results/ni_al_lookup')
args=parser.parse_args()
OUT=ROOT/args.output
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

def read(path):
    with netCDF4.Dataset(path) as ds:
        names=list(netCDF4.chartostring(ds['name_nod_var'][:]))
        x,y=np.array(ds['coordx'][:]),np.array(ds['coordy'][:])
        order=np.lexsort((x,y)); nx=len(np.unique(x)); ny=len(np.unique(y))
        t=np.array(ds['time_whole'][:]); assert np.all(np.diff(t)>0)
        fields={}
        for name in ['c','eta','mu']:
            a=ds[f'vals_nod_var{names.index(name)+1}'][:]
            assert not np.ma.getmaskarray(a).any()
            fields[name]=np.array(a)[:,order].reshape(-1,ny,nx)
        for a in fields.values(): assert np.isfinite(a).all()
        elem={}
        if 'name_elem_var' in ds.variables:
            enames=list(netCDF4.chartostring(ds['name_elem_var'][:]))
            conn=np.array(ds['connect1'][:])-1
            cx=x[conn].mean(axis=1); cy=y[conn].mean(axis=1)
            eo=np.lexsort((cx,cy))
            for name in ['c_gamma','c_gamma_prime']:
                if name in enames:
                    ev=ds[f'vals_elem_var{enames.index(name)+1}eb1'][:]
                    assert not np.ma.getmaskarray(ev).any()
                    elem[name]=np.array(ev)[:,eo].reshape(-1,ny-1,nx-1)
                    assert np.isfinite(elem[name]).all()
        return t,fields,elem

def track(t,f):
    eta,c=f['eta'],f['c']; h=eta**3*(6*eta**2-15*eta+10)
    n=eta.shape[-1]; w=np.ones((n,n)); w[[0,-1],:]*=.5; w[:,[0,-1]]*=.5
    step=96/(n-1); m=(n-1)//2
    slices=[(slice(0,m+1),slice(0,m+1)),(slice(0,m+1),slice(m,n)),(slice(m,n),slice(0,m+1)),(slice(m,n),slice(m,n))]
    rows=[]
    for i,time in enumerate(t):
        row={'time_star':time,'mean_c':float((c[i]*w).sum()*step**2/96**2),'area_fraction':float((h[i]*w).sum()*step**2/96**2), 'count':label(eta[i]>=.5)[1]}
        for j,(ys,xs) in enumerate(slices,1):
            q=h[i,ys,xs]; qw=np.ones_like(q); qw[[0,-1],:]*=.5; qw[:,[0,-1]]*=.5
            row[f'R{j}_nm']=np.sqrt(max(0,(q*qw).sum()*step**2)/np.pi)
            row[f'present{j}']=int(eta[i,ys,xs].max()>=.5)
        rows.append(row)
    return pd.DataFrame(rows)

runs={}
metrics={}
stems={'quadratic_control':'quadratic_matched' if (OUT/'quadratic_matched.e').exists() else 'quadratic_control','lookup':'lookup'}
metrics['result_file_stems']=stems
for key in ['quadratic_control','lookup']:
    t,f,elem=read(OUT/f'{stems[key]}.e'); tr=track(t,f)
    post=pd.read_csv(OUT/f'{stems[key]}.csv')
    # Live runs may have CSV entries beyond the last saved field; final runs must match.
    if not args.partial: assert abs(t[-1]-60000)<1e-6 and abs(post.time.iloc[-1]-60000)<1e-6
    tr.to_csv(OUT/f'{key}_particle_history.csv',index=False)
    err=float(abs(tr.mean_c-.14).max()); assert err<1e-8
    dE=np.diff(post.total_energy)
    assert max(0,dE.max(initial=0))<1e-5
    eta=f['eta']; assert eta.min()>-.01 and eta.max()<1.01
    m={'frames':len(t),'steps':len(post)-1,'last_time_star':float(t[-1]),'max_mass_error':err,'energy_initial':float(post.total_energy.iloc[0]),'energy_final':float(post.total_energy.iloc[-1]),'max_energy_increase':float(max(0,dE.max(initial=0))), 'radii_initial_nm':tr[[f'R{j}_nm' for j in range(1,5)]].iloc[0].tolist(),'radii_final_nm':tr[[f'R{j}_nm' for j in range(1,5)]].iloc[-1].tolist(), 'final_count':int(tr['count'].iloc[-1]),'eta_range':[float(eta.min()),float(eta.max())],'phase_composition_ranges':{},'dissolution_brackets_star':{}}
    for j in range(1,5):
        zero=np.flatnonzero(tr[f'present{j}'].to_numpy()==0)
        m['dissolution_brackets_star'][str(j)]=[float(t[zero[0]-1]),float(t[zero[0]])] if len(zero) else None
    center_eta=(eta[:,:-1,:-1]+eta[:,1:,:-1]+eta[:,:-1,1:]+eta[:,1:,1:])/4
    core=pd.DataFrame({'time_star':t})
    for name,a in elem.items():
        m['phase_composition_ranges'][name]=[float(a.min()),float(a.max())]
        mask=center_eta<.01 if name=='c_gamma' else center_eta>.99
        m[name+'_final_core_mean']=float(a[-1][mask[-1]].mean())
        m[name+'_core_range']=[float(a[mask].min()),float(a[mask].max())]
        core[name+'_mean']=[float(a[i][mask[i]].mean()) for i in range(len(t))]
    core.to_csv(OUT/f'{key}_core_compositions.csv',index=False)
    for j in range(1,5):
        for lo,hi in [(5000,30000),(10000,40000)]:
            if t[-1]>=hi:
                tx=np.linspace(lo,hi,251)
                r=np.interp(tx,t,tr[f'R{j}_nm'])
                slope,intercept=np.polyfit(tx,r,1)
                m[f'P{j}_dR_dt_{lo}_{hi}']=float(slope)
                slope3,intercept3=np.polyfit(tx,r**3,1)
                m[f'P{j}_dR3_dt_{lo}_{hi}']=float(slope3)
    metrics[key]=m
    runs[key]=(t,f,elem,tr,post)

bt,bf,_=read(ROOT/'results/ni_al_coarsening/coarsening.e')
qt,qf,_,qtr,qpost=runs['quadratic_control']
lt,lf,_,ltr,lpost=runs['lookup']
control={}
original_qt,original_qf,_=read(OUT/'quadratic_control.e')
for name in ['c','eta','mu']:
    errors=[]
    for i,time in enumerate(original_qt):
        j=int(np.argmin(abs(bt-time)))
        if abs(bt[j]-time)<1e-7: errors.append(float(abs(original_qf[name][i]-bf[name][j]).max()))
    control[name+'_max_abs_error']=max(errors)
metrics['baseline_control']=control
metrics['identical_initial_c_eta']=all(np.array_equal(qf[n][0],lf[n][0]) for n in ['c','eta'])
assert metrics['identical_initial_c_eta']
metrics['same_saved_time_grid']=bool(len(qt)==len(lt) and np.allclose(qt,lt,rtol=0,atol=1e-7))
reference_times=pd.read_csv(OUT/'lookup.csv').time.to_numpy()
metrics['time_sequence_exact']=all(len(runs[k][4])<=len(reference_times) and np.allclose(runs[k][4].time,reference_times[:len(runs[k][4])],rtol=0,atol=1e-7) for k in runs)
if not args.partial:
    qearly=qpost.time[qpost.time<=40000].to_numpy();learly=lpost.time[lpost.time<=40000].to_numpy()
    assert np.allclose(qearly,learly,rtol=0,atol=1e-7), 'Rate comparison windows require identical timesteps'
    metrics['identical_time_steps_through_40000']=True
    assert max(control.values())<1e-6, 'Nested quadratic control differs from the preserved baseline'
    for key in runs:
        for phase,bounds in [('c_gamma',(.05,.30)),('c_gamma_prime',(.18,.30))]:
            actual=metrics[key]['phase_composition_ranges'][phase]
            assert bounds[0]<=actual[0]<=actual[1]<=bounds[1]
    metrics['validation_status']='PASS'
else:
    metrics['validation_status']='PARTIAL'
thermo=json.loads((OUT/'thermodynamic_validation.json').read_text())
metrics['thermodynamics']=thermo
for k in ['quadratic_control','lookup']:
    if not args.partial: metrics[k]['exodus_sha256']=hashlib.sha256((OUT/f'{stems[k]}.e').read_bytes()).hexdigest()
if (OUT/'matched_time_control.json').exists():
    metrics['matched_time_control']=json.loads((OUT/'matched_time_control.json').read_text())
    metrics['quadratic_time_refinement_final_max_field_change']={n:float(abs(qf[n][-1]-original_qf[n][-1]).max()) for n in qf}
    original_track=track(original_qt,original_qf)
    metrics['quadratic_time_refinement_final_radius_change_nm']=[float(qtr[f'R{j}_nm'].iloc[-1]-original_track[f'R{j}_nm'].iloc[-1]) for j in range(1,5)]
if (OUT/'assembled_control.json').exists(): metrics['assembled_control']=json.loads((OUT/'assembled_control.json').read_text())
if not args.partial:
    original_hash=hashlib.sha256((ROOT/'results/ni_al_coarsening/coarsening.e').read_bytes()).hexdigest()
    recorded=json.loads((ROOT/'results/ni_al_coarsening/exodus_sha256.json').read_text(encoding='utf-8-sig'))['Hash'].lower()
    assert original_hash==recorded
    metrics['original_exodus_unchanged_sha256']=original_hash
(OUT/'comparison_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
summary=[]
for i,phase in enumerate(['gamma','gamma_prime']):
    a=100*thermo['quadratic_equilibrium_x_Al'][i]; b=100*thermo['common_tangent_x_Al'][i]
    summary.append({'metric':phase+'_planar_equilibrium_Al','unit':'at.%','quadratic':a,'lookup':b,'lookup_minus_quadratic':b-a,'relative_difference_percent':100*(b-a)/a})
for key in metrics['quadratic_control']:
    if ('_dR_dt_' in key or '_dR3_dt_' in key) and key in metrics['lookup']:
        a=metrics['quadratic_control'][key]; b=metrics['lookup'][key]
        summary.append({'metric':key,'unit':'nm^3/t*' if '_dR3_' in key else 'nm/t*','quadratic':a,'lookup':b,'lookup_minus_quadratic':b-a,'relative_difference_percent':100*(b-a)/a})
pd.DataFrame(summary).to_csv(OUT/'comparison_summary.csv',index=False)

def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=180)
    fig.savefig(OUT/f'{name}.pdf')
    plt.close(fig)

# Show the actual separate phase branches, not the equilibrium envelope.
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
for j,(phase,ceq,A,limits) in enumerate([('gamma',.124644,347973,(.123,.140)),('gamma_prime',.237162,424347.4,(.235,.25))]):
    df=pd.read_csv(ROOT/f'thermodynamics/ni_al_1000K/{phase}_1000K.csv')
    s=CubicHermiteSpline(df.x_Al,df.Gm_reduced_J_mol,df.dGm_reduced_dx_J_mol)
    xx=np.linspace(*limits,801)
    axs[0,j].plot(xx*100,.5*A*(xx-ceq)**2,'--',label='Quadratic')
    axs[0,j].plot(xx*100,s(xx),label='Table / cubic Hermite')
    sel=df.x_Al.between(*limits)
    axs[0,j].scatter(df.x_Al[sel]*100,df.Gm_reduced_J_mol[sel],s=10,c='black',alpha=.4,label='TC exported nodes')
    axs[1,j].plot(xx*100,A*(xx-ceq),'--',label='Quadratic')
    axs[1,j].plot(xx*100,s(xx,1),label='Table')
    for ax in axs[:,j]: ax.axvline(ceq*100,c='gray',lw=.7); ax.set(xlabel='Al (at.%)',title=phase); ax.legend(fontsize=8)
    axs[0,j].set_ylabel('Reduced molar free energy (J/mol)')
    axs[1,j].set_ylabel('Reduced diffusion potential (J/mol)')
save(fig,'free_energy_comparison')

# Expose the scale of equilibrium differences without truncating a bar axis.
fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
eq_q=np.array(thermo['quadratic_equilibrium_x_Al'])*100
eq_l=np.array(thermo['common_tangent_x_Al'])*100
xx=np.arange(2)
axs[0].bar(xx-.18,eq_q,.36,label='Quadratic')
axs[0].bar(xx+.18,eq_l,.36,label='Table')
axs[0].set(xticks=xx,xticklabels=['gamma','gamma-prime'],ylabel='Equilibrium Al (at.%)',ylim=(0,27))
axs[0].legend()
axs[1].bar(xx,(eq_l-eq_q)*1e6,color=['#0072b2','#d55e00'])
axs[1].axhline(0,color='black',lw=.7)
axs[1].set(xticks=xx,xticklabels=['gamma','gamma-prime'],ylabel='Table minus quadratic (10^-6 at.% points)')
fig.suptitle('Planar common tangent: differences are within export rounding precision')
save(fig,'equilibrium_comparison')

fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for j,(name,title) in enumerate([('c_gamma','Matrix core: eta < 0.01'),('c_gamma_prime','Precipitate core: eta > 0.99')]):
    for key,style,lab in [('quadratic_control','--','Quadratic'),('lookup','-','Table')]:
        core=pd.read_csv(OUT/f'{key}_core_compositions.csv')
        axs[j].plot(core.time_star,core[name+'_mean']*100,style,label=lab)
    axs[j].axhline(eq_l[j],color='gray',lw=.8,label='Planar equilibrium')
    axs[j].set(xlabel='t*',ylabel='Mean core Al (at.%)',title=title)
    axs[j].legend(fontsize=8)
save(fig,'core_compositions_comparison')

colors=['#0072b2','#d55e00','#009e73','#cc79a7']
fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for j,ax in enumerate(axs.flat,1):
    ax.plot(qt,qtr[f'R{j}_nm'],'--',color=colors[j-1],label='Quadratic')
    ax.plot(lt,ltr[f'R{j}_nm'],color=colors[j-1],label='Table')
    ax.set(xlabel='Reference-model time t*',ylabel='Equivalent radius (nm)',title=f'Particle {j}')
    ax.legend()
save(fig,'particle_radii_comparison')

fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for j in range(1,5):
    tx=lt[lt<=qt[-1]]
    delta=ltr[f'R{j}_nm'].to_numpy()[:len(tx)]-np.interp(tx,qt,qtr[f'R{j}_nm'])
    axs[0,0].plot(tx,delta,label=f'P{j}',color=colors[j-1])
axs[0,0].set(xlabel='t*',ylabel='R(table) - R(quadratic) (nm)');axs[0,0].legend()
for key,style,title in [('quadratic_control','--','Quadratic'),('lookup','-','Table')]:
    t,_,_,tr,post=runs[key]
    axs[0,1].plot(t,tr.R4_nm**3,style,label=title)
    axs[1,0].plot(t,tr.area_fraction,style,label=title)
    axs[1,1].plot(post.time,post.total_energy/post.total_energy.iloc[0],style,label=title)
axs[0,1].set(xlabel='t*',ylabel='Largest particle R4^3 (nm^3)')
axs[1,0].set(xlabel='t*',ylabel='Gamma-prime area fraction')
axs[1,1].set(xlabel='t*',ylabel='F(t) / F(0)')
for ax in [axs[0,1],axs[1,0],axs[1,1]]: ax.legend()
save(fig,'kinetics_comparison')

if min(qt[-1],lt[-1])>=40000:
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,(lo,hi) in zip(axs,[(5000,30000),(10000,40000)]):
        x=np.arange(4)
        for key,offset,title in [('quadratic_control',-.18,'Quadratic'),('lookup',.18,'Table')]:
            rates=[metrics[key][f'P{j}_dR_dt_{lo}_{hi}']*1e5 for j in range(1,5)]
            ax.bar(x+offset,rates,width=.36,label=title)
        ax.set(xticks=x,xticklabels=['P1','P2','P3','P4'],ylabel='dR/dt* (10^-5 nm)',title=f'Linear fit: {lo:,} <= t* <= {hi:,}')
        ax.axhline(0,color='black',lw=.7);ax.legend()
    save(fig,'growth_rates_comparison')

fig,axs=plt.subplots(2,4,figsize=(12,6),layout='constrained')
end=min(qt[-1],lt[-1])
for row,(t,f,_,_,_) in enumerate([runs['quadratic_control'],runs['lookup']]):
    for col,wanted in enumerate([0,end/3,2*end/3,end]):
        k=int(np.argmin(abs(t-wanted)))
        im=axs[row,col].imshow(f['eta'][k],origin='lower',extent=[0,96,0,96],vmin=0,vmax=1,cmap='magma')
        axs[row,col].set_title(f't*={t[k]:,.0f}')
        axs[row,col].set_xlabel('x (nm)')
        if col==0:axs[row,col].set_ylabel(('Quadratic' if row==0 else 'Table')+'\ny (nm)')
fig.colorbar(im,ax=list(axs.flat),label='eta',shrink=.85)
save(fig,'microstructure_comparison')

if not args.no_gif and not args.partial:
    fig,axs=plt.subplots(1,2,figsize=(8,4.6),layout='constrained')
    ims=[]
    for ax,title in zip(axs,['Quadratic','Table']):
        ims.append(ax.imshow(qf['eta'][0],origin='lower',extent=[0,96,0,96],vmin=0,vmax=1,cmap='magma'))
        ax.set(title=title,xlabel='x (nm)',ylabel='y (nm)')
    st=fig.suptitle('')
    sample=np.linspace(0,60000,101)
    def update(i):
        actual=[]
        for im,t,f in zip(ims,[qt,lt],[qf,lf]):
            k=int(np.argmin(abs(t-sample[i])));im.set_data(f['eta'][k]);actual.append(t[k])
        st.set_text('Actual saved frames | t* = '+', '.join(f'{v:,.0f}' for v in actual))
        return ims+[st]
    FuncAnimation(fig,update,frames=len(sample),interval=100).save(OUT/'comparison.gif',writer=PillowWriter(fps=10))
    plt.close(fig)

def pct(new,old): return 100*(new-old)/old
q,l=metrics['quadratic_control'],metrics['lookup']
if not args.partial:
    rkey='P4_dR_dt_5000_30000'; kkey='P4_dR3_dt_5000_30000'
    eqtable='<table><tr><th>量</th><th>二次近似</th><th>查表</th><th>差值</th></tr>'
    for i,name in enumerate(['γ 平衡 Al / at.%','γ′ 平衡 Al / at.%']):
        a=100*thermo['quadratic_equilibrium_x_Al'][i]; b=100*thermo['common_tangent_x_Al'][i]
        eqtable+=f'<tr><td>{name}</td><td>{a:.9f}</td><td>{b:.9f}</td><td>{b-a:+.9f} 个百分点</td></tr>'
    eqtable+='</table>'
    ratetable='<table><tr><th>指标</th><th>二次近似</th><th>查表</th><th>相对差异</th></tr>'
    rate_rows=[(f'P{j}_dR_dt_5000_30000',f'P{j} 的 dR/dt* / (nm/t*)，5000–30000') for j in range(1,5)]
    rate_rows += [(kkey,'P4 的 d(R³)/dt* / (nm³/t*)，5000–30000'),('P4_dR_dt_10000_40000','P4 的 dR/dt* / (nm/t*)，10000–40000')]
    for key,title in rate_rows:
        ratetable+=f'<tr><td>{title}</td><td>{q[key]:.9g}</td><td>{l[key]:.9g}</td><td>{pct(l[key],q[key]):+.4f}%</td></tr>'
    for j in range(4):
        a=q['radii_final_nm'][j];b=l['radii_final_nm'][j]
        ratetable+=f'<tr><td>P{j+1} 最终半径 / nm</td><td>{a:.6f}</td><td>{b:.6f}</td><td>{b-a:+.6f} nm</td></tr>'
    ratetable+='</table>'
    coretable='<table><tr><th>终点内部平均 Al / at.%</th><th>二次近似</th><th>查表</th></tr>'
    for name,title in [('c_gamma','γ 基体，η&lt;0.01'),('c_gamma_prime','γ′ 颗粒，η&gt;0.99')]:
        coretable+=f'<tr><td>{title}</td><td>{100*q[name+"_final_core_mean"]:.7f}</td><td>{100*l[name+"_final_core_mean"]:.7f}</td></tr>'
    coretable+='</table>'
    body=f'''<h1>Ni–14 at.%Al、1000 K：查表自由能与二次近似的对照</h1>
<p>两版完整 MOOSE 模拟均运行至 t*=60000。查表自由能与原二次近似的平面界面平衡成分差异小于 0.000014 个 at.% 百分点；在 5000–30000 的共同时间窗内，最大颗粒 P4 的半径增长速率相对变化为 <b>{pct(l[rkey],q[rkey]):+.4f}%</b>。这些是当前参考参数二维算例的计算结果，不是实验时效速率。</p>
<h2>对照设计</h2><p>保持 Ni–14 at.%Al、1000 K、γ/γ′ 两相、96×96 nm² 周期区域、96×96 Q1 网格、四个初始颗粒、相同节点 η 与 c 初值、共同摩尔体积 7×10⁻⁶ m³/mol、参考界面能 0.020 J/m²、界面宽度 4 nm、参考 D=10⁻¹⁸ m²/s、守恒迁移率 M=0.002873786184560299 与界面迁移率 L=10 不变。二次模型和查表模型都用局部嵌套 KKS，Newton + MUMPS、隐式 Euler、4 MPI。共同重放原算例接受的时间序列。<b>保持 M 不变</b>：自由能曲率改变时，化学扩散系数随热力学因子改变，不能再调 M 来强行保持每个成分处的 D，否则会同时改变动力学模型。</p>
<p>完整二次控制计算与原解析消元结果比较：c 最大绝对差 {control['c_max_abs_error']:.3g}，η 最大绝对差 {control['eta_max_abs_error']:.3g}，μ 最大绝对差 {control['mu_max_abs_error']:.3g}。两版初始节点 c、η 完全一致：{metrics['identical_initial_c_eta']}；时间序列一致：{metrics['time_sequence_exact']}。辅助相浓度由积分点材料属性输出为单元平均量。</p>
<h2>查表的具体实现</h2><p>使用 Thermo-Calc Educational 2026b / NIDEMO 2.0 已导出的两条单相分支：γ 共 502 点、xAl=0.05–0.30；γ′ 共 242 点、xAl=0.18–0.30。两相减去相同仿射能量参考后，以每个表点的 G 与 dG/dc 构造三次 Hermite 分段插值。预处理器输出每一小段系数，并生成二分区间选择表达式；MOOSE 运行时按当前相浓度选择区间并计算能量。此处是预处理后的表系数查找，不是全区间多项式拟合，也没有在线调用 Thermo-Calc。</p>
<p>自由能的一、二阶导数由同一个分段能量求导，满足 μ=dG/dc。Hermite 插值保证 G 与一阶导连续，二阶导允许在节点处跳变；原始数据有限有效位数仍会影响局部曲率。阻尼局部 Newton 通过表范围条件拒绝越界解。表格节点能量及一阶导重现误差见 thermodynamic_validation.json。当前数据密度未做 Thermo-Calc 重新加密收敛试验，因此不把很小的差异解释为材料本征的高精度预测。</p>
<h2>平衡成分</h2>{eqtable}<p>以上由两条自由能分支的独立公切线求解得到，是平直界面、无弹性条件下的两相平衡成分。原二次模型本来就锚定该平衡成分，所以接近是预期结果。有限尺寸颗粒具有曲率效应，时效终点仍在粗化，不应把颗粒内部浓度当作平面界面的平衡成分。</p>
<p>查表在公切线处的曲率为 γ：{thermo['curvatures_at_tangent_J_mol'][0]:.3f} J/mol，γ′：{thermo['curvatures_at_tangent_J_mol'][1]:.3f} J/mol；二次模型分别为 347973 和 424347.4 J/mol。远离平衡点时，表中的非二次形状改变扩散势和界面驱动力。</p><img src="free_energy_comparison.png">
<h2>颗粒演化与粗化速率</h2><p>每个颗粒的等效半径定义为 R=√(∫h(η)dA/π)，在其初始所在象限积分；使用 Q1 梯形权重。四颗粒在本算例中保持分离，适用该跟踪方式。速率在明确的公共时间窗内，对等间隔时间采样的 R(t) 作线性拟合；同时给出 R³ 的拟合斜率作为有效比较指标。由于是二维、四颗粒、有限区域，<b>R³ 斜率不等同于三维稀释极限 LSW 材料常数</b>。改变拟合窗口也会改变该有效速率。</p>{ratetable}
<img src="particle_radii_comparison.png"><img src="growth_rates_comparison.png"><img src="kinetics_comparison.png">
<p>P1 首次消失的保存帧区间：二次 {q['dissolution_brackets_star']['1']}；查表 {l['dissolution_brackets_star']['1']}。判据为该象限内不再存在 η≥0.5；只能确定区间，不能把保存帧时间当作精确消失时刻。两版最终颗粒数为 {q['final_count']} 与 {l['final_count']}。</p>
<h2>实际组织与质量检查</h2><img src="microstructure_comparison.png"><p>二次/查表的全程平均 Al 分数最大偏差分别为 {q['max_mass_error']:.3g} / {l['max_mass_error']:.3g}。总自由能最大正增量分别为 {q['max_energy_increase']:.3g} / {l['max_energy_increase']:.3g}。自由能总值受共同参考与表数据精度影响，比较演化时展示各自归一化值。输出场均来自 MOOSE Exodus，没有人为删除颗粒或修改轨迹。</p>
<h2>结论边界与复现</h2><p>查表保留已导出分支的非二次形状，但仍假定 γ′ 子晶格在每个成分处局部弛豫，不描述独立的有序化动力学。当前比较控制了空间网格、时间序列及材料参数，尚未完成空间、时间、界面宽度的系统收敛研究。界面能和扩散参数仍为参考值，没有加入相干弹性或三维效应，因此不能转换成可靠的实验小时数。</p>
<p>输入：examples/ni_al_lookup/lookup.i 和 quadratic_control.i；生成器：scripts/prepare_lookup_comparison.py；分析器：scripts/analyze_lookup_comparison.py。输入使用相对项目路径，启动脚本固定项目根目录。完整 Exodus 可直接在 ParaView 打开，坐标按 nm、时间按参考无量纲 t* 解释。本次 MOOSE 为官方预编译 2026.10.14，PETSc 3.25.4，WSL2 Ubuntu 24.04。每个结果文件的 SHA256 与所有统计值保存在 comparison_metrics.json。</p>
<pre>bash scripts/run_lookup_comparison.sh</pre>
<p><a href="lookup.e">查表 Exodus</a> · <a href="quadratic_control.e">二次控制 Exodus</a> · <a href="comparison_metrics.json">数值与校验</a> · <a href="comparison.gif">对照动画</a></p>
<p>实现参考：<a href="https://mooseframework.inl.gov/source/materials/KKSPhaseConcentrationMaterial.html">MOOSE 嵌套 KKS</a>；<a href="https://mooseframework.inl.gov/source/materials/DerivativeParsedMaterial.html">自由能导数材料</a>；<a href="https://mooseframework.inl.gov/source/timesteppers/TimeSequenceStepper.html">指定时间序列</a>。</p>'''
    core_explanation=f'''<p>模拟终点的内部平均成分如下。以四节点 η 平均值选择单元，对单元平均相浓度作等面积平均；边界掩膜变化可能使曲线出现小台阶。</p>{coretable}<img src="core_compositions_comparison.png"><p>全程查表 γ′ 内部成分范围为 {100*l['c_gamma_prime_core_range'][0]:.5f}–{100*l['c_gamma_prime_core_range'][1]:.5f} at.%。实际颗粒内部成分靠近平衡展开点，因此即使远离该范围时自由能形状差异明显，本算例的速率差异也可能很小。公切线数值的额外小数位只用于显示两种表示的数值差别，不代表数据库具有该数量级的物理精度。</p>'''
    body=body.replace('<h2>颗粒演化与粗化速率</h2>','<img src="equilibrium_comparison.png">'+core_explanation+'<h2>颗粒演化与粗化速率</h2>')
    if 'matched_time_control' in metrics:
        mt=metrics['matched_time_control']
        note=f"首先重放原算例的 {mt['original_steps']} 个已接受时间步。查表在颗粒消失附近（首次额外时刻 t*={mt['first_additional_time']:.6f}）触发自动减步，最终接受 {mt['lookup_steps']} 步。为消除这项差异，二次模型从 t*={mt['checkpoint_time_star']:.6f} 的检查点恢复，严格采用查表版的完整已接受时间序列。下文比较使用这个匹配版本；0–40000 的速率拟合窗口完全不受补算影响。"
        body=body.replace('共同重放原算例接受的时间序列。',note)
        body=body.replace('的检查点恢复，严格采用','的保存场通过 SolutionIC 重新启动（隐式 Euler），指定采用')
        body=body.replace('二次模型和查表模型都用局部嵌套 KKS，','先用局部嵌套 KKS 对二次模型作整段数值控制；最终匹配时间步的二次终段采用已验证等价的解析消元，以避免局部 Newton 在小颗粒消失时的额外减步。查表模型采用局部嵌套 KKS。共同使用 ')
        body=body.replace('quadratic_control.e">二次控制 Exodus','quadratic_matched.e">二次控制 Exodus')
        body=body.replace('<h2>结论边界与复现</h2>','<p>额外时间细分对二次模型终点字段的最大绝对影响：'+html.escape(json.dumps(metrics['quadratic_time_refinement_final_max_field_change']))+'。原时间序列的完整二次控制结果另存 quadratic_control.e，用于与原算例核对。</p><h2>结论边界与复现</h2>')
        body=body.replace('<h2>结论边界与复现</h2>','<p>匹配控制的前段和新算出的终段保存为独立来源，随后逐变量原样拼接为完整 Exodus 序列；未插值或修改场值。初始保存场重载误差与拼接记录见 assembled_control.json。直接使用检查点恢复会带回旧时间步控制器，因此最终流程采用保存场重新启动。</p><h2>结论边界与复现</h2>')
        assembly=metrics.get('assembled_control',{})
        extra=assembly.get('additional_control_times',[])
        if extra:
            warning=f'<p><b>终段仍存在时间步差异：</b>二次解析控制自动增加 {len(extra)} 个时刻：{extra}。因此这是一项终段时间细分敏感性检查，未实现终段时间网格完全相同。两个主要速率窗口（至 t*=40000）的时间步逐项一致，严格对照结论限于这些窗口。消失区间及终点半径只能作为带有时间步影响的辅助结果。二次模型终点半径相对原时间序列的变化 / nm：{metrics["quadratic_time_refinement_final_radius_change_nm"]}。</p>'
            body=body.replace('<h2>颗粒演化与粗化速率</h2>',warning+'<h2>颗粒演化与粗化速率</h2>')
    qa=q['dissolution_brackets_star']['1']; la=l['dissolution_brackets_star']['1']
    if qa and la and max(qa[0],la[0])<min(qa[1],la[1])-1e-7:
        body=body.replace('<h2>实际组织与质量检查</h2>','<p><b>两个消失时间区间重叠。</b>现有保存帧不能分辨哪个模型的颗粒确定先消失，不能把首次无颗粒帧的时间差当作精确的溶解时间差。</p><h2>实际组织与质量检查</h2>')
    elif qa and la and la[1]<=qa[0]+1e-7:
        body=body.replace('<h2>实际组织与质量检查</h2>',f'<p>补充控制在共同保存时刻 t*={la[1]:.6f} 提供了进一步证据：查表版已无 η≥0.5 的 P1 区域，而二次版仍存在。以这一判据看查表版先消失，但时间差只能约束在 0 到 {qa[1]-la[0]:.3f} 个 t* 单位之间，不能用两次首次无颗粒帧之差作为精确时间差。区间按（最后存在时刻，首次不存在时刻] 理解。</p><h2>实际组织与质量检查</h2>')
    doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>查表与二次近似对照报告</title><style>body{max-width:1080px;margin:40px auto;padding:0 28px;font-family:system-ui,"Microsoft YaHei",sans-serif;line-height:1.8;color:#182c3b}h1{font-size:28px}h2{margin-top:38px;color:#125c78}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #cbd5df;padding:9px;text-align:left}th{background:#eaf2f6}img{max-width:100%;display:block;margin:24px auto}pre{background:#edf2f5;padding:16px;overflow:auto}a{color:#006e9b}@media print{body{margin:0}img,table{break-inside:avoid}}</style><body>'+body+'</body></html>'
    # Keep the report portable: figures are embedded, while Exodus and CSV stay separate.
    def embedded_image(match):
        encoded=base64.b64encode((OUT/match.group(1)).read_bytes()).decode('ascii')
        return '<img src="data:image/png;base64,'+encoded+'">'
    doc=re.sub(r'<img src="([^"/]+\.png)">',embedded_image,doc)
    (OUT/'REPORT.html').write_text(doc,encoding='utf-8')
print(json.dumps(metrics,indent=2))
