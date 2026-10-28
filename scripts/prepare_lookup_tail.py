"""Fresh backward-Euler restart via saved nodal fields, so no old timestepper state is restored."""
from pathlib import Path
import argparse,re,json
import netCDF4,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',default='results/ni_al_lookup');a=p.parse_args()
OUT=ROOT/a.output;rel=OUT.relative_to(ROOT).as_posix()
q=pd.read_csv(OUT/'quadratic_control.csv'); t0=float(q.time.iloc[1000])
l=pd.read_csv(OUT/'lookup.csv'); times=l.time[l.time>=t0-1e-7].to_numpy().copy();times[0]=t0
with netCDF4.Dataset(OUT/'quadratic_control.e') as d:
    index=int(np.argmin(abs(d['time_whole'][:]-t0))); assert abs(d['time_whole'][index]-t0)<1e-7
text=(ROOT/'examples/ni_al_coarsening/coarsening.i').read_text()
v='[Variables]\n'
for name in ['eta','c','mu']:
    v+=f'  [{name}]\n    [InitialCondition]\n      type = SolutionIC\n      solution_uo = saved\n      from_variable = {name}\n    []\n  []\n'
v+='[]\n'
text=re.sub(r'\[Variables\].*?(?=\[Functions\])',lambda m:v,text,flags=re.S)
text=re.sub(r'\[Functions\].*?(?=\[BCs\])','',text,flags=re.S)
text+=f'\n[UserObjects]\n  [saved]\n    type = SolutionUserObject\n    mesh = {rel}/quadratic_control.e\n    timestep = {index+1}\n    system_variables = \'eta c mu\'\n  []\n[]\n'
text=text.replace('  type = Transient','  type = Transient\n  start_time = '+format(t0,'.17g'))
text=re.sub(r'  \[TimeStepper\].*?  \[\]',"  [TimeStepper]\n    type = TimeSequenceStepper\n    time_sequence = '"+' '.join(format(t,'.17g') for t in times)+"'\n  []",text,flags=re.S)
text=text.replace('results/ni_al_coarsening/coarsening',rel+'/quadratic_analytic_tail')
# The replacement above must not change the source-solution filename.
text=text.replace('mesh = '+rel+'/quadratic_matched_tail.e','mesh = '+rel+'/quadratic_control.e')
text=text.replace('mesh = '+rel+'/quadratic_control.e','mesh = ../../'+rel+'/quadratic_control.e')
for var in ['c_gamma','c_gamma_prime']:
    text=text.replace(f'  [{var}]\n  []',f'  [{var}]\n    family = MONOMIAL\n    order = CONSTANT\n  []')
for aux,var,prop in [('matrix_composition','c_gamma','cm'),('precipitate_composition','c_gamma_prime','cd')]:
    text=re.sub(r'  \['+aux+r'\].*?  \[\]',f'  [{aux}]\n    type = MaterialRealAux\n    variable = {var}\n    property = {prop}\n    execute_on = \'INITIAL TIMESTEP_END\'\n  []',text,flags=re.S)
extra_mat='[Materials]\n'
for prop,ceq,A in [('cm',.124644,347.9729999999997),('cd',.237162,424.3473999999995)]:
    expr=f'h:=eta^3*(6*eta^2-15*eta+10);q:=(1-h)/347.9729999999997+h/424.3473999999995;{ceq}+(c-0.124644-h*0.112518)/q/{A}'
    extra_mat+=f'  [{prop}_analytic]\n    type = ParsedMaterial\n    property_name = {prop}\n    coupled_variables = \'c eta\'\n    expression = \'{expr}\'\n  []\n'
text=text.replace('[Materials]',extra_mat,1)
(ROOT/'examples/ni_al_lookup/quadratic_matched_tail.i').write_text(text)
extra=[float(t) for t in l.time if min(abs(q.time.to_numpy()-t))>1e-7]
(OUT/'matched_time_control.json').write_text(json.dumps({'restart_method':'SolutionIC fresh implicit-Euler restart; no recovered timestepper','checkpoint_time_star':t0,'source_frame_index_zero_based':index,'lookup_steps':len(l)-1,'original_steps':len(q)-1,'extra_lookup_times':extra,'first_additional_time':min(extra) if extra else None},indent=2)+'\n')
print('Start',t0,'Exodus frame (1-based)',index+1,'tail intervals',len(times)-1)
