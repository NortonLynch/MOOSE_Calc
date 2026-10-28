"""Diagnostic only: checkpoint recovery restores the old time sequence.

This attempted approach is retained for audit. Use prepare_lookup_tail.py and
assemble_lookup_control.py for the validated matched-time comparison instead.
"""
from pathlib import Path
import re, shutil, json, argparse
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--output',default='results/ni_al_lookup')
args=parser.parse_args()
OUT=ROOT/args.output
lookup=pd.read_csv(OUT/'lookup.csv')
quad=pd.read_csv(OUT/'quadratic_control.csv')
assert lookup.time.iloc[-1]==quad.time.iloc[-1]==60000
prefix=quad.iloc[:1001]
assert np.allclose(prefix.time,lookup.time.iloc[:1001],atol=1e-7,rtol=0)
target=OUT/'quadratic_matched.e'
assert not target.exists(), 'Preserve existing matched control'
shutil.copy2(OUT/'quadratic_control.e',target)
prefix.to_csv(OUT/'quadratic_matched.csv',index=False,float_format='%.17g')
text=(ROOT/'examples/ni_al_lookup/quadratic_control.i').read_text()
times=' '.join(format(t,'.17g') for t in lookup.time)
text=re.sub(r"time_sequence = '[^']*'",f"time_sequence = '{times}'\n    reset_dt = true",text)
text=text.replace('results/ni_al_lookup/quadratic_control',OUT.relative_to(ROOT).as_posix()+'/quadratic_matched')
(ROOT/'examples/ni_al_lookup/quadratic_matched.i').write_text(text)
lookup[['time']].to_csv(ROOT/'examples/ni_al_lookup/lookup_accepted_time_sequence.csv',index=False)
extra=[float(t) for t in lookup.time if np.min(abs(quad.time.to_numpy()-t))>1e-7]
meta={'restart_checkpoint':'quadratic_control_cp/1000','checkpoint_time_star':float(prefix.time.iloc[-1]),'lookup_steps':len(lookup)-1,'original_steps':len(quad)-1,'extra_lookup_times':extra,'first_additional_time':min(extra) if extra else None}
(OUT/'matched_time_control.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
