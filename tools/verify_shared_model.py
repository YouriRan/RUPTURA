#!/usr/bin/env python3
"""Reproduce input/smoke checks; --references also runs Marx and Haghpanah cases.

Use --executable /path/to/ruptura. Results go to verification/integration;
complete rerun outputs are in verification/integration/runs and can be large.
"""
from pathlib import Path
import argparse,concurrent.futures,json,os,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,default=ROOT/'build/ruptura');p.add_argument('--smoke',action='store_true');p.add_argument('--references',action='store_true');p.add_argument('--jobs',type=int,default=2);args=p.parse_args();exe=args.executable.resolve();report=ROOT/'verification/integration';report.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 inputs=sorted((ROOT/'examples/new_examples/integrated').rglob('simulation.json'))
 def validate(path):
  try:r=subprocess.run([str(exe),str(path),'--validate'],capture_output=True,text=True,timeout=60);return {'input':str(path.relative_to(ROOT)),'code':r.returncode,'error':r.stderr[-2000:]}
  except subprocess.TimeoutExpired:return {'input':str(path.relative_to(ROOT)),'code':'TIMEOUT'}
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:rows=list(pool.map(validate,inputs))
 (report/'input_validation.json').write_text(json.dumps(rows,indent=2)+'\n');failed=[r for r in rows if r['code']!=0];print('Input checks:',len(rows),'passed',len(rows)-len(failed),flush=True)
 for row in failed:print(row,flush=True)
 def run(item):
  label,source,changes,limit=item;directory=report/'runs'/label;directory.mkdir(parents=True,exist_ok=True);doc=json.loads(source.read_text());doc.update(changes);path=directory/'simulation.json';path.write_text(json.dumps(doc,indent=2)+'\n');t=time.monotonic()
  try:
   r=subprocess.run([str(exe),str(path),'--no-plots'],capture_output=True,text=True,timeout=limit,env=env);(directory/'run.log').write_text(r.stdout+r.stderr);row={'case':label,'source':str(source.relative_to(ROOT)),'overrides':changes,'code':r.returncode,'seconds':round(time.monotonic()-t,3),'error':r.stderr[-2500:]}
   if r.returncode==0 and (directory/'results/summary.json').is_file():
    s=json.loads((directory/'results/summary.json').read_text());m=s['final_metrics'];row.update(cycles=s['cycles'],converged=s['converged'],mass_balance_error=m['mass_balance_error'],energy_balance_error=m['energy_balance_error'],final_metrics=m)
  except subprocess.TimeoutExpired:row={'case':label,'code':'TIMEOUT','seconds':limit}
  print(label,row['code'],row.get('seconds'),row.get('error','')[-250:],flush=True);return row
 if args.smoke:
  cases=[(f'{kind}_{fluid}_{geometry}',(ROOT/f'examples/new_examples/unified/fvm/{fluid}_{"breakthrough" if kind=="breakthrough" else "cycle"}/simulation.json' if geometry=='fixed_bed' else ROOT/'examples/new_examples/integrated'/kind/f'{fluid}_{geometry}/simulation.json'),{},240) for kind in ['breakthrough','cycles'] for fluid in ['gas','liquid'] for geometry in ['fixed_bed','monolith']]
  cases += [(f'particle_{fluid}_{shape}',ROOT/f'examples/new_examples/integrated/particle/{fluid}_hsdm_{shape}/simulation.json',{},240) for fluid in ['gas','liquid'] for shape in ['sphere','cylinderradial']]
  cases += [(f'thermal_{fluid}',ROOT/f'examples/new_examples/integrated/thermal/{fluid}_transient_coolant/simulation.json',{},240) for fluid in ['gas','liquid']]
  cases += [('automatic_saturation',ROOT/'examples/new_examples/integrated/runtime/automatic_saturation/simulation.json',{},300)]
  cases += [('native_bdf',ROOT/'examples/new_examples/unified/fvm/gas_breakthrough/simulation.json',{'BreakthroughIntegrator':'BDF'},240)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:rows=list(pool.map(run,cases))
  (report/'smoke_runs.json').write_text(json.dumps(rows,indent=2)+'\n')
 if args.references:
  cases=[('marx_T15_original',ROOT/'examples/new_examples/integrated/validation/marx_T15_original/simulation.json',{},360),('marx_T15_fvm',ROOT/'examples/new_examples/integrated/validation/marx_T15_fvm/simulation.json',{},360),('haghpanah_2013_local_pressure',ROOT/'examples/new_examples/integrated/validation/haghpanah_2013_equations/dcb_local_pressure/simulation.json',{},360),('haghpanah_2013_reservoir',ROOT/'examples/new_examples/integrated/validation/haghpanah_2013_equations/dcb_local_pressure/simulation.json',{'InletDensityBasis':'Reservoir'},360),('haghpanah_2014_adsorption',ROOT/'examples/new_examples/integrated/validation/haghpanah_2014/01_breakthrough/adsorption/simulation.json',{'BreakthroughIntegrator':'CVODE','CVODELinearSolver':'SPGMR'},360)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:rows=list(pool.map(run,cases))
  (report/'reference_runs.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
