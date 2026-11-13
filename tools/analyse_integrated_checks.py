#!/usr/bin/env python3
"""Summarize saved reference/particle checks; never classify them as validated fits."""
from pathlib import Path
import csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'verification/integration'
RUNS=REPORT/'runs'

def recorded(name):
    with np.load(RUNS/name/'results/states.npz',allow_pickle=False) as data:
        t=data['time'].copy();c=data['A_c'][:,:,-1];y=np.divide(c,c.sum(axis=1,keepdims=True),out=np.zeros_like(c),where=c.sum(axis=1,keepdims=True)>0)
        return t,y,data['A_Tg'][:,-1].copy()

def crossing(t,y,level,after=0.,up=True):
    for i in range(1,len(t)):
        if t[i]<after:continue
        if (up and y[i-1]<=level<=y[i]) or (not up and y[i-1]>=level>=y[i]):
            return float(t[i-1]+(t[i]-t[i-1])*(level-y[i-1])/(y[i]-y[i-1])) if y[i]!=y[i-1] else float(t[i])
    return None

def text_outlet(path,length=1.2):
    rows=[]
    with path.open() as stream:
        for line in stream:
            if not line.strip() or line.startswith('#'):continue
            tokens=line.split()
            if math.isclose(float(tokens[2]),length,abs_tol=1e-8):rows.append([float(tokens[1])*60.,float(tokens[5])])
    return np.asarray(rows)

def main():
    summary={};plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for column,(name,label) in enumerate([('haghpanah_2013_local_pressure','2013 source-derived case · local-pressure inlet'),('haghpanah_2014_adsorption','2014 conditional reconstruction')]):
        if not (RUNS/name/'results/summary.json').exists():continue
        t,y,temp=recorded(name);axes[0,column].plot(t/60.,y[:,0],label='Shared FVM');axes[1,column].plot(t/60.,temp)
        summary[name]={'duration_s':float(t[-1]),'maximum_outlet_temperature_K':float(max(temp)),
            'final_CO2_mole_fraction':float(y[-1,0]),'first_front_half_plateau_s':crossing(t,y[:,0],.1275/2.),
            'cooling_front_midpoint_s':crossing(t,temp,(max(temp)+298.15)/2.,after=2000.,up=False)}
        axes[0,column].set_title(label);axes[0,column].set_ylabel('Outlet CO₂ mole fraction');axes[1,column].set_ylabel('Outlet fluid temperature [K]');axes[1,column].set_xlabel('Time [min]')
    reservoir=RUNS/'haghpanah_2013_reservoir/results/summary.json'
    if reservoir.exists():
        t,y,temp=recorded('haghpanah_2013_reservoir');axes[0,0].plot(t/60.,y[:,0],ls='--',label='Matched reservoir-density inlet');axes[1,0].plot(t/60.,temp,ls='--')
        summary['haghpanah_2013_reservoir']={'duration_s':float(t[-1]),'maximum_outlet_temperature_K':float(max(temp)),
            'first_front_half_plateau_s':crossing(t,y[:,0],.1275/2.),'cooling_front_midpoint_s':crossing(t,temp,(max(temp)+298.15)/2.,after=2000.,up=False)}
    markers=ROOT/'examples/new_examples/integrated/validation/haghpanah_2014/reference/experimental_markers.csv'
    if markers.exists() and (RUNS/'haghpanah_2014_adsorption/results/summary.json').exists():
        rows=[r for r in csv.DictReader(markers.open()) if r['experiment']=='breakthrough_adsorption' and r['variable']=='outlet_CO2']
        x=np.array([float(r['phase_time_s']) for r in rows]);v=np.array([float(r['value']) for r in rows]);axes[0,1].scatter(x/60.,v,s=25,facecolors='none',edgecolors='#bd5545',label='Approximate donor digitization')
        t,y,temp=recorded('haghpanah_2014_adsorption');error=np.interp(x,t,y[:,0])-v
        summary['haghpanah_2014_adsorption'].update(digitized_marker_RMSE_mole_fraction=float(np.sqrt(np.mean(error**2))),marker_origin='Approximate manual digitization in supplied donor, 50 s / 0.002 stated resolution; no acceptance threshold or refit')
    for ax in axes.flat:ax.grid(alpha=.2)
    for ax in axes[0]:ax.legend(fontsize=8)
    fig.suptitle('Haghpanah breakthrough comparisons · literature validation remains open',fontsize=13)
    fig.savefig(REPORT/'haghpanah_comparison.png',dpi=180);plt.close(fig)

    fvm=RUNS/'marx_T15_fvm/results/summary.json';original=RUNS/'marx_T15_original'
    if fvm.exists() and (original/'component_1_CO2.data').exists():
        t,y,temp=recorded('marx_T15_fvm');fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for j,(name,fraction) in enumerate([('CO2',.4),('N2',.3),('H2',.3)],start=1):
            old=text_outlet(original/f'component_{j}_{name}.data');axes[0].plot(old[:,0]/60.,old[:,1]/fraction,label=name+' original');axes[0].plot(t/60.,y[:,j]/fraction,ls='--',label=name+' FVM')
        oldtemp=text_outlet(original/'column.data');axes[1].plot(oldtemp[:,0]/60.,oldtemp[:,1],label='Original main');axes[1].plot(t/60.,temp,ls='--',label='Shared FVM')
        axes[0].set(ylabel='Outlet y/yfeed',xlabel='Time [min]');axes[1].set(ylabel='Outlet fluid temperature [K]',xlabel='Time [min]')
        for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle('Marx T15 · model comparison, without a fit or equivalence claim');fig.savefig(REPORT/'marx_T15_comparison.png',dpi=180);plt.close(fig)
        old_co2=text_outlet(original/'component_1_CO2.data')
        summary['marx_T15']={'horizon_s':float(t[-1]),'FVM_maximum_outlet_temperature_K':float(max(temp)),
            'original_maximum_outlet_temperature_K':float(max(oldtemp[:,1])),
            'FVM_CO2_half_feed_s':crossing(t,y[:,1],.2),'original_CO2_half_feed_s':crossing(old_co2[:,0],old_co2[:,1],.2),
            'status':'Full simulations completed; equations and thermal inputs differ; this is a comparison, not matched-model validation.'}

    particles={};fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for fluid in ['gas','liquid']:
        for shape,power in [('sphere',2),('cylinderradial',1)]:
            name=f'particle_{fluid}_{shape}';directory=RUNS/name/'results'
            if not (directory/'summary.json').exists():continue
            meta=json.loads((directory/'metadata.json').read_text());record=meta['column_records']['columns']['A'];nz=len(record['z']);nc=len(meta['component_names']);groups=record['history_group_names'];nr=sum('/HSDM/' in x for x in groups);nload=len(record['loading_group_names'])
            with np.load(directory/'states.npz',allow_pickle=False) as data:
                states=data['state'];offset=(1+nload)*nc*nz;shells=np.stack([states[:,offset+r*nc*nz:offset+(r+1)*nc*nz].reshape(-1,nc,nz) for r in range(nr)],axis=1)
                edges=np.linspace(0.,1.,nr+1);weights=np.diff(edges**(power+1));means=np.einsum('r,trcn->tcn',weights,shells);physical=states[:,nc*nz:2*nc*nz].reshape(-1,nc,nz);err=float(np.max(np.abs(physical-means)));assert err<2e-6,(name,err)
                particles[name]={'maximum_mean_shell_loading_difference_mol_kg':err,'shells':nr,'completed':True}
                ax.plot((edges[:-1]+edges[1:])/2.,shells[-1,:,0,0],label=f'{fluid} · {shape}')
    if particles:
        ax.set(xlabel='Relative radius r/R',ylabel='A loading at inlet cell [mol/kg]',title='Finite-film HSDM · final radial loading');ax.legend();ax.grid(alpha=.2);fig.savefig(REPORT/'particle_profiles.png',dpi=180)
    plt.close(fig);summary['particle_conservation']=particles
    (REPORT/'comparison_metrics.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
