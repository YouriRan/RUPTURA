"""Plot the shared model's recorded arrays, including compact cycle outputs."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def figure(directory, config, options):
    output=Path(config.get('Output',{}).get('Directory','results'))
    if not output.is_absolute():output=Path(directory)/output
    metadata=json.loads((output/'metadata.json').read_text())
    names=metadata['column_names'];name=options.get('column',names[0])
    if name not in names:raise ValueError('Choose a column from the saved run.')
    with np.load(output/'states.npz',allow_pickle=False) as archive:
        data={key:archive[key] for key in archive.files}
    t=np.asarray(data['time']);c=data[name+'_c'];q=data[name+'_q'];nz=c.shape[-1]
    record=metadata['column_records']['columns'][name]
    z=np.asarray(record['z']);mode=options.get('plot','breakthrough');metric=options.get('metric','C');metric={'C':'concentration','Q':'loading','Qtotal':'loading','Y':'mole_fraction','P':'partial_pressure','Pt':'pressure'}.get(metric,metric)
    components=config['Components'];selected=options.get('components')
    if selected is None:selected=[j for j,x in enumerate(components) if not x.get('CarrierGas') or options.get('carrier')]
    selected=[j for j in selected if 0<=j<c.shape[1]]
    factor={'s':1.,'min':60.,'hr':3600.}.get(options.get('xUnits','min'),60.)
    units=options.get('xUnits','min');liquid=config.get('FluidPhase','Gas')=='Liquid'
    feed=np.array([x.get('LiquidPhaseConcentration',0.) if liquid else x.get('GasPhaseMolFraction',0.) for x in components],float)
    if not liquid and feed.sum()>0:feed/=feed.sum()
    total=c.sum(axis=1,keepdims=True)
    fractions=np.divide(c,total,out=np.zeros_like(c),where=total>0)
    temperature=data[name+'_Tg'];pressure=data[name+'_P'];frame=int(options.get('frame',len(t)-1));frame=max(0,min(frame,len(t)-1))
    grid=int(options.get('grid',-1));grid=grid if grid>=0 else nz+grid
    if not 0<=grid<nz:raise ValueError('Choose a recorded axial position.')
    profile=mode=='profile' and options.get('xAxis','position')=='position'
    x=z if profile else t/factor
    def values(array):return array[frame] if profile else array[:,grid]
    if mode=='particle':
        groups=record.get('history_group_names',[]);radial=[(i,label) for i,label in enumerate(groups) if '/HSDM/' in label]
        if not radial:raise ValueError('This column has no HSDM radial states.')
        offset=metadata['column_slices'][name][0];nloading=len(record['loading_group_names']);nc=len(components)
        fig=go.Figure();materials={label.split('/HSDM/')[0] for i,label in radial}
        for material in sorted(materials):
            cells=[i for i,label in radial if label.startswith(material+'/HSDM/')]
            for j in selected:
                values=[data['state'][frame,offset+(1+nloading+i)*nc*nz+j*nz+grid] for i in cells]
                fig.add_trace(go.Scatter(x=(np.arange(len(cells))+.5)/len(cells),y=values,name=f'{material} · {components[j]["Name"]}',mode='lines+markers'))
        fig.update_layout(xaxis_title='Relative particle radius r/R',yaxis_title='Loading [mol/kg]',template='plotly_dark' if options.get('theme')=='dark' else 'plotly_white')
        return fig
    traces=[];ylabel='Concentration [mol/m³]' 
    if mode=='temperature' or metric in ('gas_temperature','solid_temperature','wall_temperature','temperature','Tg','Ts','Tw','Tj'):
        fields=['Tg','Ts','Tw'] if mode=='temperature' else [{'gas_temperature':'Tg','solid_temperature':'Ts','wall_temperature':'Tw','temperature':'Tg'}.get(metric,metric)]
        for field in fields:
            if name+'_'+field in data:traces.append((field,values(data[name+'_'+field])))
        ylabel='Temperature [K]'
        if options.get('temperatureUnits')=='celsius':traces=[(label,np.asarray(y)-273.15) for label,y in traces];ylabel='Temperature [°C]'
    elif metric in ('pressure','total_pressure'):
        traces=[('Pressure',values(pressure))];ylabel='Absolute pressure [Pa]'
    else:
        for j in selected:
            array=c[:,j,:]
            if metric in ('physisorption','loading','total_loading'):
                array=q[:,j,:];ylabel='Loading [mol/kg]'
            elif metric in ('mole_fraction','gas_mole_fraction'):
                array=fractions[:,j,:];ylabel='Tracked-solute fraction [-]' if liquid else 'Mole fraction [-]'
            elif metric in ('partial_pressure',):
                array=fractions[:,j,:]*pressure;ylabel='Partial pressure [Pa]'
            elif mode=='breakthrough':
                ymode=options.get('yUnits','normalized_concentration')
                if ymode in ('normalized','normalised','ratio','normalized_concentration'):
                    numerator=c[:,j,:] if liquid else fractions[:,j,:]
                    if feed[j]<=0:continue
                    array=numerator/feed[j];ylabel='C/Cfeed [-]' if liquid else 'y/yfeed [-]'
                elif ymode in ('molefraction','mole_fraction','fraction'):
                    array=fractions[:,j,:];ylabel='Mole fraction [-]'
            traces.append((components[j]['Name'],values(array)))
    stacked=False
    fig=make_subplots(rows=2 if stacked else 1,cols=1,shared_xaxes=stacked,vertical_spacing=.1)
    for label,y in traces:fig.add_trace(go.Scatter(x=x,y=y,name=label,mode='lines+markers' if options.get('markers') else 'lines'),row=1,col=1)
    fig.update_yaxes(title_text=ylabel,row=1,col=1)
    if stacked:
        for field in ('Tg','Ts','Tw'):
            fig.add_trace(go.Scatter(x=t/factor,y=data[name+'_'+field][:,grid],name=field),row=2,col=1)
        fig.update_yaxes(title_text='Temperature [K]',row=2,col=1)
    xlabel='Position [m]' if profile else f'Time [{units}]'
    fig.update_xaxes(title_text=xlabel,row=2 if stacked else 1,col=1)
    if not profile and config.get('SimulationType')=='SwingAdsorption':
        path=output/'phase_metrics.csv'
        if path.is_file():
            rows=list(csv.DictReader(path.open()))
            for row in rows[-100:]:
                boundary=float(row['start_time_s'])+float(row['duration_s'])
                fig.add_vline(x=boundary/factor,line_dash='dot',line_color='#9da8a4')
    fig.update_layout(title=f'{name} · {"axial profile" if profile else "recorded history"}',template='plotly_dark' if options.get('theme')=='dark' else 'plotly_white')
    return fig


def metadata(directory,config):
    output=Path(config.get('Output',{}).get('Directory','results'))
    if not output.is_absolute():output=Path(directory)/output
    record=json.loads((output/'metadata.json').read_text())
    with np.load(output/'states.npz',allow_pickle=False) as data:times=data['time'].copy()
    names=record['column_names'];positions={name:record['column_records']['columns'][name]['z'] for name in names}
    components=[{'index':i,'name':c['Name'],'carrier':bool(c.get('CarrierGas')),'adsorbing':bool(c.get('PhysisorptionSites') or c.get('ChemisorptionSites'))} for i,c in enumerate(config['Components'])]
    choices=[{'id':'breakthrough','label':'Breakthrough curves'},{'id':'temperature','label':'Temperature history'},{'id':'profile','label':'Column profiles'}]
    if any('/HSDM/' in group for name in names for group in record['column_records']['columns'][name].get('history_group_names',[])):choices.append({'id':'particle','label':'Particle loading profiles'})
    fields=[('C','Concentration [mol/m³]','component'),('Q','Loading [mol/kg]','component'),('Y','Tracked component fraction [-]','component'),('P','Partial pressure [Pa]','component'),('Pt','Absolute pressure [Pa]','column'),('Tg','Fluid temperature [K]','column'),('Ts','Solid temperature [K]','column'),('Tw','Wall temperature [K]','column')]
    if config.get('FluidPhase')=='Liquid':fields=[x for x in fields if x[0] not in ('Y','P')]
    return {'kind':config.get('SimulationType','Breakthrough'),'components':components,'columns':names,'columnPositions':positions,'positions':positions[names[0]],'energyBalance':bool(config.get('energyBalance')),'liquid':config.get('FluidPhase')=='Liquid','plots':choices,'metrics':[{'id':key,'label':label,'source':source} for key,label,source in fields],'frames':len(times),'totalFrames':len(times),'sampled':False,'times':(times/60.).tolist()}
