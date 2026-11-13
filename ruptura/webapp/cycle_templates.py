"""Complete, editable demonstration cycles using the public main JSON schema."""
from copy import deepcopy
import json
from pathlib import Path

def templates():
    base=json.loads((Path(__file__).with_name('cycle_template_base.json')).read_text())
    # Installed wheels may not contain the source example tree.
    return {name: build(base,n) for name,n in [('single-psa',1),('two-psa-pe',2),('four-psa',4)]}

def build(base,n):
    d=deepcopy(base);names=list('ABCD'[:n]);d['Columns']={k:{} for k in names}
    d['DisplayName']=f'{n}-column PSA'+(' with pressure equalization' if n>1 else '')
    d['SpatialDiscretization']='FVM';d['NumberOfCycles']=2
    d['CyclicSteadyState']={'Enabled':True,'MaximumCycles':100,'StateRelativeTolerance':1e-4,'StateAbsoluteTolerance':1e-7,'InventoryRelativeTolerance':1e-5,'MetricRelativeTolerance':1e-4,'ConsecutiveCycles':3}
    d['NumberOfGridPoints']=12;d['Output']['SamplesPerPhase']=100
    def valve(name,src,dst,tag='',conductance=1e-8):
        return {'Name':name,'From':src,'To':dst,'FlowModel':'linear','Conductance':conductance,'MaximumMolarFlow':.05,'CheckValve':True,'Opening':1.,'Tag':tag or 'equalization'}
    def phase(name,duration,roles=None,connections=None):
        links=list(connections or [])
        for col,role in (roles or {}).items():
            if role in ('PRESS','ADS'):links.append(valve('V-F'+col,'source:feed',col+'.left','feed' if role=='ADS' else 'repressurization_gas'))
            if role=='ADS':links.append(valve('V-L'+col,col+'.right','source:light_receiver','light_product'))
            if role in ('BD','PURGE'):links.append(valve('V-H'+col,col+'.left','source:heavy_receiver','heavy_product'))
            if role=='PURGE':links.append(valve('V-P'+col,'source:external_purge',col+'.right','purge'))
        return {'Name':name,'Duration':duration,'Columns':{k:{'Name':(roles or {}).get(k,'WAIT'),'LeftPort':{'Mode':'closed'},'RightPort':{'Mode':'closed'}} for k in names},'Connections':links}
    if n==1:
        phases=[phase(r,t,{'A':r}) for r,t in [('PRESS',30),('ADS',60),('BD',30),('PURGE',60)]]
    elif n==2:
        phases=[phase('Startup / repressurize A',30,{'A':'PRESS','B':'BD'}),phase('A adsorption / B purge',60,{'A':'ADS','B':'PURGE'}),phase('A to B equalization',5,connections=[valve('V-E','A.right','B.right')]),phase('A blowdown / B pressurization',30,{'A':'BD','B':'PRESS'}),phase('B adsorption / A purge',60,{'B':'ADS','A':'PURGE'}),phase('B to A equalization',5,connections=[valve('V-E','B.right','A.right')])]
    else:
        phases=[]
        for i in range(4):
            roles={names[(i+j)%4]:r for j,r in enumerate(['ADS','BD','PURGE','PRESS'])}
            phases.append(phase(f'Rotation {i+1}',30,roles))
            phases.append(phase(f'Equalization {i+1}',5,connections=[valve('V-E'+str(i+1),names[i]+'.right',names[(i+3)%4]+'.right')]))
    d['SwingAdsorptionPhases']=phases
    d.setdefault('Performance',{})['FreshFeedTags']=['feed','repressurization_gas']
    d['Notes']=['Editable demonstration recipe. Timings and valve conductances are assumptions, not optimized designs or literature validation.','CSS must be demonstrated by the calculated state, inventory and metric tests.']
    return d
