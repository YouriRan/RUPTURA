#!/usr/bin/env python3
"""Generate matched public gas/liquid chemistry examples; no literature data."""
from copy import deepcopy
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ruptura.webapp.config import KINETICS


def main():
    output=ROOT/'examples/new_examples/shared_chemisorption';manifest=[]
    for method in ['FDM','FVM']:
        for liquid in [False,True]:
            phase='liquid' if liquid else 'gas'
            base=json.loads((ROOT/f'ruptura/webapp/presets/fvm-{phase}-breakthrough.json').read_text())
            base.update(SpatialDiscretization=method,NumberOfGridPoints=8,SimulationTime=10.,WriteInterval=.1,
                        PhysisorptionEnabled=False,BreakthroughIntegrator='CVODE',ConvectionScheme='WENO3')
            if liquid:
                base["MixturePredictionMethod"]="SPI"
                for component in base['Components']:
                    component.pop('CarrierGas',None)
                    if not component.get('PhysisorptionSites'):component['PhysisorptionSites']=[{'Type':'Henry','Parameters':[0.]}]
            base['Solver']={'RelativeTolerance':1e-6,'AbsoluteTolerance':1e-9,'MaximumTimeStep':.1}
            base['CVODEAbsoluteToleranceLoading']=1e-10
            base['CVODEAbsoluteToleranceConcentration']=1e-10
            base['Output']={'Directory':'results','Plots':True,'MainFormatFiles':True,'Checkpoint':True}
            for model,fields in KINETICS.items():
                d=deepcopy(base);params={key:deepcopy(value) for key,_,value in fields}
                if model in ['General','Elovich']:params['UsePoreSurfaceTransport']=False
                if model=='Avrami':params.update(HistoryPolicy='ContinuousErasure',Order=1.2,HistoryDrivingRegularizer=1e-7)
                if model=='DualNth':params['Branches']=[
                    {'Fraction':.35,'AdsorptionRateCoefficient':.4,'DesorptionRateCoefficient':.2,'AdsorptionOrder':1.,'DesorptionOrder':1.},
                    {'Fraction':.65,'AdsorptionRateCoefficient':.1,'DesorptionRateCoefficient':.06,'AdsorptionOrder':1.5,'DesorptionOrder':1.5}]
                d['Components'][0]['ChemisorptionSites']=[{'Type':'Langmuir','Parameters':[1.,.4 if liquid else 1e-5],
                    'MaximumLoading':1.,'HeatOfAdsorption':25000.,'RateEquation':model,'RateEquationParameters':params}]
                cases=[('breakthrough',d)]
                if model in ['Avrami','DualNth']:
                    c=deepcopy(d);c.update(SimulationType='SwingAdsorption',NumberOfCycles=2);c.pop('SimulationTime',None)
                    names=[x['Name'] for x in c['Components']]
                    concentration={n:float(c['Components'][i].get('LiquidPhaseConcentration',0.)) for i,n in enumerate(names)}
                    fractions={n:float(c['Components'][i].get('GasPhaseMolFraction',0.)) for i,n in enumerate(names)}
                    p=float(c.get('InletPressure',1e5));T=float(c.get('Temperature',298.15))
                    feed={'Pressure':p,'Temperature':T,('LiquidPhaseConcentrations' if liquid else 'GasPhaseMolFractions'):(concentration if liquid else fractions)}
                    regenerant=deepcopy(feed)
                    regenerant['LiquidPhaseConcentrations' if liquid else 'GasPhaseMolFractions']={n:(0. if liquid or i==0 else 1.) for i,n in enumerate(names)}
                    c['Sources']={'feed':feed,'regenerant':regenerant,'receiver':deepcopy(feed)}
                    velocity=float(c.get('ColumnEntranceVelocity',.01))
                    c['SwingAdsorptionPhases']=[
                        {'Name':'Adsorption','Duration':2.,'LeftPort':{'Mode':'Velocity','Source':'feed','Velocity':velocity},
                         'RightPort':{'Mode':'Pressure','Source':'receiver'},'StreamTags':{'LeftIn':'feed','RightOut':'light_product'}},
                        {'Name':'Regeneration','Duration':4.,'LeftPort':{'Mode':'Pressure','Source':'receiver'},
                         'RightPort':{'Mode':'Velocity','Source':'regenerant','Velocity':-velocity},'StreamTags':{'RightIn':'regenerant','LeftOut':'heavy_product'}}]
                    cases.append(('cycle',c))
                for operation,config in cases:
                    folder=output/method.lower()/phase/model.lower()/operation;folder.mkdir(parents=True,exist_ok=True)
                    (folder/'simulation.json').write_text(json.dumps(config,indent=2,allow_nan=False)+'\n')
                    manifest.append(str((folder/'simulation.json').relative_to(ROOT)))
    (output/'manifest.txt').write_text('\n'.join(manifest)+'\n')
    (output/'README.md').write_text('''# Shared gas/liquid chemisorption examples

These 32 synthetic examples cover FirstOrder, PseudoNth, Avrami, General,
Elovich and DualNth in gas/liquid columns with FDM/FVM. Avrami/DualNth also
include two-cycle adsorption/countercurrent-regeneration schedules. Only the
spatial selector differs between matched FDM/FVM pairs. All cases use WENO3. Rates and material
properties are illustrative; these are not fitted materials or validated recipes.

Run `build/ruptura PATH/simulation.json --validate`, then the same command without
`--validate`. The same JSON imports into Lab or ProcessSimulation. General and
Elovich examples disable additional surface/pore transport to isolate their
kinetic law. See docs/SHARED_CHEMISORPTION.md for units and scope. Regenerate with
`python tools/generate_chemisorption_examples.py`.
''')
    print(f'Generated {len(manifest)} shared chemisorption examples')

if __name__=='__main__':main()
