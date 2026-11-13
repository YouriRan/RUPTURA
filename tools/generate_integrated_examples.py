#!/usr/bin/env python3
"""Create independent examples using the public RUPTURA-main JSON conventions."""
from pathlib import Path
import copy,json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'examples'/'new_examples' / 'integrated'
def base(liquid=False,monolith=False):
    a={'Name':'A','MolecularWeight':0.044,'GasHeatCapacity':37.0,'GasPhaseMolFraction':0.1,'LiquidPhaseConcentration':0.5,'InitialLiquidPhaseConcentration':0.0,'MassTransferCoefficient':0.2,'AxialDispersionCoefficient':1e-6,'HeatOfAdsorption':20000.,'PhysisorptionSites':[{'Type':'Langmuir','Parameters':[2.,.4 if liquid else 1e-5]}]}
    b={'Name':'B','MolecularWeight':.028,'GasHeatCapacity':29.,'GasPhaseMolFraction':.9,'CarrierGas':True,'AxialDispersionCoefficient':1e-6}
    geometry={'Type':'Monolith' if monolith else 'PackedBed','InternalDiameter':.01,'OuterDiameter':.012,'ColumnVoidFraction':.4,'ParticleDiameter':.001}
    if monolith:geometry.update(ChannelDiameter=.0005,NumberOfChannels=100,ChannelShape='circular',WashcoatThickness=.00005)
    return {'SimulationType':'Breakthrough','SpatialDiscretization':'FVM','FluidPhase':'Liquid' if liquid else 'Gas','DisplayName':f'{"Liquid" if liquid else "Gas"} {"monolith" if monolith else "fixed bed"}', 'Components':[a,b],'MixturePredictionMethod':'SCI','Geometry':geometry,'ColumnLength':.1,'ParticleDensity':1000.,'NumberOfGridPoints':20,'ColumnEntranceVelocity':.01,'Temperature':298.15,'InletPressure':100000.,'OutletPressure':100000.,'BoundaryCondition':'InletVelocityOutletPressure','BreakthroughIntegrator':'CVODE','CVODELinearSolver':'Dense','ConvectionScheme':'WENO3','SimulationTime':40.,'WriteInterval':.5,'Solver':{'MaximumTimeStep':.5,'RelativeTolerance':1e-6},'Performance':{'LightComponent':'A','HeavyComponent':'B'},'Output':{'Directory':'results','Profile':'compact','Plots':True,'MainFormatFiles':True}}
def write(name,doc,description):
    doc['ConvectionScheme']='WENO3'
    doc.setdefault('Output',{}).update(Plots=True)
    path=OUT/name;path.mkdir(parents=True,exist_ok=True)
    (path/'simulation.json').write_text(json.dumps(doc,indent=2)+'\n')
    run='''#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$case_dir"
root="$case_dir"
while [[ ! -f "$root/CMakeLists.txt" || ! -d "$root/ruptura" ]]; do
  parent="$(dirname -- "$root")"
  [[ "$parent" != "$root" ]] || { echo "Cannot locate the RUPTURA project" >&2; exit 1; }
  root="$parent"
done
engine="${RUPTURA_EXECUTABLE:-$root/build/ruptura}"
[[ -x "$engine" ]] || { echo "Build RUPTURA first, or set RUPTURA_EXECUTABLE" >&2; exit 1; }
exec "$engine" "$@"
'''
    (path/'run').write_text(run);(path/'run').chmod(0o755)
    (path/'submit_sim_hal.sh').write_text('''#!/usr/bin/env bash
#SBATCH --partition=serial
#SBATCH --exclude=c171
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname -- "${BASH_SOURCE[0]}")}"
exec ./run "$@"
''');(path/'submit_sim_hal.sh').chmod(0o755)
    (path/'README.md').write_text(description+'\n\nAll pressures are absolute Pa; times are seconds; concentrations are mol/m³; loadings are mol/kg adsorbent. Velocity is interstitial. Run `./run --validate`, then `./run`. The parameters are illustrative, not experimental validation.\n')
    return path
def main():
    manifest=[]
    for liquid in [False,True]:
      # The four packed-bed cases have one canonical home in unified/fvm.
      for monolith in [True]:
        phase='liquid' if liquid else 'gas';kind='monolith' if monolith else 'fixed_bed'
        doc=base(liquid,monolith);name=f'breakthrough/{phase}_{kind}';write(name,doc,'Single-column FVM breakthrough with conservative species transport.');manifest.append(name)
        doc=base(liquid,monolith);doc['SimulationType']='SwingAdsorption';doc.pop('SimulationTime');doc['NumberOfCycles']=3
        feed={'Pressure':100000.,'Temperature':298.15,'GasPhaseMolFractions':{'A':.1,'B':.9}}
        regeneration=copy.deepcopy(feed);regeneration['GasPhaseMolFractions']={'A':0.,'B':1.}
        if liquid:feed['LiquidPhaseConcentrations']={'A':.5,'B':0.};regeneration['LiquidPhaseConcentrations']={'A':0.,'B':0.}
        doc['Sources']={'feed':feed,'regenerant':regeneration,'product':copy.deepcopy(feed)}
        doc['SwingAdsorptionPhases']=[{'Name':'Adsorption','Duration':10.,'LeftPort':{'Mode':'Velocity','Source':'feed','Velocity':.01},'RightPort':{'Mode':'Pressure','Source':'product'},'StreamTags':{'LeftIn':'feed','RightOut':'light_product'}},{'Name':'Regeneration','Duration':20.,'LeftPort':{'Mode':'Pressure','Source':'product'},'RightPort':{'Mode':'Velocity','Source':'regenerant','Velocity':-.005},'StreamTags':{'RightIn':'regenerant','LeftOut':'heavy_product'}}]
        name=f'cycles/{phase}_{kind}';write(name,doc,'Adsorption followed by countercurrent regeneration. Liquid regeneration changes the desorbent concentration; the solvent is implicit.');manifest.append(name)
    for liquid in [False,True]:
     for shape in ['Sphere','CylinderRadial']:
        doc=base(liquid,False);doc['CVODELinearSolver']='SPGMR';doc['ParticleTransport']={'Model':'HSDM','Geometry':shape,'Radius':.0005,'RadialCells':12,'SurfaceDiffusivity':{'A':5e-11,'B':0.},'referenceTemperature':298.15,'FilmResistance':True,'FilmCoefficient':{'A':5e-5,'B':0.}}
        if liquid:doc['MixturePredictionMethod']='SPI';doc['Components'][0]['PhysisorptionSites']=[{'Type':'Freundlich','Parameters':[1.,.3],'ExtraParameters':{'local_c1':True}}]
        name=f'particle/{"liquid" if liquid else "gas"}_hsdm_{shape.lower()}';write(name,doc,'Physisorption HSDM with coupled finite-film surface equilibrium. Radial cylinder transport excludes the cylinder ends.');manifest.append(name)
    for liquid in [False,True]:
        doc=base(liquid,False);doc['energyBalance']=True;doc.update(heatCapacitySolid=900.,wallDensity=7800.,heatCapacityWall=450.,heatTransferGasSolid=20.,heatTransferGasWall=10.,wallThermalConductivity=10.)
        doc['Coolant']={'Model':'Transient','Volume':.0001,'Density':1000.,'HeatCapacity':4180.,'InitialTemperature':298.15,'InletTemperature':293.15,'MassFlowRate':.01,'HeatTransferCoefficient':50.,'FlowDirection':'CounterCurrent'}
        name=f'thermal/{"liquid" if liquid else "gas"}_transient_coolant';write(name,doc,'Separate fluid, solid, wall and distributed coolant temperatures.');manifest.append(name)
    doc=base();doc['SimulationTime']='auto';doc['AutomaticStopping']={'Components':['A','B'],'Threshold':.999,'LoadingRelativeTolerance':.001,'MaximumTime':2000.,'MinimumTime':1.,'ExtraFraction':.1};write('runtime/automatic_saturation',doc,'Automatic ending requires both outlet arrival and loading equilibrium for every selected fed component, followed by an observation interval.');manifest.append('runtime/automatic_saturation')
    # Retain the legacy formulation with the supported boundary and explicit WENO3.
    t15=json.loads((ROOT/'marx2/ternary/T15/simulation.json').read_text());t15['BoundaryCondition']='FixedVelocity';write('validation/marx_T15_original',t15,'Legacy-formulation Marx T15 comparison with FixedVelocity and WENO3. The supplied source input remains under marx2/ternary/T15.');manifest.append('validation/marx_T15_original')
    t15['SpatialDiscretization']='FVM';t15['ConvectionScheme']='WENO3';t15['OutletPressure']=t15['InletPressure'];t15['Solver']={'MaximumTimeStep':2.,'NonnegativeConstraints':True,'IsothermPowerRegularizer':0.0};t15['CVODEMaximumOrder']=2;t15['Output']={'Plots':True,'Directory':'results','Profile':'compact','MainFormatFiles':True}
    cp={'Helium':20.786156545,'CO2':37.,'N2':29.1,'H2':28.8}
    for c in t15['Components']:c['GasHeatCapacity']=cp[c['Name']]
    write('validation/marx_T15_fvm',t15,'Explicit FVM variant of Marx T15. Component heat capacities and an outlet reservoir complete the shared compressible model; it is a comparison case, not a claim of equation identity with the original main formulation.');manifest.append('validation/marx_T15_fvm')
    (OUT/'cases.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Generated {len(manifest)} independent JSON examples')

if __name__=="__main__":main()
