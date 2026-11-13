"""Conservation, regeneration, migration and stateful kinetics for both fluids."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = Path(os.environ.get('RUPTURA_CLI_EXECUTABLE', ROOT/'build/ruptura'))
pytestmark = pytest.mark.skipif(not CLI.is_file(), reason='Build the current CLI first')


def site(model, liquid=False):
    parameters = {
        'FirstOrder': {'MassTransferCoefficient': .4},
        'PseudoNth': {'RateCoefficient': .4, 'Order': 1.5},
        'Avrami': {'RateCoefficient': .4, 'Order': 1.2, 'HistoryPolicy': 'ContinuousErasure'},
        'General': {'AdsorptionRateCoefficient': .4, 'AdsorptionActivationEnergy': 0.,
                    'DesorptionRateCoefficient': .1, 'DesorptionActivationEnergy': 0.,
                    'PoreConcentrationOrder': 1, 'CapacityOrder': 1, 'DesorptionOrder': 1,
                    'FilmMassTransferCoefficient': 0., 'PoreDiffusivity': 0., 'UsePoreSurfaceTransport': False},
        'Elovich': {'Alpha': .1, 'Beta': 2., 'FilmMassTransferCoefficient': 0.,
                    'PoreDiffusivity': 0., 'UsePoreSurfaceTransport': False},
        'DualNth': {'ReferenceTemperature': 298.15, 'Branches': [
            {'Fraction': .35, 'AdsorptionRateCoefficient': .4, 'DesorptionRateCoefficient': .2,
             'AdsorptionOrder': 1., 'DesorptionOrder': 1., 'HeatOfAdsorption': 20000.},
            {'Fraction': .65, 'AdsorptionRateCoefficient': .1, 'DesorptionRateCoefficient': .06,
             'AdsorptionOrder': 1.5, 'DesorptionOrder': 1.5, 'HeatOfAdsorption': 35000.}]},
    }[model]
    return {'Type':'Langmuir', 'Parameters':[1., 1. if liquid else 1./(8.31446261815324*298.15)],
            'MaximumLoading': 1., 'HeatOfAdsorption': 25000., 'RateEquation':model,
            'RateEquationParameters':parameters}


def document(spatial='FDM', liquid=False, model='DualNth'):
    d={'SimulationType':'SwingAdsorption', 'SpatialDiscretization':spatial,
       'FluidPhase':'Liquid' if liquid else 'Gas', 'MixturePredictionMethod':'SPI',
       'Temperature':298.15, 'ColumnLength':.02, 'NumberOfGridPoints':4, 'ParticleDensity':10.,
       'InletPressure':2.*8.31446261815324*298.15, 'OutletPressure':2.*8.31446261815324*298.15,
       'BoundaryCondition':'InletVelocityOutletPressure', 'ColumnEntranceVelocity':.01,
       'Geometry':{'Type':'PackedBed','InternalDiameter':.01,'OuterDiameter':.012,
                   'ColumnVoidFraction':.4,'ParticleDiameter':.001},
       'Components':[{'Name':'A','MolecularWeight':.044,'GasHeatCapacity':37.,'GasPhaseMolFraction':.2,
                      'LiquidPhaseConcentration':.4,'InitialLiquidPhaseConcentration':.4,
                      'MassTransferCoefficient':.4, 'AxialDispersionCoefficient':0.,
                      'PhysisorptionSites':[{'Type':'Henry','Parameters':[0.]}],
                      'ChemisorptionSites':[site(model,liquid)]},
                     {'Name':'B','MolecularWeight':.028,'GasHeatCapacity':29.,'GasPhaseMolFraction':.8,
                      'LiquidPhaseConcentration':1.6,'InitialLiquidPhaseConcentration':1.6,
                      'PhysisorptionSites':[{'Type':'Henry','Parameters':[0.]}]}],
       'PhysisorptionEnabled':False, 'InitialConditions':{'Type':'Custom','Loading':'zero',
          'GasPhaseMolFractions':{'A':.2,'B':.8},'LiquidPhaseConcentrations':{'A':.4,'B':1.6}},
       'NumberOfCycles':1, 'SwingAdsorptionPhases':[{'Name':'Closed relaxation','Duration':.1,
          'Columns':{'A':{'LeftPort':{'Mode':'closed'},'RightPort':{'Mode':'closed'}}}}],
       'BreakthroughIntegrator':'CVODE','ConvectionScheme':'TVD','WriteInterval':.01,
       'CVODEAbsoluteToleranceLoading':1e-12,'CVODEAbsoluteToleranceConcentration':1e-12,
       'Solver':{'RelativeTolerance':1e-9,'AbsoluteTolerance':1e-12,'MaximumTimeStep':.001},
       'Output':{'Directory':'results','Plots':False,'Checkpoint':False,'MainFormatFiles':False}}
    return d


def execute(d, directory, *flags):
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/'simulation.json';path.write_text(json.dumps(d,allow_nan=False))
    run=subprocess.run([str(CLI),str(path),'--no-plots',*flags],capture_output=True,text=True,timeout=45)
    assert run.returncode==0,run.stderr+run.stdout[-2000:]
    with np.load(directory/'results/states.npz',allow_pickle=False) as f:
        data={k:f[k].copy() for k in f.files}
    summary=json.loads((directory/'results/summary.json').read_text())
    metadata=json.loads((directory/'results/metadata.json').read_text())
    assert np.isfinite(data['state']).all()
    assert abs(summary['final_metrics']['mass_balance_error'])<1e-6
    return data,summary,metadata


def chemical_groups(data, metadata):
    record=metadata['column_records']['columns']['A']
    nc=len(record['component_names']);nz=len(record['z']);block=nc*nz
    groups=[i for i,name in enumerate(record['loading_group_names']) if '/chemisorption/' in name]
    return [data['state'][:,(1+i)*block:(2+i)*block].reshape(-1,nc,nz) for i in groups]


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('model',['FirstOrder','PseudoNth','Avrami','General','Elovich','DualNth'])
def test_closed_adsorption_conserves_and_uses_one_interface(tmp_path,spatial,liquid,model):
    d=document(spatial,liquid,model)
    data,summary,metadata=execute(d,tmp_path/'run')
    groups=chemical_groups(data,metadata)
    assert sum(g[-1,0].mean() for g in groups)>0.
    assert len(groups)==(2 if model=='DualNth' else 1)
    assert 'ChemicalKinetics' not in metadata['input_document']


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('model',['FirstOrder','PseudoNth','Avrami','DualNth'])
def test_regeneration_releases_loading_for_integer_and_fractional_orders(tmp_path,spatial,liquid,model):
    d=document(spatial,liquid,model);d['InitialConditions']['Loading']={'A':.8,'B':0.}
    if model=='PseudoNth':d['Components'][0]['ChemisorptionSites'][0]['RateEquationParameters']['Order']=2.
    data,_,metadata=execute(d,tmp_path/'run');groups=chemical_groups(data,metadata)
    assert sum(g[-1,0].mean() for g in groups)<.8
    assert all(g[-1,0].mean()<g[0,0].mean() for g in groups)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('model',['FirstOrder','PseudoNth','Avrami','DualNth'])
def test_matching_gas_and_liquid_targets_use_identical_kinetics(tmp_path,spatial,model):
    gas,_,_=execute(document(spatial,False,model),tmp_path/'gas')
    liquid,_,_=execute(document(spatial,True,model),tmp_path/'liquid')
    np.testing.assert_allclose(gas['state'],liquid['state'],rtol=2e-7,atol=2e-10)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('model',['FirstOrder','PseudoNth','Avrami','General','Elovich','DualNth'])
def test_branch_heats_conserve_adiabatic_energy(tmp_path,spatial,liquid,model):
    d=document(spatial,liquid,model);d.update(energyBalance=True,ThermalModel='Adiabatic',
       heatCapacitySolid=900.,wallDensity=7800.,heatCapacityWall=450.,
       heatTransferGasSolid=20.,heatTransferGasWall=0.,wallThermalConductivity=0.)
    _,summary,_=execute(d,tmp_path/'run')
    assert abs(summary['final_metrics']['energy_balance_error'])<2e-6


@pytest.mark.parametrize('model',['pfo','pso','avrami','dual_nth'])
@pytest.mark.parametrize('liquid',[False,True])
def test_legacy_import_and_canonical_input_execute_identically(tmp_path,model,liquid):
    d=document(liquid=liquid);d['Components'][0].pop('ChemisorptionSites')
    d['ChemicalKinetics']={'Model':model,'KAds':[.4,0.], 'KDes':[.2,0.],
        'Equilibrium':{'Model':'Henry','K':[1. if liquid else 1./(8.31446261815324*298.15),0.]}}
    if model=='dual_nth':d['ChemicalKinetics'].update(Fractions=[.35,.65],OrderAds=1.5,OrderDes=1.2)
    if model=='avrami':d['ChemicalKinetics'].update(NAds=1.2,NDes=1.5)
    path=tmp_path/'legacy.json';path.write_text(json.dumps(d))
    p=subprocess.run([str(CLI),str(path),'--canonicalize'],capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    canonical=json.loads(p.stdout);assert 'ChemicalKinetics' not in canonical
    old,_,_=execute(d,tmp_path/'old');new,_,_=execute(canonical,tmp_path/'new')
    np.testing.assert_array_equal(old['state'],new['state'])


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('model',['Avrami','DualNth'])
def test_stage_splits_and_restart_preserve_kinetic_state(tmp_path,spatial,liquid,model):
    d=document(spatial,liquid,model)
    whole,_,_=execute(d,tmp_path/'whole')
    split=deepcopy(d);phase=split['SwingAdsorptionPhases'][0];phase['Duration']=.05
    split['SwingAdsorptionPhases'].append({**deepcopy(phase),'Name':'Continued relaxation'})
    divided,_,_=execute(split,tmp_path/'split')
    np.testing.assert_allclose(whole['state'][-1],divided['state'][-1],rtol=2e-6,atol=2e-9)
    d['Output']['Checkpoint']=True
    execute(d,tmp_path/'restart')
    d['NumberOfCycles']=2
    restarted,_,_=execute(d,tmp_path/'restart','--restart')
    uninterrupted,_,_=execute(d,tmp_path/'uninterrupted')
    np.testing.assert_allclose(restarted['state'][-1],uninterrupted['state'][-1],rtol=2e-6,atol=2e-9)


@pytest.mark.parametrize('change,message',[
    (lambda p:p['Branches'][0].update(Fraction=.9),'sum to one'),
    (lambda p:p['Branches'][0].update(AdsorptionOrder=0.),'order'),
    (lambda p:p['Branches'][0].update(AdsorptionRateCoefficient=-1.),'rate'),
    (lambda p:p.update(Branches=p['Branches'][:1]),'exactly two'),
    (lambda p:p.update(RateBasis='nonsense'),'RateBasis'),
    (lambda p:p.update(ReferenceTemperature=0.),'reference temperature'),
])
def test_invalid_branch_models_fail_clearly(tmp_path,change,message):
    d=document(liquid=True);change(d['Components'][0]['ChemisorptionSites'][0]['RateEquationParameters'])
    path=tmp_path/'simulation.json';path.write_text(json.dumps(d))
    p=subprocess.run([str(CLI),str(path),'--validate'],capture_output=True,text=True,timeout=30)
    assert p.returncode!=0
    assert message.lower() in (p.stderr+p.stdout).lower()


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('model',['FirstOrder','PseudoNth','Avrami','DualNth'])
def test_breakthrough_and_one_phase_cycle_use_same_chemical_states(tmp_path,spatial,liquid,model):
    cycle=document(spatial,liquid,model)
    cycle['SwingAdsorptionPhases']=[{'Name':'Adsorption','Duration':.02}]
    bt=deepcopy(cycle);bt.pop('SwingAdsorptionPhases');bt.pop('NumberOfCycles')
    bt.update(SimulationType='Breakthrough',SimulationTime=.02)
    a,_,_=execute(bt,tmp_path/'bt');b,_,_=execute(cycle,tmp_path/'cycle')
    np.testing.assert_allclose(a['state'],b['state'],rtol=1e-11,atol=1e-11)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('liquid',[False,True])
@pytest.mark.parametrize('solver',['BDF1','BDF','CVODE','SSPRK3','SemiImplicitSSPRK3'])
def test_dual_branches_execute_with_every_shared_integrator(tmp_path,spatial,liquid,solver):
    d=document(spatial,liquid);d['BreakthroughIntegrator']=solver
    d['SwingAdsorptionPhases'][0]['Duration']=.005;d['WriteInterval']=.001
    d['Solver']['MaximumTimeStep']=.0001
    data,_,metadata=execute(d,tmp_path/'run')
    assert all(g[-1,0].mean()>0. for g in chemical_groups(data,metadata))


@pytest.mark.parametrize('liquid',[False,True])
def test_imported_capacity_partition_preserves_total_equilibrium(tmp_path,liquid):
    from ruptura.process import canonical_chemisorption
    from ruptura.webapp.cluster import cluster_bundle
    import io,zipfile
    d=document(liquid=liquid);d['Components'][0].pop('ChemisorptionSites')
    d['Components'][0]['PhysisorptionSites']=[{'Type':'Langmuir','Parameters':[1.,1. if liquid else 1./(8.31446261815324*298.15)]}]
    d['ChemicalKinetics']={'Model':'dual_nth','KAds':.4,'KDes':.2,'Fractions':[.35,.65],'CapacityFraction':.3}
    d['InitialConditions']['Loading']='equilibrium'
    canonical=canonical_chemisorption(d)
    assert canonical['PhysisorptionFraction']==.7
    data,_,metadata=execute(canonical,tmp_path/'run')
    np.testing.assert_allclose(sum(g[0,0] for g in chemical_groups(data,metadata)),.3*.4/1.4,rtol=1e-12)
    with zipfile.ZipFile(io.BytesIO(cluster_bundle(d))) as archive:
        assert json.loads(archive.read('simulation.json'))==canonical
    assert 'ChemicalKinetics' in d


@pytest.mark.parametrize('liquid',[False,True])
def test_nonstandard_reference_temperatures_match_closed_henry_solution(tmp_path,liquid):
    d=document(liquid=liquid);d['Components'][0].pop('ChemisorptionSites')
    R=8.31446261815324;T=298.15;gamma=15.
    d['ChemicalKinetics']={'Model':'pfo','KAds':[.4,0.],'KDes':[.2,0.],
        'EaAds':[20000.,0.],'ReferenceTemperature':314.,'HeatOfAdsorption':[25000.,0.],
        'Equilibrium':{'Model':'Henry','K':[.5 if liquid else .5/(R*T),0.],
            'ReferenceTemperature':320.,'HeatOfAdsorption':[12000.,0.]}}
    data,_,metadata=execute(d,tmp_path/'run')
    chemical=metadata['input_document']['Components'][0]['ChemisorptionSites'][0]
    assert chemical['ReferenceTemperature']==320.
    assert chemical['AdsorptionReferenceTemperature']==320.
    assert chemical['RateEquationParameters']['ReferenceTemperature']==314.
    alpha=.5*np.exp(12000./R*(1./T-1./320.));rate=.4*np.exp(20000./R*(1./314.-1./T))
    target=alpha*.4/(1.+alpha*gamma)*(1.-np.exp(-rate*(1.+alpha*gamma)*.1))
    actual=sum(g[-1,0] for g in chemical_groups(data,metadata))
    np.testing.assert_allclose(actual,target,rtol=2e-6,atol=2e-10)
