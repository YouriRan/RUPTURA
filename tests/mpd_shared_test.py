"""Shared MPD checks against a synthetic, exactly equivalent Langmuir model.

The distribution is a one-site occupancy model: states (0,0), (0,1),
(1,0), (1,1) have weights 1, 2, 3, 0 at 300 K and 1 bar. This is synthetic
verification data; no experimental MPD or literature validation is implied.
"""
from copy import deepcopy
import io
import json
import os
from pathlib import Path
import subprocess
import zipfile

import numpy as np
import pytest

from ruptura.process import ProcessSimulation
from ruptura.webapp.cluster import cluster_bundle
from ruptura.webapp.config import normalize_input, SCHEMA
from ruptura.webapp.runs import RunManager

ROOT=Path(__file__).resolve().parents[1]
CLI=Path(os.environ.get('RUPTURA_CLI_EXECUTABLE',ROOT/'build/ruptura'))
pytestmark=pytest.mark.skipif(not CLI.is_file(),reason='Build the native CLI first')


def document(directory,spatial='FDM',solver='CVODE'):
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'distribution.data').write_text('0.16666666666666666\n0.33333333333333331\n0.5\n0\n')
    return {'SimulationType':'Breakthrough','SpatialDiscretization':spatial,
        'MixturePredictionMethod':'MPD','Temperature':300.,'ColumnLength':.02,
        'NumberOfGridPoints':4,'ParticleDensity':1000.,'InletPressure':2.5e5,'OutletPressure':2.5e5,
        'ColumnEntranceVelocity':.05,'BoundaryCondition':'FixedVelocity',
        'Geometry':{'Type':'PackedBed','InternalDiameter':.01,'OuterDiameter':.012,
                    'ColumnVoidFraction':.4,'ParticleDiameter':.001},
        'Components':[
            {'Name':'A','GasPhaseMolFraction':.2,'MolecularWeight':.044,'GasHeatCapacity':37.,
             'MassTransferCoefficient':.2,'AxialDispersionCoefficient':1e-6,'HeatOfAdsorption':20000.},
            {'Name':'He','GasPhaseMolFraction':.2,'MolecularWeight':.004,'GasHeatCapacity':20.8,'CarrierGas':True},
            {'Name':'B','GasPhaseMolFraction':.6,'MolecularWeight':.028,'GasHeatCapacity':29.,
             'MassTransferCoefficient':.1,'AxialDispersionCoefficient':1e-6,'HeatOfAdsorption':18000.}],
        'MPDSettings':{'FileName':'distribution.data','ReferenceTemperature':300.,'ReferenceFugacity':1e5,
                       'ReferenceFrameworkMass':1./6.02214076e23,'ComponentBounds':[
                           {'Component':n,'NMin':0,'NMax':1,'DeltaN':1} for n in ['A','B']]},
        'InitialConditions':{'Type':'CarrierGas','Loading':'zero'},
        'BreakthroughIntegrator':solver,'ConvectionScheme':'TVD','SimulationTime':.01,'WriteInterval':.002,
        'Solver':{'MaximumTimeStep':.0005,'FirstTimeStep':1e-5,'RelativeTolerance':1e-7,'AbsoluteTolerance':1e-10},
        'Output':{'Directory':'results','Plots':False,'Checkpoint':False,'MainFormatFiles':False}}


def analytical(d):
    d=deepcopy(d);d['MixturePredictionMethod']='SCI';d.pop('MPDSettings')
    d['Components'][0]['PhysisorptionSites']=[{'Type':'Langmuir','Parameters':[1.,3e-5]}]
    d['Components'][2]['PhysisorptionSites']=[{'Type':'Langmuir','Parameters':[1.,2e-5]}]
    return d


def run(d,directory,*options):
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/'simulation.json';path.write_text(json.dumps(d))
    p=subprocess.run([str(CLI),str(path),'--no-plots',*options],capture_output=True,text=True,timeout=60)
    assert p.returncode==0,p.stderr+p.stdout[-2000:]
    with np.load(directory/'results/states.npz',allow_pickle=False) as f:
        arrays={k:f[k].copy() for k in f.files}
    summary=json.loads((directory/'results/summary.json').read_text())
    assert np.isfinite(arrays['state']).all()
    assert abs(summary['final_metrics']['mass_balance_error'])<1e-6
    return arrays,summary


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('solver',['BDF1','BDF','CVODE','SSPRK3','SemiImplicitSSPRK3'])
def test_mpd_executes_with_each_shared_integrator(tmp_path,spatial,solver):
    directory=tmp_path/'mpd';d=document(directory,spatial,solver)
    actual,_=run(d,directory);expected,_=run(analytical(d),tmp_path/'analytical')
    np.testing.assert_allclose(actual['state'],expected['state'],rtol=2e-6,atol=2e-8)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_mpd_monolith_matches_analytic_equilibrium(tmp_path,spatial):
    directory=tmp_path/'mpd';d=document(directory,spatial)
    d['Geometry'].update(Type='Monolith',ChannelShape='Square',InternalChannelDimension=.001,
                         NumberOfChannels=16,WashcoatThickness=.0001)
    actual,_=run(d,directory);expected,_=run(analytical(d),tmp_path/'analytical')
    np.testing.assert_allclose(actual['state'],expected['state'],rtol=2e-6,atol=2e-8)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('shape',['Sphere','CylinderRadial'])
@pytest.mark.parametrize('film',[False,True])
def test_mpd_hsdm_with_and_without_film_matches_analytic_equilibrium(tmp_path,spatial,shape,film):
    directory=tmp_path/'mpd';d=document(directory,spatial)
    d['InitialConditions']={'Type':'FeedPressurized','Loading':'zero'}
    d['ParticleTransport']={'Model':'HSDM','Geometry':shape,'Radius':.0005,'RadialCells':4,
        'SurfaceDiffusivity':{'A':1e-10,'He':0.,'B':2e-10},'FilmResistance':film,
        'FilmCoefficient':{'A':2e-5,'He':0.,'B':1e-5}}
    actual,_=run(d,directory);expected,_=run(analytical(d),tmp_path/'analytical')
    np.testing.assert_allclose(actual['state'],expected['state'],rtol=2e-5,atol=2e-8)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_complete_synthetic_mpd_breakthrough_matches_analytic_curve(tmp_path,spatial):
    directory=tmp_path/'mpd';d=document(directory,spatial)
    d.update(NumberOfGridPoints=10,SimulationTime=60.,WriteInterval=.5)
    d['Solver']['MaximumTimeStep']=.1
    actual,_=run(d,directory);expected,_=run(analytical(d),tmp_path/'analytical')
    np.testing.assert_allclose(actual['state'],expected['state'],rtol=3e-5,atol=3e-6)
    final=actual['state'][-1,:33].reshape(3,11)
    outlet=final[:,-1]/final[:,-1].sum()
    np.testing.assert_allclose(outlet,[.2,.2,.6],rtol=.01,atol=.002)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_mpd_shared_cycles_use_the_same_equilibrium(tmp_path,spatial):
    directory=tmp_path/'mpd';d=document(directory,spatial)
    d['SimulationType']='SwingAdsorption';duration=d.pop('SimulationTime')
    d['NumberOfCycles']=2
    d['SwingAdsorptionPhases']=[{'Name':'Adsorption','Duration':duration}]
    actual,_=run(d,directory);expected,_=run(analytical(d),tmp_path/'analytical')
    np.testing.assert_allclose(actual['state'],expected['state'],rtol=2e-6,atol=2e-8)


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_changed_distribution_invalidates_restart(tmp_path,spatial):
    d=document(tmp_path,spatial);d['Output']['Checkpoint']=True
    d['SimulationType']='SwingAdsorption';duration=d.pop('SimulationTime')
    d['SwingAdsorptionPhases']=[{'Name':'Adsorption','Duration':duration}];d['NumberOfCycles']=1
    run(d,tmp_path)
    path=tmp_path/'simulation.json'
    unchanged=subprocess.run([str(CLI),str(path),'--restart','--cycles','2','--no-plots'],
                             capture_output=True,text=True,timeout=30)
    assert unchanged.returncode==0,unchanged.stderr
    (tmp_path/'distribution.data').write_text('0.1\n0.2\n0.7\n0\n')
    p=subprocess.run([str(CLI),str(path),'--restart','--cycles','3','--no-plots'],capture_output=True,text=True,timeout=30)
    assert p.returncode!=0
    assert 'checkpoint incompatible' in p.stderr.lower()


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_lab_json_and_cluster_bundle_execute_the_same_mpd_configuration(tmp_path,spatial):
    directory=tmp_path/'source';d=document(directory,spatial);before=deepcopy(d)
    manager=RunManager(directory,CLI)
    try:canonical=manager.validate(d,directory)
    finally:manager.close()
    assert canonical['MixturePredictionMethod']=='MPD'
    assert 'MPD' in next(f for f in SCHEMA['operating'] if f['key']=='MixturePredictionMethod')['options']
    assert Path(canonical['MPDSettings']['FileName']).is_absolute()
    assert normalize_input(json.loads(json.dumps(canonical)))==canonical
    archive=cluster_bundle(canonical)
    target=tmp_path/'cluster';target.mkdir()
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        exported=json.loads(z.read('simulation.json'))
        assert z.read(exported['MPDSettings']['FileName'])==(directory/'distribution.data').read_bytes()
        z.extractall(target)
    compare=deepcopy(exported);compare['MPDSettings']['FileName']=canonical['MPDSettings']['FileName']
    assert compare==canonical
    p=subprocess.run(['bash','run.sh','--no-plots'],cwd=target,env={**os.environ,'RUPTURA_EXECUTABLE':str(CLI)},
                     capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr+p.stdout[-2000:]
    expected,_=run(d,directory)
    with np.load(target/'results/states.npz',allow_pickle=False) as f:
        np.testing.assert_allclose(f['state'],expected['state'],rtol=1e-11,atol=1e-11)
    assert d==before


@pytest.mark.parametrize('spatial',['FDM','FVM'])
def test_mpd_python_interface_matches_cli(tmp_path,spatial):
    pytest.importorskip('ruptura._ruptura')
    directory=tmp_path/'cli';d=document(directory,spatial)
    expected,_=run(d,directory)
    d['MPDSettings']['FileName']=str(directory/'distribution.data')
    simulation=ProcessSimulation(d,directory=tmp_path/'python')
    assert simulation.validate()['valid']
    simulation.run(plots=False)
    np.testing.assert_allclose(simulation.results()['state'],expected['state'],rtol=1e-11,atol=1e-11)
    assert 'MPD' in simulation.capabilities()['mixture_methods']


@pytest.mark.parametrize('spatial',['FDM','FVM'])
@pytest.mark.parametrize('invalid',['missing_settings','wrong_bounds_order','wrong_rows','negative_probability','liquid'])
def test_invalid_mpd_data_are_rejected(tmp_path,spatial,invalid):
    d=document(tmp_path,spatial)
    if invalid=='missing_settings':d.pop('MPDSettings')
    elif invalid=='wrong_bounds_order':d['MPDSettings']['ComponentBounds'][0]['Component']='B'
    elif invalid=='wrong_rows':(tmp_path/'distribution.data').write_text('0.1\n0.2\n0.3\n')
    elif invalid=='negative_probability':(tmp_path/'distribution.data').write_text('0.1\n-0.2\n0.3\n0\n')
    elif invalid=='liquid':d['FluidPhase']='Liquid'
    path=tmp_path/'simulation.json';path.write_text(json.dumps(d))
    p=subprocess.run([str(CLI),str(path),'--validate'],capture_output=True,text=True,timeout=30)
    assert p.returncode!=0 and 'MPD' in p.stderr
