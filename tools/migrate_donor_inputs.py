#!/usr/bin/env python3
"""Convert the supplied standalone C++ packages to the public main JSON style.

Usage: python tools/migrate_donor_inputs.py CYCLE_PROJECT BREAKTHROUGH_PROJECT
Generated inputs are self-contained. PyYAML is needed for conversion only.
The native solver always reads simulation.json.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
import shutil
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
R = 8.31446261815324
PORT_KEYS = {
    'name':'Name','columns':'Columns','connections':'Connections','from':'From','to':'To',
    'left_port':'LeftPort','right_port':'RightPort','mode':'Mode','source':'Source','sink':'Sink',
    'pressure_Pa':'Pressure','pressure':'Pressure','temperature_K':'Temperature','temperature':'Temperature',
    'velocity_m_s':'Velocity','velocity':'Velocity','molar_flow_mol_s':'MolarFlowRate',
    'volume_flow_m3_s':'VolumetricFlowRate','composition':'GasPhaseMolFractions',
    'concentration':'LiquidPhaseConcentrations','composition_source':'CompositionSource',
    'fallback_composition':'FallbackComposition','stream_tags':'StreamTags',
    'left_in':'LeftIn','left_out':'LeftOut','right_in':'RightIn','right_out':'RightOut',
    'duration_s':'Duration','duration':'Duration','stop':'StopCriterion',
    'maximum_duration_s':'MaximumDuration','minimum_duration_s':'MinimumDuration',
    'maximum_duration':'MaximumDuration','minimum_duration':'MinimumDuration',
    'failure_maximum_simulated_time_s':'FailureMaximumSimulatedTime',
    'timeout_policy':'TimeoutPolicy','thermal':'ThermalControl','model':'Model',
    'initial_value':'InitialValue','final_value':'FinalValue','ramp_time_s':'RampTime',
    'time_offset_s':'TimeOffset','decay_rate_s_inv':'DecayRate', 'times_s':'Times','values':'Values',
    'mass_flow_rate_kg_s':'MassFlowRate','heat_capacity_J_kg_K':'HeatCapacity',
    'heat_transfer_coefficient_W_m2_K':'HeatTransferCoefficient','flow_direction':'FlowDirection',
    'power_W':'Power','heat_flux_W_m2':'HeatFlux','flow_model':'FlowModel',
    'conductance_mol_s_Pa':'Conductance','conductance':'Conductance','opening':'Opening',
    'check_valve':'CheckValve','maximum_molar_flow':'MaximumMolarFlow',
    'area_m2':'Area','Cd':'DischargeCoefficient','gamma':'HeatCapacityRatio',
    'tag':'Tag','in_tag':'InTag','out_tag':'OutTag','compressor':'Compressor',
    'efficiency':'Efficiency','discharge_pressure':'DischargePressure',
    'metric':'Metric','operator':'Operator','value':'Value','column':'Column','component':'Component',
    'side':'Side','port':'Port','logic':'Logic','criteria':'Criteria',
    'reference_concentration':'ReferenceConcentration','reference_pressure':'ReferencePressure',
    'reference_temperature':'ReferenceTemperature','buffer':'Buffer','volume_m3':'Volume',
    'volume':'Volume','initial_pressure_Pa':'InitialPressure','initial_pressure':'InitialPressure',
    'initial_temperature_K':'InitialTemperature','initial_temperature':'InitialTemperature',
    'initial_composition':'InitialComposition','initial_concentration':'InitialConcentration',
    'thermal_mode':'ThermalMode','UA_W_K':'UA','ambient_temperature':'AmbientTemperature',
    'pressure_control':'PressureControl','setpoint_pressure':'SetpointPressure',
    'conductance_mol_s_pa':'ControllerConductance','deadband_pa':'Deadband',
    'export_tag':'ExportTag','import_tag':'ImportTag','supply_source':'SupplySource',
}

def process_tree(x, *, named=False):
    if isinstance(x, list): return [process_tree(v) for v in x]
    if not isinstance(x, dict): return x
    result={}
    for key,value in x.items():
        target=key if named else PORT_KEYS.get(key)
        if target is None: raise ValueError(f'Unmapped process setting {key}')
        if key in ('composition','fallback_composition','initial_composition','concentration','initial_concentration'):
            result[target]=copy.deepcopy(value)
        else: result[target]=process_tree(value,named=key=='columns')
    return result

CSS_KEYS={'state_rtol':'StateRelativeTolerance','state_atol':'StateAbsoluteTolerance',
    'inventory_rtol':'InventoryRelativeTolerance','metric_rtol':'MetricRelativeTolerance',
    'consecutive_cycles':'ConsecutiveCycles','diagnostic_window':'DiagnosticWindow',
    'diagnostic_minimum_cycles':'DiagnosticMinimumCycles','maximum_period':'MaximumPeriod',
    'plot_cycle_selection':'PlotCycleSelection','mass_balance_tolerance':'MassBalanceTolerance'}
UTILITY_KEYS={'ambient_pressure_Pa':'AmbientPressure','feed_supply_pressure_Pa':'FeedSupplyPressure',
    'product_delivery_pressure_Pa':'ProductDeliveryPressure','product_buffer_pressure_Pa':'ProductBufferPressure',
    'plant_fresh_feed_target_Nm3_h':'PlantFreshFeedTarget','standard_temperature_K':'StandardTemperature',
    'standard_pressure_Pa':'StandardPressure','charge_feed_compression':'ChargeFeedCompression',
    'charge_recycle_compression':'ChargeRecycleCompression','compressor_efficiency':'CompressorEfficiency',
    'vacuum_efficiency':'VacuumEfficiency','compression_model':'CompressionModel','heat_capacity_ratio':'HeatCapacityRatio',
    'heating_efficiency':'HeatingEfficiency','cooling_COP':'CoolingCOP','heat_recovery_fraction':'HeatRecoveryFraction',
    'annual_operating_hours':'AnnualOperatingHours','availability':'Availability',
    'electricity_cost_per_kWh':'ElectricityCostPerKWh','heat_cost_per_kWh':'HeatCostPerKWh',
    'cooling_cost_per_kWh':'CoolingCostPerKWh','cooling_cost_basis':'CoolingCostBasis'}
CHEM_KEYS={'model':'Model','k_ads':'KAds','k_des':'KDes','Ea_ads':'EaAds','Ea_des':'EaDes',
    'T_ref':'ReferenceTemperature','heat_release_J_mol':'HeatOfAdsorption','heat_affinity_J_mol':'HeatAffinity',
    'equilibrium':'Equilibrium','qsat':'QSat','b':'Affinity','K':'HenryCoefficient','fractions':'Fractions',
    'order_ads':'OrderAds','order_des':'OrderDes','order':'Order','rate_basis':'RateBasis',
    'q_reference':'QReference','driving_regularizer':'DrivingRegularizer','n_ads':'NAds','n_des':'NDes',
    'n':'N','age_regularizer':'AgeRegularizer','history_policy':'HistoryPolicy',
    'memory_erasure_rate_ads':'MemoryErasureRateAds','memory_erasure_rate_des':'MemoryErasureRateDes'}

def mapped(data, keys):
    unknown=set(data)-keys.keys()
    if unknown: raise ValueError(f'Unmapped settings: {sorted(unknown)}')
    return {keys[k]:v for k,v in data.items()}

def chemistry(x):
    if isinstance(x,dict):return {CHEM_KEYS[k]:chemistry(v) for k,v in x.items()}
    return copy.deepcopy(x)

def sites(eq, name):
    entry=eq['components'][name];model=entry['model'];model='gab_water' if model=='gab' else model;out=[]
    if model in ('zero','none','nonadsorbing'):return out
    for old in entry['sites']:
        item=copy.deepcopy(old);heat=item.pop('heat_J_mol',item.pop('dH',0.));q=item.pop('capacity_mol_kg',item.pop('qsat',0.))
        if 'affinity_Pa_inv' in item:
            affinity=item.pop('affinity_Pa_inv')*math.exp(-heat/(R*eq.get('reference_temperature_K',298.15)))
        else:affinity=item.pop('k0')
        parameters=[q,affinity]
        if model in ('langmuir_freundlich','sips','toth','modified_redlich','freundlich'):
            nu=item.pop('nu0',item.pop('nu',1.))
            parameters.append(1./nu if model=='sips' else nu)
        if model=='freundlich':parameters=[q,parameters[2]]
        site={'Type':{'langmuir':'Langmuir','langmuir_freundlich':'Langmuir-Freundlich','sips':'Sips','toth':'Toth','freundlich':'Freundlich'}.get(model,model),'Parameters':parameters,'EquilibriumHeat':heat}
        if 'capacity_temperature' in item:
            law=item.pop('capacity_temperature');site['CapacityTemperature']={'Model':law['model'],'referenceTemperature':law['reference_temperature_K'],'Coefficient':law['coefficient']}
        if item:site['ExtraParameters']=item
        out.append(site)
    return out

def from_cycle(old):
    bed=old['bed'];thermal=old.get('thermal',{});feed=old['feed'];init=old['initial'];solver=old['solver'];out=old.get('output',{})
    doc={'SimulationType':'SwingAdsorption','SpatialDiscretization':'FVM','FluidPhase':'Gas',
        'DisplayName':old['simulation']['name'],'Components':[],
        'ColumnLength':bed['length_m'],'NumberOfGridPoints':bed['nodes']-1,
        'ParticleDensity':bed['adsorbent_density_kg_m3'],'Temperature':feed['temperature_K'],
        'InletPressure':feed['pressure_Pa'],'OutletPressure':feed['pressure_Pa'],
        'DynamicViscosity':old['gas']['viscosity_Pa_s'],'gasThermalConductivity':old['gas']['thermal_conductivity_W_m_K'],
        'Geometry':{'Type':'Monolith' if bed['type']=='monolith' else 'PackedBed',
            'InternalDiameter':bed['inner_diameter_m'],'OuterDiameter':bed['outer_diameter_m'],
            'ColumnVoidFraction':bed['void_fraction'],'ParticleDiameter':bed['particle_diameter_m']},
        'InitialConditions':{'Pressure':init['pressure_Pa'],'GasTemperature':init['gas_temperature_K'],
            'SolidTemperature':init['solid_temperature_K'],'WallTemperature':init['wall_temperature_K'],
            'GasPhaseMolFractions':init['composition'],'Loading':init['loading']},
        'MixturePredictionMethod':old['equilibrium']['method'],
        'energyBalance':thermal.get('mode','isothermal')!='isothermal',
        'heatCapacitySolid':bed['solid_heat_capacity_J_kg_K'],'heatTransferGasSolid':bed['gas_solid_heat_transfer_W_m2_K'],
        'NumberOfCycles':old['simulation'].get('cycles',1),
        'CyclicSteadyState':{'Enabled':old['simulation'].get('run_mode')=='css',
            'MaximumCycles':old['simulation'].get('maximum_cycles',200),**mapped(old.get('css',{}),CSS_KEYS)},
        'BreakthroughIntegrator':solver['method'],'ConvectionScheme':'WENO3',
        'Solver':{'RelativeTolerance':solver['rtol'],'AbsoluteTolerance':solver['atol'],
            'MaximumTimeStep':solver['maximum_step_s'],'JacobianMode':solver.get('jacobian_mode','structured_fd'),
            'CacheRHS':solver.get('cache_rhs',True),'EventTimeTolerance':solver.get('event_time_tolerance_s',1e-6), 'DispersionClosure':solver.get('dispersion_closure','concentration')},
        'WriteInterval':out.get('sample_interval_s',1.),
        'Output':{'Directory':'results','Profile':out.get('profile','compact'),'Plots':out.get('plots',True),
            'Checkpoint':out.get('checkpoint',True),'MainFormatFiles':True,'SamplesPerPhase':out.get('samples_per_phase',1000000),
            'CheckpointInterval':solver.get('checkpoint_interval_s',600.)},
        'Columns':copy.deepcopy(old.get('columns',{'A':{}})),
        'Sources':{name:process_tree(source) for name,source in old.get('sources',{}).items()},
        'SwingAdsorptionPhases':process_tree(old['cycle']['phases']),
        'Utilities':mapped(old.get('utilities',{}),UTILITY_KEYS),'Notes':old.get('notes',[])}
    for k,v in doc['Columns'].items():
        if v:
            if set(v)!={'initial'}:raise ValueError(f'Unmapped column override: {k}')
            initial_override=mapped(v['initial'],{'pressure_Pa':'Pressure','gas_temperature_K':'GasTemperature','solid_temperature_K':'SolidTemperature','wall_temperature_K':'WallTemperature','composition':'GasPhaseMolFractions','loading':'Loading'})
            doc['Columns'][k]={'InitialConditions':{**doc['InitialConditions'],**initial_override}}
    for key,public in {'wall_density_kg_m3':'wallDensity','wall_heat_capacity_J_kg_K':'heatCapacityWall','wall_conductivity_W_m_K':'wallThermalConductivity','gas_wall_heat_transfer_W_m2_K':'heatTransferGasWall','ambient_heat_transfer_W_m2_K':'heatTransferWallExternal','ambient_temperature_K':'AmbientTemperature'}.items():
        if key in thermal:doc[public]=thermal[key]
    if 'gas_storage_void_fraction' in bed:doc['GasStorageVoidFraction']=bed['gas_storage_void_fraction']
    if 'htf' in thermal:
        h=thermal['htf'];doc['Coolant']={'Model':h['model'],'Volume':h['volume_m3'],'Density':h['density_kg_m3'],'HeatCapacity':h['heat_capacity_J_kg_K'],'InitialTemperature':h['initial_temperature_K']}
        # Per-phase thermal controls contain the original inlet and flow settings.
    if 'monolith' in old:
        doc['Geometry'].update(mapped(old['monolith'],{'shape':'ChannelShape','channel_size_m':'ChannelSize','coating_thickness_m':'WashcoatThickness','channel_count':'NumberOfChannels','open_fraction':'ColumnVoidFraction','coating_solid_fraction':'WashcoatSolidFraction','substrate_density_kg_m3':'SubstrateDensity','substrate_heat_capacity_J_kg_K':'SubstrateHeatCapacity'}))
    kinetics=old.get('kinetics',{})
    for j,c in enumerate(old['components']):
        name=c['name'];rate=kinetics.get('ldf_rate_s_inv',0.);heat=kinetics.get('adsorption_heat_J_mol',0.)
        entry={'Name':name,'MolecularWeight':c['molar_mass_kg_mol'],'GasHeatCapacity':c['gas_heat_capacity_J_mol_K'],
            'AdsorbedHeatCapacity':c.get('adsorbed_heat_capacity_J_mol_K',0.),'AxialDispersionCoefficient':c.get('dispersion_m2_s',0.),
            'GasPhaseMolFraction':feed['composition'].get(name,0.) if isinstance(feed['composition'],dict) else feed['composition'][j],
            'MassTransferCoefficient':rate[j] if isinstance(rate,list) else rate,
            'HeatOfAdsorption':heat[j] if isinstance(heat,list) else heat,'PhysisorptionSites':sites(old['equilibrium'],name)}
        doc['Components'].append(entry)
    doc['PhysisorptionEnabled']=kinetics.get('physisorption_enabled',True)
    if 'chemisorption' in kinetics:doc['ChemicalKinetics']=chemistry(kinetics['chemisorption'])
    if old.get('buffers'):doc['Buffers']={k:process_tree(v) for k,v in old['buffers'].items()}
    doc['Performance']=mapped(old.get('performance',{}),{'feed_tag':'FeedTag','light_product_tag':'LightProductTag','heavy_product_tag':'HeavyProductTag','light_component':'LightComponent','heavy_component':'HeavyComponent'})
    if feed.get('velocity_basis','interstitial').lower() in ('superficial','darcy'):
        def convert_velocity_profile(value):
            if isinstance(value,(int,float)):return value/bed['void_fraction']
            if isinstance(value,list):return [convert_velocity_profile(x) for x in value]
            return {k:convert_velocity_profile(v) if k in ('initial_value','final_value','values') else v for k,v in value.items()}
        def convert_ports(value):
            if isinstance(value,dict):
                if value.get('mode')=='velocity' and 'velocity_m_s' in value:value['velocity_m_s']=convert_velocity_profile(value['velocity_m_s'])
                for child in value.values():convert_ports(child)
            elif isinstance(value,list):
                for child in value:convert_ports(child)
        phases=copy.deepcopy(old['cycle']['phases']);convert_ports(phases);doc['SwingAdsorptionPhases']=process_tree(phases)
        doc['Notes'].append('Superficial donor velocity was converted to the main interstitial basis using the bed void fraction.')
    advanced=old.get('advanced',{});labels={}
    advanced_keys={'R':'GasConstant','GAS_SOLID_THERMAL_EQUILIBRIUM':'GasSolidThermalEquilibrium','ADSORPTION_REFERENCE_TEMPERATURE':'AdsorptionReferenceTemperature','PRESSURE_DROP_MODEL':'MomentumModel','LDF_MODEL':'LDFModel','LDF_PORE_RATE_S_INV':'PoreLDFRate','WALL_END_BOUNDARY':'WallEndBoundary','WALL_END_TEMPERATURE_K':'WallEndTemperatures','KINETIC_CONCENTRATION_REGULARIZER':'KineticConcentrationRegularizer','EXPORT_VALIDATION_FLUXES':'ExportValidationFluxes','GAS_STORAGE_VOID_FRACTION':'GasStorageVoidFraction','BOUNDARY_HYDRAULIC_LENGTH':'BoundaryHydraulicLength','IAST_METHOD':'IASTMethod','WATER_PSAT_MODEL':'WaterSaturationModel','WATER_COMPONENT':'WaterComponent','RH_MODE':'RelativeHumidityMode','RH_VALUE':'RelativeHumidity'}
    for key,value in advanced.items():
        if key=='FRESH_FEED_TAGS':doc['Performance']['FreshFeedTags']=value
        elif key in advanced_keys:doc[advanced_keys[key]]=value
        elif key=='ALLOW_EXISTING_OUTPUT_DIR':pass
        elif key in ('MATERIAL','SCENARIO','SCENARIO_GROUP','CYCLE_LEVEL','P_HIGH','P_LOW','PURGE_TAG','STUDY_CONTROL_BASIS','VALIDATION_STATUS','LITERATURE_REFERENCE'):labels[key]=value
        elif key in ('CVODE_LINEAR_SOLVER','CVODE_KRYLOV_DIMENSION','CVODE_MAXIMUM_ORDER'):
            doc[{'CVODE_LINEAR_SOLVER':'CVODELinearSolver','CVODE_KRYLOV_DIMENSION':'CVODEKrylovDimension','CVODE_MAXIMUM_ORDER':'CVODEMaximumOrder'}[key]]=value
        else:raise ValueError(f'Unmapped advanced physics/control setting: {key}')
    if labels:doc['Provenance']={'OriginalLabels':labels}
    if "ChemicalKinetics" in doc:
        import sys
        sys.path.insert(0, str(ROOT))
        from ruptura.process import canonical_chemisorption
        doc = canonical_chemisorption(doc)
    return doc

def from_breakthrough(old):
    a=old['adsorbents'][0];o=old['operating'];t=old['thermal'];s=old['solver'];h=t['htf'];di=math.sqrt(4*o['cross_section_area_m2']/math.pi)
    doc={'SimulationType':'Breakthrough','SpatialDiscretization':'FVM','FluidPhase':'Gas','DisplayName':old['name'],
        'Components':[],'Geometry':{'Type':'PackedBed','InternalDiameter':di,'OuterDiameter':di*math.sqrt(1+t['wall_volume_fraction']),
            'ColumnVoidFraction':a['bed_void_fraction'],'ParticleDiameter':a['particle_diameter_m']},
        'ColumnLength':old['length_m'],'NumberOfGridPoints':old['layers'][0]['cells'],'ParticleDensity':a['solid_density_kg_m3'],
        'Temperature':o['inlet_temperature_k'],'InletPressure':o['inlet_pressure_pa'],'OutletPressure':o['outlet_pressure_pa'],
        'InletMolarFlowRate':o['inlet_molar_flow_mol_s'],'BoundaryCondition':'InletVelocityOutletPressure',
        'InitialConditions':{'Pressure':o['initial_pressure_pa'],'Temperature':o['initial_temperature_k'],
            'GasPhaseMolFractions':{c['name']:c['initial_fraction'] for c in old['components']},'Loading':'zero'},
        'Sources':{'outlet':{'Pressure':o['outlet_pressure_pa'],'Temperature':o['outlet_reservoir_temperature_k'],'GasPhaseMolFractions':o['outlet_reservoir_fractions']}},
        'MixturePredictionMethod':a['physisorption']['mixture_mode'],'DynamicViscosity':old['fluid']['viscosity_pa_s'],
        'gasThermalConductivity':old['fluid']['thermal_conductivity_w_m_k'],'energyBalance':True,
        'heatCapacitySolid':a['solid_heat_capacity_j_kg_k'],'heatTransferGasSolid':a['fluid_solid_heat_transfer_w_m2_k'],
        'wallDensity':t['wall_density_kg_m3'],'heatCapacityWall':t['wall_heat_capacity_j_kg_k'],
        'wallThermalConductivity':t['wall_thermal_conductivity_w_m_k'],'heatTransferGasWall':t['fluid_wall_heat_transfer_w_m2_k'],
        'heatTransferWallExternal':t['wall_ambient_heat_transfer_w_m2_k'],'AmbientTemperature':t['ambient_temperature_k'],
        'Coolant':{'Model':h['model'],'Volume':h['volume_m3'],'Density':h['density_kg_m3'],'HeatCapacity':h['heat_capacity_j_kg_k'],
            'InitialTemperature':h['initial_temperature_k'],'InletTemperature':h['inlet_temperature_k'],'MassFlowRate':h['mass_flow_rate_kg_s'],
            'HeatTransferCoefficient':h['heat_transfer_coefficient_w_m2_k'],'FlowDirection':h['flow_direction']},
        'BreakthroughIntegrator':s['integrator'],'ConvectionScheme':'WENO3','SimulationTime':'auto',
        'AutomaticStopping':{'Components':old['auto_end']['components'],'Threshold':old['auto_end']['threshold'],
            'MaximumTime':old['end_time_s'],'MinimumTime':old['auto_end']['minimum_time_s'],'ExtraFraction':old['auto_end']['extra_fraction'],
            'LoadingRelativeTolerance':.001},
        'Solver':{'RelativeTolerance':s['relative_tolerance'],'AbsoluteTolerance':s['absolute_tolerance'],'MaximumTimeStep':s['max_step_s'],
            'FirstTimeStep':s['first_step_s'],'JacobianMode':'structured_fd','CacheRHS':s['cache_rhs']},
        'WriteInterval':old['output']['interval_s'],'Output':{'Directory':'results','Plots':old['output']['make_plots'],'Checkpoint':old['output']['write_checkpoint'],'Profile':'compact','MainFormatFiles':True},
        'Notes':['Converted from the supplied breakthrough C++ package. Shared conservative energy equations replace the donor-specific energy formulation; compare thermal predictions before scientific reuse.',
            'Automatic stopping additionally checks loading equilibrium; this is an intentional correction to outlet-only stopping. Fixed molar feed and pure-helium outlet backflow are preserved.']}
    for c in old['components']:
        name=c['name'];eq=a['physisorption']['components'][name]
        eq={'method':doc['MixturePredictionMethod'],'components':{name:{'model':eq['model'],'sites':[{'capacity_mol_kg':q['qsat_mol_kg'],'k0':q['k0'],'heat_J_mol':q.get('dh_j_mol',0.)} for q in eq['sites']]}}}
        doc['Components'].append({'Name':name,'MolecularWeight':c['molecular_weight_kg_mol'],'GasHeatCapacity':c['gas_heat_capacity_j_mol_k'],
            'AdsorbedHeatCapacity':c['adsorbed_heat_capacity_j_mol_k'],'GasPhaseMolFraction':c['feed_fraction'],
            'CarrierGas':name=='He','AxialDispersionCoefficient':c['dispersion_m2_s'],'MassTransferCoefficient':a['physical_kinetics']['rate_coefficient'][name],
            'HeatOfAdsorption':a['heat_physisorption_j_mol'][name],'PhysisorptionSites':sites(eq,name)})
    if "ChemicalKinetics" in doc:
        import sys
        sys.path.insert(0, str(ROOT))
        from ruptura.process import canonical_chemisorption
        doc = canonical_chemisorption(doc)
    return doc

def write_case(target, doc, source, source_root):
    from generate_integrated_examples import write
    # Reuse the tested portable launcher, but put the converted case at its campaign path.
    relative=target.relative_to(ROOT/'examples/new_examples/integrated')
    source_relative=str(source.relative_to(source_root))
    doc.setdefault('Provenance',{}).update(SourceInput=source_relative,SourceSHA256=hashlib.sha256(source.read_bytes()).hexdigest())
    write(str(relative),doc,'Converted main-style input. Original source: '+source_relative+'. Read the campaign README for model correspondence and validation status.')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cycle_project',type=Path);parser.add_argument('breakthrough_project',type=Path);args=parser.parse_args()
    manifest=[];errors=[]
    cycle_root=args.cycle_project.resolve();bt_root=args.breakthrough_project.resolve()
    groups=[(cycle_root/'examples/studies/co2_ch4_mfi_sei','campaigns/co2_ch4_cycles','simulation.yaml',from_cycle),
        (bt_root/'examples/co2_ch4_mfi_sei_breakthrough','campaigns/co2_ch4_breakthrough','input.yaml',from_breakthrough),
        (cycle_root/'examples/validation/haghpanah_2014','validation/haghpanah_2014','simulation.yaml',from_cycle),
        (cycle_root/'examples/literature/mocca','validation/haghpanah_2013_candidates','simulation.yaml',from_cycle),
        (cycle_root/'examples/validation/julia_references/mocca','validation/haghpanah_2013_equations','simulation.yaml',from_cycle),
        (cycle_root/'examples/equipment','equipment','simulation.yaml',from_cycle),
        (cycle_root/'examples/chemistry','chemistry','simulation.yaml',from_cycle),
        (cycle_root/'examples/moisture','moisture','simulation.yaml',from_cycle)]
    for base,destination,pattern,convert in groups:
        for source in sorted(base.rglob(pattern)):
            if destination.endswith('haghpanah_2013_candidates') and 'haghpanah' not in str(source):continue
            target=ROOT/'examples/new_examples/integrated'/destination/source.relative_to(base).parent
            try:
                doc=convert(yaml.safe_load(source.read_text()));write_case(target,doc,source,cycle_root if convert==from_cycle else bt_root)
                manifest.append(str((target/'simulation.json').relative_to(ROOT)))
            except Exception as e:errors.append({'source':str(source.relative_to(base)),'group':destination,'reason':str(e)})
    (ROOT/'examples/new_examples/integrated/migrated_cases.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (ROOT/'verification/integration/migration.json').write_text(json.dumps({'converted':len(manifest),'errors':errors},indent=2)+'\n')
    print(json.dumps({'converted':len(manifest),'errors':errors},indent=2))
    for folder in [ROOT/'examples/new_examples/integrated/validation/haghpanah_2014',ROOT/'examples/new_examples/integrated/validation/haghpanah_2013_candidates']:
        original=cycle_root/('examples/validation/haghpanah_2014' if folder.name=='haghpanah_2014' else 'examples/literature/mocca')
        if folder.exists() and (original/'reference').exists():shutil.copytree(original/'reference',folder/'reference',dirs_exist_ok=True)

if __name__=='__main__':main()
