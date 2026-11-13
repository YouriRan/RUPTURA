"use strict";
// This editor mutates the same input object used by Lab forms and the JSON tab.
window.CycleEditor=(()=>{
  const esc=x=>String(x??"").replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const get=(o,p)=>p.reduce((a,k)=>a?.[k],o);
  const put=(o,p,v)=>{let a=o;p.slice(0,-1).forEach((k,i)=>a=a[k]??=(typeof p[i+1]==='number'?[]:{}));if(v===undefined)delete a[p.at(-1)];else a[p.at(-1)]=v;};
  const cp=x=>structuredClone(x),options=(arr,value)=>arr.map(x=>`<option value="${esc(x)}" ${String(value)===String(x)?'selected':''}>${esc(x)}</option>`).join('');
  let host,input,layout,save,fail,phase=0,selected='column:A',templates={};
  const phases=()=>input.SwingAdsorptionPhases||[],columns=()=>Object.keys(input.Columns||{A:{}}),current=()=>phases()[phase];
  function control(label,path,choices=null,type='number',fallback=''){
    const value=get(input,path)??fallback,attr=`data-ce-path='${esc(JSON.stringify(path))}'`;
    const element=choices?`<select ${attr}>${options(choices,value)}</select>`:type==='checkbox'?`<input type="checkbox" ${attr} ${value?'checked':''}>`:`<input type="${type}" step="any" ${attr} value="${esc(value)}">`;
    return `<label class="field ce-field"><span>${esc(label)}</span>${element}</label>`;
  }
  const button=(action,label,value='')=>`<button class="button" data-ce-action="${action}" data-ce-value="${esc(value)}">${label}</button>`;
  const panel=(title,body)=>`<section class="card"><div class="card-header"><h2>${title}</h2></div><div class="card-body">${body}</div></section>`;
  const composition=(prefix,key='GasPhaseMolFractions')=>`<div class="form-grid">${input.Components.map(c=>control(c.Name+(input.FluidPhase==='Liquid'?' [mol/m³]':' mole fraction'),[...prefix,key,c.Name])).join('')}</div>`;
  function nodes(){return [...columns().map((n,i)=>({id:'column:'+n,name:n,kind:'column',x:260+i*185,y:170})),...Object.keys(input.Sources||{}).map((n,i)=>({id:'source:'+n,name:n,kind:'source',x:30+(i%2)*710,y:40+Math.floor(i/2)*260})),...Object.keys(input.Buffers||{}).map((n,i)=>({id:'buffer:'+n,name:n,kind:'buffer',x:100+i*170,y:420}))].map(n=>({...n,...layout[n.id]}));}
  function endpoints(){return [...columns().flatMap(n=>[n+'.left',n+'.right']),...Object.keys(input.Sources||{}).map(n=>'source:'+n),...Object.keys(input.Buffers||{}).map(n=>'buffer:'+n)];}
  function identity(c,i){return c.Name||`${c.From} → ${c.To} (${i+1})`;}
  function valves(){const m=new Map();phases().forEach(p=>(p.Connections||[]).forEach((c,i)=>{const k=identity(c,i);if(!m.has(k))m.set(k,c);}));return m;}
  function connectionPoint(ref,ns){const [name,side]=ref.includes('.')?ref.split('.'):[];const n=ns.find(n=>n.id===(name?'column:'+name:ref));return n?{x:n.x+(n.kind==='column'?45:65),y:n.y+(n.kind==='column'?(side==='right'?0:160):32)}:null;}
  function flowsheet(){
    const ns=nodes();const p=current()||{},active=new Map((p.Connections||[]).map((c,i)=>[identity(c,i),c]));
    let links='';for(const [name,c] of valves()){
      const present=active.get(name),actual=present||c,a=connectionPoint(actual.From,ns),b=connectionPoint(actual.To,ns);if(!a||!b)continue;
      const on=present&&(present.Opening??1)>0,mx=(a.x+b.x)/2,my=(a.y+b.y)/2;
      links+=`<g data-ce-select="valve:${esc(name)}" class="ce-link ${on?'on':''}"><path d="M${a.x},${a.y} L${a.x},${my} L${b.x},${my} L${b.x},${b.y}"/><path class="ce-valve" d="M${mx-9},${my-7} L${mx+9},${my+7} L${mx+9},${my-7} L${mx-9},${my+7} Z"/><text x="${mx}" y="${my-14}" text-anchor="middle">${esc(name)}</text></g>`;
    }
    // Display legacy/direct port controls without inventing explicit valves.
    for(const col of columns())for(const side of ['Left','Right']){
      const port=p.Columns?.[col]?.[side+'Port'];if(!port||port.Mode==='closed')continue;
      const a=connectionPoint(col+'.'+side.toLowerCase(),ns),ref=port.Source?.includes(':')?port.Source:'source:'+port.Source,b=connectionPoint(ref,ns);if(a&&b)links+=`<g class="ce-link on direct"><path d="M${a.x},${a.y} L${b.x},${b.y}"/><text x="${(a.x+b.x)/2}" y="${(a.y+b.y)/2-8}">${esc(port.Mode)}</text></g>`;
    }
    const equipment=ns.map(n=>`<g tabindex="0" role="button" aria-label="Edit ${esc(n.name)}" data-ce-select="${esc(n.id)}" data-ce-drag="${esc(n.id)}" transform="translate(${n.x},${n.y})" class="ce-node ${selected===n.id?'selected':''}"><rect width="${n.kind==='column'?90:130}" height="${n.kind==='column'?160:64}" rx="${n.kind==='column'?18:8}"/><text x="${n.kind==='column'?45:65}" y="26" text-anchor="middle">${esc(n.name)}</text><text x="${n.kind==='column'?45:65}" y="48" text-anchor="middle" class="small">${n.kind==='column'?esc(p.Columns?.[n.name]?.Name||'Column'):n.kind==='buffer'?'Finite buffer':'Reservoir'}</text>${n.kind==='column'?'<circle cx="45" cy="0" r="5"/><circle cx="45" cy="160" r="5"/><text x="45" y="-12" text-anchor="middle" class="small">right / product</text><text x="45" y="183" text-anchor="middle" class="small">left / feed</text>':''}</g>`).join('');
    return `<div class="ce-canvas"><svg viewBox="0 0 940 540" aria-label="Editable cycle flowsheet">${links}${equipment}</svg></div><p class="helper">Select equipment or a valve to edit it. Drag equipment to arrange the drawing. Highlighted paths are open in the selected phase; actual flow follows the pressure difference and check-valve setting.</p>`;
  }
  function propertyPanel(){
    let [kind,...rest]=selected.split(':'),name=rest.join(':'),body='';
    if(kind==='column'&&columns().includes(name)){
      const p=['Columns',name];body=`<p>Blank values inherit the shared Column and Components settings.</p><div class="form-grid">${control('Column length [m]',[...p,'ColumnLength'])}${control('Axial intervals',[...p,'NumberOfGridPoints'])}${control('Particle density [kg/m³]',[...p,'ParticleDensity'])}${control('Initial pressure [Pa]',[...p,'InitialConditions','Pressure'])}${control('Initial loading',[...p,'InitialConditions','Loading'],['','zero','equilibrium'])}${control('Common fluid/solid temperature',[...p,'FluidSolidThermalEquilibrium'],null,'checkbox',input.FluidSolidThermalEquilibrium??false)}</div>`;
      if(columns().length>1)body+=button('delete-equipment','Remove column',selected);
    }else if(kind==='source'){
      const p=['Sources',name];body=`<div class="form-grid">${control('Pressure [Pa]',[...p,'Pressure'])}${control('Temperature [K]',[...p,'Temperature'])}</div>${composition(p,input.FluidPhase==='Liquid'?'Concentrations':'GasPhaseMolFractions')}${button('delete-equipment','Remove reservoir',selected)}`;
    }else if(kind==='buffer'){
      const p=['Buffers',name];body=`<div class="form-grid">${control('Volume [m³]',[...p,'Volume'])}${control('Initial pressure [Pa]',[...p,'InitialPressure'])}${control('Initial temperature [K]',[...p,'InitialTemperature'])}${control('Thermal mode',[...p,'ThermalMode'],['isothermal','adiabatic'])}${input.FluidPhase==='Liquid'?control('Maximum volume [m³]',[...p,'MaximumVolume']):control('Export setpoint [Pa]',[...p,'PressureControl','SetpointPressure'])+control('Controller conductance [mol/s/Pa]',[...p,'PressureControl','ControllerConductance'])+control('Export stream tag',[...p,'PressureControl','ExportTag'],null,'text')}</div>${composition(p,input.FluidPhase==='Liquid'?'InitialConcentration':'InitialComposition')}${button('delete-equipment','Remove buffer',selected)}`;
    }else if(kind==='valve'){
      const index=(current()?.Connections||[]).findIndex((c,i)=>identity(c,i)===name);
      if(index<0)body='<p>This valve is closed in this phase. Open it in the valve schedule to edit its phase settings.</p>';
      else {const p=['SwingAdsorptionPhases',phase,'Connections',index];body=`<div class="form-grid">${control('Valve name',[...p,'Name'],null,'text')}${control('From',[...p,'From'],endpoints())}${control('To',[...p,'To'],endpoints())}${control('Flow model',[...p,'FlowModel'],input.FluidPhase==='Liquid'?['linear','incompressible']:['linear','orifice'])}${control(input.FluidPhase==='Liquid'?'Conductance [m³/s/Pa]':'Conductance [mol/s/Pa]',[...p,'Conductance'])}${input.FluidPhase==='Liquid'?control('Imposed volume flow [m³/s]',[...p,'VolumetricFlowRate']):control('Maximum flow [mol/s]',[...p,'MaximumMolarFlow'])}${get(input,[...p,'FlowModel'])==='orifice'?control('Orifice area [m²]',[...p,'Area'])+control('Discharge coefficient',[...p,'DischargeCoefficient'])+control('Heat capacity ratio',[...p,'HeatCapacityRatio']):''}${control('Opening (0–1)',[...p,'Opening'])}${control('Check valve',[...p,'CheckValve'],null,'checkbox')}${control('Stream tag',[...p,'Tag'],null,'text')}${(input.FluidPhase==='Liquid'?'':`<label class="field">Compressor<input type="checkbox" data-ce-compressor="${index}" ${get(input,[...p,'Compressor'])?'checked':''}></label>`)}${input.FluidPhase==='Liquid'?'':control('Compressor efficiency',[...p,'Compressor','Efficiency'])+control('Discharge pressure [Pa]',[...p,'Compressor','OutletPressure'])+control('Imposed compressor flow [mol/s]',[...p,'Compressor','FlowRate'])}</div><p class="helper">Phase-specific valve values. Use JSON for additional equipment laws supported by the native engine.</p>${button('remove-valve','Delete this valve in all phases',name)}`;}
    }
    return panel('Equipment · '+esc(name||'select in diagram'),body||'<p>Select equipment in the flowsheet.</p>');
  }
  function schedule(){
    const rows=[...valves()].map(([name])=>`<tr><th>${esc(name)}</th>${phases().map((p,i)=>{const c=(p.Connections||[]).find((c,k)=>identity(c,k)===name);return `<td><input type="checkbox" aria-label="${esc(name)} in ${esc(p.Name)}" data-ce-valve="${esc(name)}" data-phase="${i}" ${c&&(c.Opening??1)>0?'checked':''}></td>`;}).join('')}</tr>`).join('');
    return panel('Valve schedule',`<div class="ce-table-scroll"><table class="ce-table"><thead><tr><th>Path</th>${phases().map((p,i)=>`<th><button class="button ${phase===i?'primary':''}" data-ce-action="phase" data-ce-value="${i}">${esc(p.Name||'Phase '+(i+1))}</button><small>${p.Duration??p.MaximumDuration??'event'} s</small></th>`).join('')}</tr></thead><tbody>${rows}${columns().map(col=>`<tr><th>Column ${esc(col)}</th>${phases().map(p=>`<td>${esc(p.Columns?.[col]?.Name||'Inherited ports')}</td>`).join('')}</tr>`).join('')}</tbody></table></div>${rows?'':'<p>Add a connection to create a valve schedule. Existing direct port controls appear in the phase editor.</p>'}`);
  }
  function phasePanel(){
    if(!current())return '';const p=['SwingAdsorptionPhases',phase];
    let body=`<div class="form-grid">${control('Phase name',[...p,'Name'],null,'text')}${control('Duration [s] (blank for event)',[...p,'Duration'])}${control('Minimum duration [s]',[...p,'MinimumDuration'])}${control('Maximum duration [s]',[...p,'MaximumDuration'])}${control('Timeout policy',[...p,'TimeoutPolicy'],['error','advance'])}</div>`;
    const stop=current().StopCriterion;
    body+=`<h3>Stopping condition</h3><p class="helper">A duration or a stopping condition with a maximum duration is required. Nested all/any conditions are preserved and editable in Input JSON.</p>`;
    if(stop?.criteria||stop?.Criteria||stop?.all||stop?.any||stop?.All||stop?.Any)body+='<p>Compound condition configured. Use Input JSON to edit its branches.</p>';
    else if(stop){body+=`<div class="form-grid">${control('Metric',[...p,'StopCriterion','Metric'],['bed_average_pressure','left_pressure','right_pressure','pressure_drop','outlet_mole_fraction','bed_average_temperature','cumulative_component_moles'])}${control('Column',[...p,'StopCriterion','Column'],columns())}${control('Component',[...p,'StopCriterion','Component'],['',...input.Components.map(c=>c.Name)])}${control('Port (blank selects outlet)',[...p,'StopCriterion','Port'],['','left','right'])}${control('Operator',[...p,'StopCriterion','Operator'],['>=','<='])}${control('Target',[...p,'StopCriterion','Value'])}</div>${button('remove-stop','Remove condition')}`;}
    else body+=button('add-stop','Add stopping condition');
    for(const col of columns()){
      const q=[...p,'Columns',col];body+=`<details><summary>Column ${esc(col)} · ports, temperature and readiness</summary><div class="form-grid">${control('Action name',[...q,'Name'],null,'text')}${control('Action duration [s]',[...q,'Duration'])}${control('Wait thermal mode',[...q,'WaitAction','ThermalControl','Mode'],['','adiabatic','ambient','jacket','htf'])}${control('Thermal mode',[...q,'ThermalControl','Mode'],['','adiabatic','ambient','jacket','htf','electric'])}${control('Jacket / coolant temperature [K]',[...q,'ThermalControl','Temperature'])}${control('Coolant mass flow [kg/s]',[...q,'ThermalControl','MassFlowRate'])}${control('Heat transfer coefficient [W/m²/K]',[...q,'ThermalControl','HeatTransferCoefficient'])}${control('Flow control',[...q,'FlowControl','Model'],['','fixed_velocity',...(input.FluidPhase==='Liquid'?['fixed_pressure_gradient']:[])])}</div>`;
      for(const side of ['Left','Right'])body+=`<div class="form-grid">${control(side+' port mode',[...q,side+'Port','Mode'],['','closed','pressure','velocity'])}${control(side+' source',[...q,side+'Port','Source'],['',...Object.keys(input.Sources||{})])}${control(side+' pressure [Pa]',[...q,side+'Port','Pressure'])}${control(side+' interstitial velocity [m/s]',[...q,side+'Port','Velocity'])}</div>`;
      body+='</details>';
    }
    body+=`<div class="actions">${button('duplicate-phase','Duplicate phase')}${button('phase-up','Move earlier')}${button('phase-down','Move later')}${button('delete-phase','Remove phase')}</div>`;
    return panel('Phase '+(phase+1)+' · '+esc(current().Name),body);
  }
  function cssPanel(){const p=['CyclicSteadyState'];return panel('Cyclic steady state',`<div class="form-grid">${control('Run until CSS',[...p,'Enabled'],null,'checkbox')}${control('Maximum cycles',[...p,'MaximumCycles'])}${control('Fixed cycle count',['NumberOfCycles'])}${control('Consecutive passing cycles',[...p,'ConsecutiveCycles'])}${control('State relative tolerance',[...p,'StateRelativeTolerance'])}${control('State absolute tolerance',[...p,'StateAbsoluteTolerance'])}${control('Inventory relative tolerance',[...p,'InventoryRelativeTolerance'])}${control('Performance relative tolerance',[...p,'MetricRelativeTolerance'])}${control('Mass balance tolerance',[...p,'MassBalanceTolerance'])}${control('Maximum diagnostic period',[...p,'MaximumPeriod'])}</div><p class="helper">CSS checks the complete column and finite-buffer state, inventories and cycle performance. A completed fixed-cycle run does not establish CSS.</p>`);}
  function draw(){
    phase=Math.min(phase,Math.max(phases().length-1,0));
    host.innerHTML=`<div class="ce-toolbar"><label class="field">Start from a complete example<select id="ce-template">${options(Object.keys(templates),'single-psa')}</select></label>${button('template','Load example')}${button('add-column','＋ Column')}${button('add-source','＋ Reservoir')}${button('add-buffer','＋ Buffer')}${button('add-phase','＋ Phase')}${button('cluster','↧ Cluster bundle')}</div>${panel('Cycle flowsheet',`<div class="ce-phase-label">Selected phase <strong>${esc(current()?.Name||'—')}</strong></div>${flowsheet()}<div class="ce-connect"><label class="field">From<select id="ce-from">${options(endpoints(),endpoints()[0])}</select></label><label class="field">To<select id="ce-to">${options(endpoints(),endpoints()[1])}</select></label>${button('connect','＋ Connect with valve')}</div>`)}${schedule()}<div class="two-column"><div>${phasePanel()}${cssPanel()}</div><div>${propertyPanel()}${panel('One simulation configuration','<p>Cycle controls and the Input JSON tab edit the same simulation. Validate before running. The cluster bundle contains that JSON and portable launch scripts.</p><p>Sources are infinite reservoirs. Buffers store finite inventory and affect CSS and recovery.</p>')}</div></div>`;
  }
  const changed=()=>{save();};
  function nextName(prefix,used){let i=1;while(used.includes(prefix+i))i++;return prefix+i;}
  async function action(name,value){
    const p=current();
    if(name==='phase'){phase=Number(value);draw();return;}
    if(name==='template'){
      const d=cp(templates[host.querySelector('#ce-template').value]);Object.keys(input).forEach(k=>delete input[k]);Object.assign(input,d);Object.keys(layout).forEach(k=>delete layout[k]);phase=0;selected='column:A';
    }
    if(name==='add-column'){
      if(columns().length>=4)throw Error('The cycle engine supports up to four columns.');
      input.Columns??={A:{}};const n='ABCD'.split('').find(n=>!input.Columns[n]);input.Columns[n]={};phases().forEach(p=>{p.Columns??={};p.Columns[n]={LeftPort:{Mode:'closed'},RightPort:{Mode:'closed'}};});selected='column:'+n;
    }
    if(name==='add-source'||name==='add-buffer'){
      const buffer=name==='add-buffer',key=buffer?'Buffers':'Sources';input[key]??={};const n=nextName(buffer?'tank':'reservoir',Object.keys(input[key]));
      const composition=Object.fromEntries(input.Components.map((c,i)=>[c.Name,input.FluidPhase==='Liquid'?0:(i===input.Components.length-1?1:0)]));
      input[key][n]=buffer?{Volume:.005,InitialPressure:1e5,InitialTemperature:input.Temperature??298.15,[input.FluidPhase==='Liquid'?'InitialConcentration':'InitialComposition']:composition,ThermalMode:'isothermal'}:{Pressure:1e5,Temperature:input.Temperature??298.15,[input.FluidPhase==='Liquid'?'Concentrations':'GasPhaseMolFractions']:composition};selected=(buffer?'buffer:':'source:')+n;
    }
    if(name==='add-phase'){input.SwingAdsorptionPhases??=[];input.SwingAdsorptionPhases.push({Name:'Phase '+(phases().length+1),Duration:10,Columns:Object.fromEntries(columns().map(n=>[n,{LeftPort:{Mode:'closed'},RightPort:{Mode:'closed'}}])),Connections:[]});phase=phases().length-1;}
    if(name==='duplicate-phase'){phases().splice(phase+1,0,cp(p));phase++;current().Name+=' copy';}
    if(name==='delete-phase'){if(phases().length<2)throw Error('Keep at least one phase.');phases().splice(phase,1);}
    if(name==='phase-up'&&phase>0){[phases()[phase-1],phases()[phase]]=[p,phases()[phase-1]];phase--;}
    if(name==='phase-down'&&phase<phases().length-1){[phases()[phase+1],phases()[phase]]=[p,phases()[phase+1]];phase++;}
    if(name==='connect'){
      const from=host.querySelector('#ce-from').value,to=host.querySelector('#ce-to').value;if(from===to)throw Error('Choose different endpoints.');
      const n=nextName('V-',Array.from(valves().keys()));p.Connections??=[];p.Connections.push({Name:n,From:from,To:to,FlowModel:'linear',Conductance:input.FluidPhase==='Liquid'?1e-10:1e-8,Opening:1,CheckValve:true,Tag:'transfer'});selected='valve:'+n;
    }
    if(name==='remove-valve')phases().forEach(p=>p.Connections=(p.Connections||[]).filter((c,i)=>identity(c,i)!==value));
    if(name==='delete-equipment'){
      const [kind,n]=value.split(':'),key={column:'Columns',source:'Sources',buffer:'Buffers'}[kind];delete input[key][n];delete layout[value];
      phases().forEach(p=>{if(kind==='column')delete p.Columns?.[n];p.Connections=(p.Connections||[]).filter(c=>!([c.From,c.To].some(e=>kind==='column'?e.startsWith(n+'.'):e===value)));for(const a of Object.values(p.Columns||{}))for(const s of ['LeftPort','RightPort'])if(kind==='source'&&a[s]?.Source===n)a[s]={Mode:'closed'};});selected='column:'+columns()[0];
    }
    if(name==='add-stop'){delete p.Duration;p.MaximumDuration??=100;p.StopCriterion={Metric:'bed_average_pressure',Column:columns()[0],Operator:'>=',Value:input.InletPressure??5e5};}
    if(name==='remove-stop'){delete p.StopCriterion;p.Duration??=p.MaximumDuration??10;delete p.MaximumDuration;delete p.MinimumDuration;}
    if(name==='cluster'){
      const r=await fetch('/api/export-cluster',{method:'POST',headers:{'Content-Type':'application/json','X-Ruptura-Session':window.rupturaSession()},body:JSON.stringify({input})});
      if(!r.ok)throw Error((await r.json()).error);const url=URL.createObjectURL(await r.blob()),a=document.createElement('a');a.href=url;a.download='ruptura-cluster.zip';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);return;
    }
    changed();draw();
  }
  function bind(){
    host.onclick=e=>{const b=e.target.closest('[data-ce-action]'),node=e.target.closest('[data-ce-select]');if(b){e.stopPropagation();action(b.dataset.ceAction,b.dataset.ceValue).catch(fail);}else if(node){selected=node.dataset.ceSelect;draw();}};
    host.onchange=e=>{try{const t=e.target;if(t.dataset.cePath){
      const path=JSON.parse(t.dataset.cePath),old=get(input,path),v=t.type==='checkbox'?t.checked:t.type==='number'?(t.value===''?undefined:Number(t.value)):(t.value===''?undefined:t.value);if(typeof v==='number'&&!Number.isFinite(v))throw Error('Enter a finite number.');
      const rename=path.includes('Connections')&&path.at(-1)==='Name'&&old!==v,oldIdentity=selected.slice('valve:'.length);
      if(rename&&(!v||valves().has(v)))throw Error('Choose a nonempty, unique valve name.');
      put(input,path,v);
      if(rename){phases().forEach(p=>(p.Connections||[]).forEach((c,i)=>{if(identity(c,i)===oldIdentity)c.Name=v;}));selected='valve:'+v;}
      // A blank flow-control selection means inheritance, not an empty model.
      if(path.at(-1)==='Model'&&path.includes('FlowControl')&&v===undefined)put(input,path.slice(0,-1),undefined);
      if(path.at(-1)==='Mode'&&v===undefined)put(input,path.slice(0,-1),undefined);
      changed();if(t.tagName==='SELECT'||t.type==='checkbox'||rename)draw();
    }else if(t.dataset.ceCompressor!==undefined){const c=current().Connections[Number(t.dataset.ceCompressor)];if(t.checked)c.Compressor={Efficiency:.75,OutletPressure:input.InletPressure??5e5};else delete c.Compressor;changed();draw();}else if(t.dataset.ceValve){const i=Number(t.dataset.phase),p=phases()[i],name=t.dataset.ceValve;p.Connections??=[];const c=p.Connections.find((c,k)=>identity(c,k)===name);if(c)c.Opening=t.checked?1:0;else if(t.checked)p.Connections.push({...cp(valves().get(name)),Opening:1});changed();draw();}}catch(error){fail(error);}e.stopPropagation();};
    let drag=null;
    host.onpointerdown=e=>{const n=e.target.closest('[data-ce-drag]');if(!n)return;const svg=n.ownerSVGElement,p=new DOMPoint(e.clientX,e.clientY).matrixTransform(svg.getScreenCTM().inverse()),item=nodes().find(x=>x.id===n.dataset.ceDrag);drag={id:item.id,x:item.x,y:item.y,start:p,node:n,svg};host.setPointerCapture(e.pointerId);};
    host.onpointermove=e=>{if(!drag)return;const p=new DOMPoint(e.clientX,e.clientY).matrixTransform(drag.svg.getScreenCTM().inverse());layout[drag.id]={x:Math.max(0,Math.min(800,drag.x+p.x-drag.start.x)),y:Math.max(25,Math.min(350,drag.y+p.y-drag.start.y))};drag.node.setAttribute('transform',`translate(${layout[drag.id].x},${layout[drag.id].y})`);};
    host.onpointerup=()=>{if(drag){selected=drag.id;drag=null;changed();draw();}};
    host.onkeydown=e=>{const n=e.target.closest('[data-ce-select]');if(n&&['Enter',' '].includes(e.key)){e.preventDefault();selected=n.dataset.ceSelect;draw();}};
  }
  async function mount(node,config,positions,onSave,onError){host=node;input=config;layout=positions;save=onSave;fail=onError;templates=await(await fetch('/api/cycle-templates',{headers:{'X-Ruptura-Session':window.rupturaSession()}})).json();draw();bind();}
  return {mount};
})();
