/* Cycle-editor component tests, without a browser or simulation server.
 * npm install --prefix .build-tools/dom-tests jsdom
 * NODE_PATH="$PWD/.build-tools/dom-tests/node_modules" node tests/cycle_editor_dom_test.cjs
 * Full browser layout/accessibility checks remain a separate validation step.
 */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const {JSDOM} = require('jsdom');
const root = path.resolve(__dirname, '..');
const templates = JSON.parse(execFileSync(process.env.PYTHON || 'python3', ['-c',
  'import json;from ruptura.webapp.cycle_templates import templates;print(json.dumps(templates()))'],
  {cwd:root, env:{...process.env,PYTHONPATH:root}, encoding:'utf8'}));
const source = fs.readFileSync(path.join(root,'ruptura/webapp/static/cycle_editor.js'),'utf8');
const checks = [];
const tick = () => new Promise(resolve => setImmediate(resolve));
async function editor(config) {
  const dom = new JSDOM('<div id="editor"></div>', {runScripts:'outside-only',url:'http://localhost'});
  const w = dom.window, host = w.document.querySelector('#editor');
  let saves=0, exported=null;const errors=[];
  w.structuredClone=structuredClone;
  w.rupturaSession=()=>'';
  w.URL.createObjectURL=()=> 'blob:cycle-test';w.URL.revokeObjectURL=()=>{};
  w.HTMLAnchorElement.prototype.click=function(){};
  w.fetch=async (url,options)=>url==='/api/cycle-templates'?{json:async()=>templates}:
    (exported=JSON.parse(options.body).input,{ok:true,blob:async()=>new w.Blob(['test'])});
  w.eval(source);
  await w.CycleEditor.mount(host,config,{},()=>{saves++;},e=>errors.push(e.message));
  const field=label=>[...host.querySelectorAll('.ce-field')].find(n=>n.querySelector('span').textContent===label)?.querySelector('input,select');
  const select=async ref=>{host.querySelector(`[data-ce-select="${ref}"]`).dispatchEvent(new w.MouseEvent('click',{bubbles:true}));await tick();};
  const action=async name=>{host.querySelector(`[data-ce-action="${name}"]`).click();await tick();};
  const edit=async(label,value)=>{const node=field(label);assert.ok(node,`Missing field ${label}`);node.value=String(value);node.dispatchEvent(new w.Event('change',{bubbles:true}));await tick();};
  return {w,host,field,select,action,edit,config,errors,saves:()=>saves,exported:()=>exported};
}
(async()=>{
  const d=structuredClone(templates['two-psa-pe']),before=JSON.stringify(d),e=await editor(d);
  assert.equal(JSON.stringify(d),before);checks.push('mount preserves the canonical simulation');
  await e.select('column:B');await e.edit('Column length [m]',.2);
  assert.equal(d.Columns.B.ColumnLength,.2);assert.equal(d.ColumnLength,JSON.parse(before).ColumnLength);
  await e.select('source:feed');await e.edit('Pressure [Pa]',550000);
  assert.equal(d.Sources.feed.Pressure,550000);checks.push('column overrides and reservoir edits');
  e.host.querySelector('[data-ce-action="phase"][data-ce-value="2"]').click();await tick();
  assert.ok(e.host.querySelector('[data-ce-select="valve:V-E"].on'));
  const toggle=e.host.querySelector('[data-ce-valve="V-E"][data-phase="2"]');
  toggle.checked=false;toggle.dispatchEvent(new e.w.Event('change',{bubbles:true}));await tick();
  assert.equal(d.SwingAdsorptionPhases[2].Connections[0].Opening,0);
  assert.equal(e.host.querySelector('[data-ce-select="valve:V-E"].on'),null);
  checks.push('valve schedule and pressure-equalization highlight');
  await e.select('valve:V-E');await e.edit('Valve name','PE-shared');
  assert.equal(d.SwingAdsorptionPhases[2].Connections[0].Name,'PE-shared');
  assert.equal(d.SwingAdsorptionPhases[5].Connections[0].Name,'PE-shared');
  checks.push('valve identity is preserved across phases when renamed');
  await e.action('add-buffer');
  assert.ok(d.Buffers.tank1.InitialComposition);
  await e.edit('Export setpoint [Pa]',100000);
  await e.edit('Controller conductance [mol/s/Pa]',1e-8);
  assert.equal(d.Buffers.tank1.PressureControl.ControllerConductance,1e-8);
  checks.push('gas buffer and controller editing');
  e.host.querySelector('#ce-from').value='A.right';e.host.querySelector('#ce-to').value='buffer:tank1';
  await e.action('connect');
  assert.ok(d.SwingAdsorptionPhases[2].Connections.some(c=>c.From==='A.right'&&c.To==='buffer:tank1'));
  const compressor=e.host.querySelector('[data-ce-compressor]');
  compressor.checked=true;compressor.dispatchEvent(new e.w.Event('change',{bubbles:true}));await tick();
  await e.edit('Compressor efficiency',.8);
  assert.equal(d.SwingAdsorptionPhases[2].Connections.at(-1).Compressor.Efficiency,.8);
  checks.push('connection and compressor editing');
  await e.action('add-column');assert.ok(d.Columns.C);
  await e.action('delete-equipment');assert.equal(d.Columns.C,undefined);
  await e.action('duplicate-phase');assert.equal(d.SwingAdsorptionPhases.length,7);
  await e.action('phase-up');await e.action('phase-down');await e.action('delete-phase');
  assert.equal(d.SwingAdsorptionPhases.length,6);checks.push('column and phase lifecycle');
  await e.edit('Maximum cycles',22);assert.equal(d.CyclicSteadyState.MaximumCycles,22);
  await e.action('cluster');assert.equal(JSON.stringify(e.exported()),JSON.stringify(d));
  const copy=JSON.parse(JSON.stringify(d)),other=await editor(copy);
  assert.equal(JSON.stringify(copy),JSON.stringify(d));
  assert.ok(other.host.querySelector('[data-ce-select="valve:PE-shared"]'));
  checks.push('CSS controls, JSON round-trip and exact cluster request');
  assert.deepEqual(e.errors,[]);assert.deepEqual(other.errors,[]);assert.ok(e.saves()>10);
  const liquid=JSON.parse(fs.readFileSync(path.join(root,'ruptura/webapp/presets/fvm-liquid-cycle.json'),'utf8'));
  const le=await editor(liquid);await le.action('add-buffer');
  assert.ok(liquid.Buffers.tank1.InitialConcentration);
  assert.equal(liquid.Buffers.tank1.InitialComposition,undefined);
  assert.ok(le.field('Maximum volume [m³]'));
  assert.equal(le.field('Export setpoint [Pa]'),undefined);
  await le.edit('A [mol/m³]',.25);
  assert.equal(liquid.Buffers.tank1.InitialConcentration.A,.25);assert.deepEqual(le.errors,[]);
  checks.push('liquid buffers use concentration and volume controls');
  liquid.SwingAdsorptionPhases[0].StopCriterion={All:[{Metric:'bed_average_pressure',Column:'A',Operator:'>=',Value:1e5}]};
  delete liquid.SwingAdsorptionPhases[0].Duration;liquid.SwingAdsorptionPhases[0].MaximumDuration=1;
  const compound=await editor(liquid);
  assert.ok(compound.host.textContent.includes('Compound condition configured'));
  assert.equal(compound.field('Metric'),undefined);checks.push('compound stop conditions are preserved');
  const report={status:'passed',scope:'DOM component behavior; full browser rendering not covered',checks};
  if(process.argv[2]){fs.mkdirSync(path.dirname(process.argv[2]),{recursive:true});fs.writeFileSync(process.argv[2],JSON.stringify(report,null,2)+'\n');}
  console.log(JSON.stringify(report,null,2));
  e.w.close();other.w.close();le.w.close();compound.w.close();
})().catch(error=>{console.error(error);process.exitCode=1;});
