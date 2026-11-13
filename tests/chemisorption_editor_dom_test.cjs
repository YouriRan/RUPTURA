/* Execute the actual Lab site editor and field handlers in a DOM. */
const assert=require('node:assert/strict');
const fs=require('node:fs');const path=require('node:path');const {execFileSync}=require('node:child_process');
const {JSDOM}=require('jsdom');const root=path.resolve(__dirname,'..');
const schema=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c','import json;from ruptura.webapp.config import SCHEMA;print(json.dumps(SCHEMA))'],{cwd:root,env:{...process.env,PYTHONPATH:root},encoding:'utf8'}));
const source=fs.readFileSync(path.join(root,'ruptura/webapp/static/app.js'),'utf8');
const dom=new JSDOM(fs.readFileSync(path.join(root,'ruptura/webapp/static/index.html'),'utf8'),{runScripts:'outside-only',url:'http://localhost'});
const w=dom.window;w.structuredClone=structuredClone;
w.fetch=async()=>({ok:true,json:async()=>({})});
w.eval(source.slice(0,source.lastIndexOf('boot().catch'))+'\nwindow.testLab={state,config,renderComponents,changeField};');
const api=w.testLab;
const input={SimulationType:'Breakthrough',FluidPhase:'Gas',Components:[{Name:'A',GasPhaseMolFraction:1.,PhysisorptionSites:[],ChemisorptionSites:[{Type:'Langmuir',Parameters:[1.,1e-5],MaximumLoading:1.,RateEquation:'FirstOrder',RateEquationParameters:{MassTransferCoefficient:.04}}]}]};
Object.assign(api.state,{schema,workspace:{selectedId:'one',simulations:[{id:'one',input}]},activeTab:'components',page:'simulation',capabilities:{sundials:true}});
const target=(tail)=>[...w.document.querySelectorAll('[data-path]')].find(n=>JSON.parse(n.dataset.path).slice(-tail.length).join('/')===tail.join('/')&&['INPUT','SELECT'].includes(n.tagName));
function edit(tail,value){const n=target(tail);assert.ok(n,`Missing field ${tail}`);n.value=String(value);api.changeField(n);}
api.renderComponents();edit(['RateEquation'],'DualNth');
let params=input.Components[0].ChemisorptionSites[0].RateEquationParameters;
assert.equal(params.Branches.length,2);assert.equal(input.SpatialDiscretization,'FDM');
edit(['Branches',1,'DesorptionOrder'],1.7);assert.equal(params.Branches[1].DesorptionOrder,1.7);
edit(['Branches',0,'Fraction'],.3);edit(['Branches',1,'Fraction'],.7);
assert.equal(params.Branches[0].Fraction+params.Branches[1].Fraction,1.);
const before=JSON.stringify(input);api.renderComponents();assert.equal(JSON.stringify(input),before);
input.FluidPhase='Liquid';api.renderComponents();assert.ok(w.document.querySelector('#content').textContent.includes('m³/mol'));
edit(['RateEquation'],'Avrami');edit(['HistoryPolicy'],'ContinuousErasure');edit(['MemoryErasureRateAds'],.12);edit(['HistoryDrivingRegularizer'],1e-7);
params=input.Components[0].ChemisorptionSites[0].RateEquationParameters;
assert.equal(params.HistoryPolicy,'ContinuousErasure');assert.equal(params.MemoryErasureRateAds,.12);assert.equal(params.HistoryDrivingRegularizer,1e-7);
edit(['RateEquation'],'General');assert.equal(target(['HistoryPolicy']),undefined);
const report={status:'passed',scope:'DOM component checks; full browser layout not covered',checks:['DualNth initialization selects the shared column','Two branch forms edit the canonical JSON','Branch fractions and independent desorption orders','Rendering preserves branch state','Liquid affinity units','Avrami history and erasure controls','Explicit history-switching scale','General kinetics retain their own controls']};
if(process.argv[2]){fs.mkdirSync(path.dirname(process.argv[2]),{recursive:true});fs.writeFileSync(process.argv[2],JSON.stringify(report,null,2)+'\n');}
console.log(JSON.stringify(report,null,2));w.close();
