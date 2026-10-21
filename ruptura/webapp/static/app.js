"use strict";

const $ = (query, root = document) => root.querySelector(query);
const $$ = (query, root = document) => [...root.querySelectorAll(query)];
const clone = value => structuredClone(value);
const escapeHtml = value => String(value ?? "").replace(/[&<>"']/g, char => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[char]));
const uid = () => crypto.randomUUID();
const format = value => value == null ? "—" : Number(value).toLocaleString(undefined, {maximumSignificantDigits: 5});
const state = {workspace: null, schema: null, presets: null, runs: [], activeTab: "setup", session: "", capabilities: {},
  page: "about", theme: "dark", selectedRun: null, metadata: null, runInfo: null, plotOptions: {}, chartRequest: 0, resultRequest: 0, playing: null};
let saveTimer, saveQueue = Promise.resolve(), toastTimer, logOffset = null, polling = false;
const simulation = () => state.workspace.simulations.find(item => item.id === state.workspace.selectedId) || state.workspace.simulations[0];
const config = () => simulation().input;
const activeRun = () => state.runs.find(item => item.id === state.selectedRun);
const ongoing = record => record && ["queued", "running", "canceling"].includes(record.status);

async function api(path, body, method = body === undefined ? "GET" : "POST") {
  const response = await fetch(path, {method, headers: {"Content-Type":"application/json", "X-Ruptura-Session": state.session},
    ...(body === undefined ? {} : {body: JSON.stringify(body)})});
  const value = await response.json();
  if (!response.ok) throw new Error(value.error || `Request failed (${response.status})`);
  return value;
}
function toast(message, error = false) {
  clearTimeout(toastTimer);
  const node = $("#toast"); node.textContent = message; node.className = `toast${error ? " error" : ""}`; node.hidden = false;
  toastTimer = setTimeout(() => { node.hidden = true; }, error ? 8000 : 3500);
}
function notice(message = "", error = false) {
  const node = $("#notice"); node.textContent = message; node.className = `notice${error ? " error" : ""}`; node.hidden = !message;
}
function fail(error) { console.error(error); notice(error.message, true); toast(error.message, true); }
function scheduleSave() {
  $("#save-status").textContent = "Unsaved changes";
  clearTimeout(saveTimer); saveTimer = setTimeout(() => saveNow().catch(fail), 500);
}
function saveNow() {
  clearTimeout(saveTimer);
  const snapshot = clone(state.workspace);
  $("#save-status").textContent = "Saving…";
  saveQueue = saveQueue.catch(() => {}).then(() => api("/api/state", snapshot)).then(() => { $("#save-status").textContent = "All changes saved"; });
  return saveQueue;
}
function getPath(object, path) { return path.reduce((value, key) => value == null ? undefined : value[key], object); }
function setPath(object, path, value) {
  let target = object;
  path.slice(0, -1).forEach((key, i) => { if (target[key] == null) target[key] = typeof path[i + 1] === "number" ? [] : {}; target = target[key]; });
  if (value === undefined) delete target[path.at(-1)]; else target[path.at(-1)] = value;
}
const pathAttr = path => `data-path="${escapeHtml(JSON.stringify(path))}"`;
const describe = (key, label, value = null, options = null, type = "number") => ({key, label, default: value, options, type: options ? "select" : type});
function helpButton(label, description) {
  return `<button type="button" class="help-button" data-help="${escapeHtml(description)}" aria-label="About ${escapeHtml(label)}">?</button>`;
}
function fieldHeader(id, label, description, reset = "") {
  return `<div class="field-label-row"><label for="${id}">${escapeHtml(label)}</label><span class="field-actions">${reset}${helpButton(label, description)}</span></div>`;
}
function fieldDescription(definition, path) {
  if (definition.description) return definition.description;
  const key = path.at(-1);
  if (key === "Name") {
    if (path[0] === "Components") return "Set the component name used in result files, reaction participants, and adsorbent overrides.";
    if (path[0] === "Adsorbents") return "Set the adsorbent name referenced by the column's bed sections.";
    if (path[0] === "SwingAdsorptionPhases") return "Set a descriptive name for this operating phase of the swing process.";
  }
  if (key === "Type" && path[0] === "Geometry") return "Choose a packed bed or a monolith with parallel flow channels.";
  if (Number.isInteger(key) && path.at(-2) === "Parameters") {
    const site = getPath(config(), path.slice(0,-2));
    const parameter = state.schema.isotherms[site.Type]?.[key]?.[0];
    return state.schema.help[parameter] || "Set this parameter of the selected equilibrium isotherm.";
  }
  if (path[0] === "Reactions" && Number.isInteger(key)) {
    if (path[2] === "Reactants" || path[2] === "Products") return "Choose the component taking part on this side of the reaction.";
    if (path[2] === "Stoichiometry") return "Set the positive stoichiometric coefficient specifying how much of this component is consumed or produced.";
    if (path[2] === "ForwardOrders" || path[2] === "BackwardOrders") return "Set this participant's rate-law exponent, which defaults to its stoichiometric coefficient.";
  }
  return state.schema.help[key] || `Set ${definition.label.toLowerCase()} for the selected simulation.`;
}
function field(definition, prefix = [], extra = "") {
  const path = [...prefix, definition.key], id = `field-${path.join("-")}`;
  let value = getPath(config(), path);
  if (path.length === 1 && path[0] === "ColumnLength" && config().ColumnSections?.length) {value=columnLength();extra+=" disabled";}
  if (path[0] === "Reactions" && ["Reactants","Products"].includes(path[2]) && typeof value === "string") value=config().Components.findIndex(component=>component.Name.toLowerCase()===value.toLowerCase());
  const label = escapeHtml(definition.label), optional = !Number.isInteger(path.at(-1));
  const keepDefault = !optional || (definition.options && definition.default != null);
  const reset = value !== undefined && (optional || definition.default != null) ? `<button type="button" class="reset" data-action="reset" ${pathAttr(path)} ${keepDefault ? `data-default="${escapeHtml(JSON.stringify(definition.default))}"` : ""} title="Use default value" aria-label="Reset ${label}">↺</button>` : "";
  const description = fieldDescription(definition, path);
  if (definition.type === "boolean") return `<div class="field checkbox-field"><input id="${id}" type="checkbox" ${pathAttr(path)} ${value ? "checked" : ""} ${extra}><label for="${id}">${label}</label><span class="field-actions">${reset}${helpButton(definition.label, description)}</span></div>`;
  let control;
  if (definition.options) {
    const options = definition.options.filter(item=>definition.key!=="BreakthroughIntegrator"||item!=="CVODE"||state.capabilities.sundials).map(item => typeof item === "object" ? item : {value: item, label: item});
    control = `<select id="${id}" ${pathAttr(path)} ${definition.default == null ? "" : `data-default="${escapeHtml(JSON.stringify(definition.default))}"`} ${extra}><option value="">Default${definition.default == null ? "" : `: ${escapeHtml(definition.default)}`}</option>${options.map(item => `<option value="${escapeHtml(item.value)}" ${String(value) === String(item.value) ? "selected" : ""}>${escapeHtml(item.label)}</option>`).join("")}</select>`;
  } else {
    control = `<input id="${id}" type="${definition.type === "text" ? "text" : "number"}" step="any" ${pathAttr(path)} value="${escapeHtml(value ?? "")}" placeholder="${definition.default == null ? "Optional" : `Default: ${escapeHtml(definition.default)}`}" ${extra}>`;
  }
  return `<div class="field">${fieldHeader(id, definition.label, description, reset)}${control}</div>`;
}
function fields(group, prefix = [], filter = () => true) { return `<div class="form-grid">${state.schema[group].filter(filter).map(item => field(item, prefix)).join("")}</div>`; }
function card(title, content, subtitle = "", action = "") {
  return `<section class="card"><div class="card-header"><div><h2>${escapeHtml(title)}</h2>${subtitle ? `<p>${escapeHtml(subtitle)}</p>` : ""}</div>${action}</div><div class="card-body">${content}</div></section>`;
}
function jsonEditor(value, name, compact = true, label = "Apply changes") {
  const title = name === "simulation" ? "Simulation input JSON" : name === "grid" ? "Grid positions JSON" : "Component overrides JSON";
  const description = name === "simulation" ? "Edit the complete simulation input as JSON and validate it before applying changes." : name === "grid" ? "Enter axial positions in metres as a JSON array, or use an empty array for automatic spacing." : "Enter component-specific adsorption and kinetic overrides as a JSON object.";
  const reset = `<button type="button" class="reset" data-action="reset-json" data-editor="${escapeHtml(name)}" title="Discard unapplied changes" aria-label="Reset ${title}">↺</button>`;
  return `${fieldHeader(`json-${name}`, title, description, reset)}<textarea class="json-editor${compact ? " compact" : ""}" id="json-${name}" spellcheck="false" aria-label="${escapeHtml(name)} JSON">${escapeHtml(JSON.stringify(value, null, 2))}</textarea><div class="editor-actions"><span class="helper" id="json-status-${name}">Changes take effect after applying.</span><button class="button" data-action="apply-json" data-editor="${escapeHtml(name)}">${label}</button></div>`;
}
function pill(status) { return `<span class="pill ${escapeHtml(status)}">${escapeHtml(status)}</span>`; }
function stopPlayback() { clearInterval(state.playing); state.playing = null; const button=$("#play-button"); if(button) button.textContent="▶"; }
function setPage(page) {
  state.page = page;
  closeSimulationMenu();
  $("#simulation-page").hidden = page !== "simulation";
  $("#about-page").hidden = page !== "about";
  $("#workspace-button").setAttribute("aria-pressed",String(page === "simulation"));
  $("#about-button").setAttribute("aria-pressed", String(page === "about"));
  $("#breadcrumb-name").textContent = page === "about" ? "About" : config().DisplayName || "Untitled simulation";
}
function selectTab(tab) { stopPlayback(); setPage("simulation"); state.activeTab = tab; header(); renderContent(); $$("[data-tab]").forEach(button => button.classList.toggle("active", button.dataset.tab === tab)); }
function header() {
  closeSimulationMenu();
  const input = config();
  $("#project-title").textContent = input.DisplayName || "Untitled simulation";
  $("#breadcrumb-name").textContent = state.page === "about" ? "About" : input.DisplayName || "Untitled simulation";
  $("#project-description").textContent = input.SimulationType === "MixturePrediction" ? "Explore equilibrium, mixture loading, and selectivity." : `${input.FluidPhase === "Liquid" ? "Liquid-phase" : "Gas-phase"} ${input.SimulationType === "SwingAdsorption" ? "swing adsorption" : "breakthrough"} · Configure, simulate, and explore your results.`;
  const list=$("#experiments"),scroll=list.scrollTop,left=list.scrollLeft;
  list.innerHTML=state.workspace.simulations.map(item=>`<div class="experiment-row ${item.id===simulation().id?"active":""}"><button class="experiment ${item.id===simulation().id?"active":""}" data-action="select-simulation" data-id="${escapeHtml(item.id)}"><strong>${escapeHtml(item.input.DisplayName||"Untitled simulation")}</strong><small>${escapeHtml(item.input.SimulationType||"Breakthrough")}</small></button><button class="simulation-menu-button" data-action="simulation-menu" data-id="${escapeHtml(item.id)}" aria-label="Actions for ${escapeHtml(item.input.DisplayName||"simulation")}" aria-haspopup="menu" aria-expanded="false" aria-controls="simulation-menu" title="Simulation actions">⋯</button></div>`).join("");
  list.scrollTop=scroll;list.scrollLeft=left;
  $("#experiment-count").textContent=String(state.workspace.simulations.length).padStart(2,"0");
  $("#component-count").textContent = input.Components?.length || 0;
  const runs = state.runs.filter(record => record.simulationId === simulation().id || record.id === simulation().lastRunId);
  $("#results-count").textContent = runs.length || "";
  const latest = runs[0], mixture = input.SimulationType === "MixturePrediction";
  const metrics = [
    ["Simulation type", mixture ? "Equilibrium" : input.SimulationType === "SwingAdsorption" ? "Swing adsorption" : "Breakthrough", ""],
    ["Temperature", format(input.Temperature ?? 433), "K"],
    [mixture ? "Mixture method" : "Column length", mixture ? (input.MixturePredictionMethod || "IAST") : format(columnLength()), mixture ? "" : "m"],
    ["Latest run", latest ? pill(latest.status) : "Ready to run", ""]];
  $("#metrics").innerHTML = metrics.map(([label, value, unit]) => `<div class="metric"><span class="metric-label">${label}</span><span class="metric-value">${label === "Latest run" && latest ? value : escapeHtml(value)}</span><span class="metric-unit">${unit}</span></div>`).join("");
  $("#engine-dot").className = `status-dot ${state.capabilities.cli && !state.capabilities.error ? "green" : "orange"}`;
}
function render() { setPage("simulation"); header(); renderContent(); }
function renderContent() {
  if (state.activeTab === "setup") renderSetup();
  else if (state.activeTab === "components") renderComponents();
  else if (state.activeTab === "column") renderColumn();
  else if (state.activeTab === "reactions") renderReactions();
  else if (state.activeTab === "notes") renderNotes();
  else if (state.activeTab === "json") renderJson();
  else if (state.activeTab === "results") renderResults().catch(fail);
  else if (state.activeTab === "log") renderLog().catch(fail);
}
function renderAbout() {
  stopPlayback(); setPage("about");
  const about = state.about,item=about.citation;
  const reference = `<article class="reference"><a class="reference-title" href="https://doi.org/${escapeHtml(item.doi)}" target="_blank" rel="noopener noreferrer">${escapeHtml(item.title)} ↗</a><p>${escapeHtml(item.authors)}</p><p>${escapeHtml(item.publication)}</p><a class="reference-doi" href="https://doi.org/${escapeHtml(item.doi)}" target="_blank" rel="noopener noreferrer">doi:${escapeHtml(item.doi)}</a></article>`;
  $("#about-page").innerHTML = `<div class="about-heading"><div class="about-brand"><img src="/ruptura-logo.png" alt="" width="70" height="70"><div><div class="eyebrow">VERSION ${escapeHtml(about.version)}</div><h1>About ${escapeHtml(about.name)}</h1><p>${escapeHtml(about.description)}</p></div></div><div class="actions"><button class="button primary" data-action="workspace">Open your workspace →</button><button class="button" data-action="back-to-simulation">Back to experiment</button></div></div><div class="two-column"><div>${card("Citing RUPTURA",`<p class="card-description">Please cite the original RUPTURA publication when using this software in scientific work.</p>${reference}<div class="actions"><button class="button" data-action="download-citation">↧ Download BibTeX</button><a class="button" href="${escapeHtml(about.repository)}" target="_blank" rel="noopener noreferrer">Project & source ↗</a></div>`,"Original publication")}${card("Contact & license",`<p><strong>${escapeHtml(about.contact.name)}</strong><br><a href="mailto:${escapeHtml(about.contact.email)}">${escapeHtml(about.contact.email)}</a></p><p class="section-note">RUPTURA is open-source software distributed under the ${escapeHtml(about.license)} license.</p><button class="button" data-action="license">Read the full MIT license</button><p class="section-note">Simulations run locally, and inputs and results remain in your workspace.</p>`)}</div><div>${card("RUPTURA 2.0 authors",`<ol class="author-list">${about.authors.map(author=>`<li><strong>${escapeHtml(author.name)}</strong><span>${escapeHtml(author.affiliation)}</span></li>`).join("")}</ol>`,"Original RUPTURA authors, with Youri Ran and Muhammed Hassan")}</div></div>`;
}
function columnLength() {
  if (config().ColumnSections?.length) return config().ColumnSections.reduce((sum, section) => sum + (section.Length ?? section.InterfaceLength ?? 0), 0);
  if (config().Adsorbents?.length && config().Adsorbents.every(bed => bed.AdsorbentLength != null)) return config().Adsorbents.reduce((sum, bed) => sum + bed.AdsorbentLength, 0);
  return config().ColumnLength ?? 0.3;
}
function diagram() {
  const input = config(), beds = input.Adsorbents || [{Name: input.Geometry?.Type === "Monolith" ? "Monolith channels" : "Packed adsorbent"}];
  return `<div class="column-visual"><div class="flow-label">FEED <span>→ → → → → → → → →</span> OUTLET</div><div class="column-tube">${beds.map(bed => `<div class="bed-segment" style="flex:1">${escapeHtml(bed.Name)}</div>`).join("")}</div><div class="column-axis"><span>z = 0</span><span>${format(columnLength())} m</span></div></div>`;
}
function renderSetup() {
  const input = config(), mixture = input.SimulationType === "MixturePrediction";
  let main = card("Experiment", `<div class="form-grid">${field(describe("DisplayName", "Simulation name", null, null, "text"))}${field(describe("SimulationType", "Simulation type", "Breakthrough", state.schema.types))}</div>`);
  main += card("Operating conditions", fields("operating"), "The conditions shared by your components and column.");
  main += mixture ? card("Pressure range", fields("pressure")) : card("Time integration", fields("integration"), "Fixed steps for reproducible runs, or enter auto to run until breakthrough.");
  if (!mixture && input.BreakthroughIntegrator === "CVODE") main += card("CVODE controls", fields("cvode"));
  if (input.FluidPhase === "Liquid") main += card("Liquid & pH", fields("liquid"), "The acid/base component is a zero-based index from the component list.");
  if (input.SimulationType === "SwingAdsorption") {
    main += card("Swing phases", (input.SwingAdsorptionPhases || []).map((phase, i) => `<div class="site"><div class="site-header"><strong>Phase ${i + 1}</strong><button class="icon-button" data-action="remove-phase" data-index="${i}" aria-label="Remove phase ${i + 1}">✕</button></div><div class="form-grid">${[
      describe("Name", "Phase name", null, null, "text"), describe("Temperature", "Temperature [K] (blank inherits)"), describe("InletPressure", "Inlet pressure [Pa] (blank inherits)"), describe("NumberOfTimeSteps", "Time steps", 1000)
    ].map(def => field(def, ["SwingAdsorptionPhases", i])).join("")}</div></div>`).join("") + `<button class="button" data-action="add-phase">＋ Add phase</button>`);
  }
  const secondary = card("Column at a glance", diagram() + `<div class="mini-spec"><span>Geometry</span><strong>${escapeHtml(input.Geometry?.Type || "PackedBed")}</strong></div><div class="mini-spec"><span>Components</span><strong>${input.Components.length}</strong></div><div class="mini-spec"><span>Inlet pressure</span><strong>${format((input.InletPressure ?? 0) / 1000)} kPa</strong></div>`, "A view of your current configuration.");
  $("#content").innerHTML = `<div class="two-column"><div>${main}</div><div class="secondary">${secondary}${card("Your experiment", `<p class="section-note">Inputs and output files stay in your workspace. Each run keeps a snapshot of its configuration and an editable analysis notebook.</p>`)}${card("Quick start", `<p class="section-note">Choose an example or paste a complete simulation JSON to start a new experiment.</p><button class="button" data-action="new">Browse examples →</button>`)}</div></div>`;
}
function siteEditor(site, path, chemical) {
  const defs = state.schema.isotherms[site.Type] || state.schema.isotherms.Langmuir;
  const isothermSelect = field(describe("Type", "Equilibrium isotherm", "Langmuir", Object.keys(state.schema.isotherms)), path);
  let body = `<div class="form-grid">${isothermSelect}${defs.map(([key,label,value], i) => field(describe(i, label, value), [...path, "Parameters"])).join("")}</div>`;
  if (chemical) {
    const kinetics = site.RateEquation || "FirstOrder";
    body += `<div class="subheading">Kinetics</div><div class="form-grid">${field(describe("RateEquation", "Rate equation", "FirstOrder", Object.keys(state.schema.kinetics)), path)}${field(describe("MaximumLoading", "Maximum loading [mol/kg]", 0.8), path)}${field(describe("HeatOfAdsorption", "Heat of chemisorption [J/mol]"), path)}${state.schema.kinetics[kinetics].map(([key,label,value]) => field(describe(key,label,value,null,typeof value === "boolean" ? "boolean" : "number"), [...path,"RateEquationParameters"])).join("")}</div>`;
  } else {
    body += `<details style="margin-top:16px"><summary class="helper">Site transfer & heat overrides</summary><div class="form-grid" style="margin-top:12px">${field(describe("MassTransferCoefficient", "Site mass transfer [1/s]"), [...path,"RateEquationParameters"])}${field(describe("HeatOfAdsorption", "Site heat [J/mol]"), path)}</div></details>`;
  }
  return `<div class="site"><div class="site-header"><strong>${chemical ? "Chemisorption" : "Physisorption"} site ${path.at(-1) + 1}</strong><button class="icon-button" data-action="remove-site" ${pathAttr(path)} aria-label="Remove adsorption site">✕</button></div>${body}</div>`;
}
function renderComponents() {
  const liquid=config().FluidPhase==="Liquid",total=config().Components.reduce((sum,item)=>sum+(item.GasPhaseMolFraction||0),0);
  $("#content").innerHTML=`<div class="section-title"><div><h2>Components & adsorption</h2><p>Set the feed, species properties, and adsorption models for this simulation.</p></div><div class="actions">${liquid?"":`<button class="button" data-action="normalize-feed">Normalize feed · Σ ${format(total)}</button>`}<button class="button" data-action="add-component">＋ Add component</button></div></div>`+config().Components.map((component,i)=>{
    const path=["Components",i];
    let body=`<div class="form-grid">${field(describe("Name","Component name",null,null,"text"),path)}${liquid?"":field(describe("CarrierGas","Carrier gas",false,null,"boolean"),path)}</div><div style="height:18px"></div>`;
    body+=fields("component",path,def=>liquid?def.key!=="GasPhaseMolFraction":!["LiquidPhaseConcentration","InitialLiquidPhaseConcentration"].includes(def.key));
    {
      for(const [kind,label] of [["PhysisorptionSites","Physisorption"],["ChemisorptionSites","Chemisorption"]]) body+=`<div class="subheading">${label}</div>${(component[kind]||[]).map((site,j)=>siteEditor(site,[...path,kind,j],kind==="ChemisorptionSites")).join("")}<button class="button" data-action="add-site" data-index="${i}" data-kind="${kind}">＋ ${label} site</button>`;
    }
    return card(`${i} · ${component.Name}`,body,component.CarrierGas?"Carrier gas":`${(component.PhysisorptionSites||[]).length} physical · ${(component.ChemisorptionSites||[]).length} chemical sites`,`<button class="icon-button" data-action="remove-component" data-index="${i}" aria-label="Remove component ${escapeHtml(component.Name)}">✕</button>`);
  }).join("");
}
function renderColumn() {
  const input=config(),mode=!input.Adsorbents?.length?"single":input.Adsorbents[0].MixFraction!==undefined?"mixed":"layered";
  let html=`<div class="two-column"><div>${card("Column & transport",fields("column"),"Boundary conditions and physical properties for this simulation.")}${card("Geometry",fields("geometry",["Geometry"],def=>input.Geometry?.Type==="Monolith"?!["ColumnVoidFraction","ParticleDiameter","InternalDiameter"].includes(def.key):["Type","ColumnVoidFraction","ParticleDiameter","InternalDiameter","OuterDiameter"].includes(def.key)))}</div><div class="secondary">${card("Column structure",diagram()+`<div class="field" style="margin-top:20px">${fieldHeader("bed-mode","Adsorbent arrangement","Choose one adsorbent, separate beds, or a uniform mixture of adsorbents for this simulation.")}<select id="bed-mode"><option value="single" ${mode==="single"?"selected":""}>Single adsorbent</option><option value="layered" ${mode==="layered"?"selected":""}>Layered beds</option><option value="mixed" ${mode==="mixed"?"selected":""}>Uniform mixture</option></select></div>`)}${card("Spatial grid & initial state",jsonEditor(input.ColumnDistances||[],"grid")+`<div style="height:18px"></div>${field(describe("ReadColumnFile","Initial column state file",null,null,"text"))}<p class="helper">Explicit grid positions are supported for multibed columns. File paths are local to the server workspace.</p>`)}</div></div>`;
  if(mode!=="single") html+=card("Adsorbent beds",input.Adsorbents.map((bed,i)=>{
    const sectionIndex=input.ColumnSections?.findIndex(section=>section.Adsorbent===bed.Name)??-1;
    const sectionPath=sectionIndex>=0?["ColumnSections",sectionIndex]:["Adsorbents",i];
    const physical=[describe("Name","Adsorbent name",null,null,"text"),describe("ColumnVoidFraction","Void fraction",input.ColumnVoidFraction??.4),describe("ParticleDensity","Particle density [kg/m³]",input.ParticleDensity??1000),describe("ParticleDiameter","Particle diameter [m]",input.ParticleDiameter??.001)];
    let body=`<div class="form-grid">${physical.map(def=>field(def,["Adsorbents",i])).join("")}${mode==="mixed"?field(describe("MixFraction","Mix fraction",1/input.Adsorbents.length),["Adsorbents",i]):field(describe(sectionIndex>=0?"Length":"AdsorbentLength","Bed length [m]",columnLength()/input.Adsorbents.length),sectionPath)}${sectionIndex>=0?field(describe("NumberOfGridPoints","Grid intervals (blank = global)"),sectionPath):""}</div>`;
    if(sectionIndex>=0&&input.ColumnSections[sectionIndex+1]?.InterfaceLength!=null) body+=`<div class="subheading">Interface after this bed</div><div class="form-grid">${field(describe("InterfaceLength","Interface length [m]",0),["ColumnSections",sectionIndex+1])}${field(describe("NumberOfGridPoints","Interface grid intervals"),["ColumnSections",sectionIndex+1])}</div>`;
    body+=`<div class="subheading">Component-specific adsorption & kinetics</div><p class="helper">Set this bed's overrides by component name, or supply a Components array.</p>${jsonEditor(Object.fromEntries(["ComponentParameters","Components"].filter(key=>bed[key]!==undefined).map(key=>[key,bed[key]])),`bed-${i}`)}`;
    return `<div class="site"><div class="site-header"><strong>Bed ${i+1} · ${escapeHtml(bed.Name)}</strong><button class="icon-button" data-action="remove-bed" data-index="${i}" aria-label="Remove bed ${i+1}">✕</button></div>${body}</div>`;
  }).join("")+`<button class="button" data-action="add-bed">＋ Add bed</button>`);
  html+=card("Energy balance",fields("thermal",[],def=>input.FluidPhase==="Liquid"?!["gasThermalConductivity","heatTransferGasSolid","heatTransferGasWall","heatCapacityGas"].includes(def.key):!def.key.toLowerCase().includes("liquid")),"Enable coupled fluid, solid, and wall temperature evolution.");
  $("#content").innerHTML=html;
}
function renderReactions() {
  let html=`<div class="section-title"><div><h2>Reaction network</h2><p>Set stoichiometry, reaction orders, and kinetics for this simulation.</p></div><button class="button" data-action="add-reaction">＋ Add reaction</button></div>`;
  if(!(config().Reactions||[]).length) html+=card("Add chemistry to your column",'<p class="section-note">Configure power-law, Langmuir–Hinshelwood, or LHHW reactions coupled to the component balances.</p>');
  html+=(config().Reactions||[]).map((reaction,i)=>{
    const path=["Reactions",i],options=config().Components.map((component,j)=>({value:j,label:`${j} · ${component.Name}`}));
    let body=`<div class="form-grid">${field(describe("Phase","Phase","Physisorbed",["Physisorbed","Chemisorbed","PoreConcentration"]),path)}${field(describe("Style","Rate law","GeneralPowerLaw",["GeneralPowerLaw","LangmuirHinshelwood","LangmuirHinshelwoodHougenWatson"]),path)}${field(describe("Site","Site index",0),path)}${field(describe("RateLimitTime","Rate limit time [s]",1e-4),path)}</div>`;
    for(const side of ["Reactants","Products"]) {
      body+=`<div class="subheading">${side}</div>`+(reaction[side]||[]).map((_,j)=>{
        const offset=side==="Reactants"?j:(reaction.Reactants?.length||0)+j;
        const stoichPath=Array.isArray(reaction.Stoichiometry)?[...path,"Stoichiometry"]:[...path,"Stoichiometry",side];
        return `<div class="participant-row">${field(describe(j,"Component",null,options),[...path,side])}${field(describe(Array.isArray(reaction.Stoichiometry)?offset:j,"Stoichiometry",1),stoichPath)}${field(describe(j,"Order (optional)",1),[...path,side==="Reactants"?"ForwardOrders":"BackwardOrders"])}<button class="icon-button" data-action="remove-participant" data-index="${i}" data-side="${side}" data-row="${j}" aria-label="Remove participant">✕</button></div>`;
      }).join("")+`<button class="button" data-action="add-participant" data-index="${i}" data-side="${side}">＋ ${side==="Reactants"?"Reactant":"Product"}</button>`;
    }
    body+=`<div class="subheading">Kinetic constants</div><div class="form-grid">${[describe("forwardRateCoefficient","Forward rate coefficient",1),describe("forwardActivationEnergy","Activation energy [J/mol]",0),describe("equilibriumConstant","Equilibrium constant",1),describe("gibbsFreeEnergy","Gibbs free energy [J/mol]",0)].map(def=>field(def,[...path,"Kinetics"])).join("")}</div>`;
    return card(`Reaction ${i+1}`,body,"",`<button class="icon-button" data-action="remove-reaction" data-index="${i}" aria-label="Remove reaction ${i+1}">✕</button>`);
  }).join("");
  $("#content").innerHTML=html;
}
function renderJson() { $("#content").innerHTML = card("Simulation input", `<p class="card-description">This is the complete input for the selected simulation. Form controls preserve fields they do not edit. Validation uses the C++ input reader.</p>${jsonEditor(config(),"simulation",false,"Validate & apply JSON")}<div class="actions" style="margin-top:20px"><button class="button" data-action="export-input">↧ Download simulation.json</button></div>`); }

function runsForCurrent() { return state.runs.filter(record => record.simulationId === simulation().id || record.id === simulation().lastRunId); }
function ensureRun() { if (!activeRun()) state.selectedRun = simulation().lastRunId || runsForCurrent()[0]?.id || state.runs[0]?.id || null; }
function runSelector() {
  return `<div class="toolbar"><div class="field"><label for="run-selector">Saved run</label><select id="run-selector">${state.runs.map(record => `<option value="${escapeHtml(record.id)}" ${record.id === state.selectedRun ? "selected" : ""}>${escapeHtml(record.name || record.id)} · ${escapeHtml(record.status)} · ${escapeHtml(record.id.slice(0,15))}</option>`).join("")}</select></div>${ongoing(activeRun()) ? `<button class="button danger" data-action="cancel">■ Cancel run</button>` : ""}<button class="button" data-action="export-input-run">↧ Run input</button></div>`;
}
function emptyResults() { return `<section class="card"><div class="empty"><div class="empty-icon">⌁</div><h2>Your results start here</h2><p>Run a simulation to explore breakthrough curves, isotherms, temperature histories, and profiles along the column.</p><button class="button primary" data-action="run">▶ Run simulation</button></div></section>`; }
async function renderResults() {
  const request = ++state.resultRequest;
  ensureRun();
  if (!state.selectedRun) { $("#content").innerHTML = emptyResults(); return; }
  const id = state.selectedRun, record = activeRun();
  $("#content").innerHTML = runSelector() + `<section class="card"><div class="empty"><div class="spinner"></div><h2>${ongoing(record) ? "Your simulation is running" : "Loading results"}</h2><p>${ongoing(record) ? "Plots will appear when the run finishes. Follow its progress in Run log." : "Reading the saved output files."}</p></div></section>`;
  if (ongoing(record)) return;
  try {
    const [meta, info] = await Promise.all([api(`/api/runs/${id}/plots`),api(`/api/runs/${id}`)]);
    if (request !== state.resultRequest || state.activeTab !== "results" || state.selectedRun !== id) return;
    state.metadata = meta; state.runInfo = info;
    const options = state.plotOptions;
    if (!meta.plots.some(item => item.id === options.plot)) options.plot = meta.plots[0].id;
    if (options.components == null || state.plotRun !== id) options.components = meta.components.filter(item => !item.carrier || item.adsorbing).map(item => item.index);
    state.plotRun = id;
    options.frame = Math.min(options.frame ?? Math.max(meta.frames-1,0),Math.max(meta.frames-1,0));
    options.grid = Math.min(options.grid ?? Math.max((meta.positions?.length || 1)-1,0),Math.max((meta.positions?.length || 1)-1,0));
    options.metric = meta.metrics?.some(item => item.id === options.metric) ? options.metric : "C";
    renderResultLayout();
    await updatePlot();
  } catch (error) {
    if (request !== state.resultRequest || state.activeTab !== "results") return;
    $("#content").innerHTML = runSelector() + `<section class="card"><div class="empty"><div class="empty-icon">⌁</div><h2>Results unavailable</h2><p>${escapeHtml(error.message)}</p><button class="button" data-tab="log">View run log</button></div></section>`;
  }
}
function plotSelect(id,label,values,value) {
  return `<div class="field"><label for="plot-${id}">${label}</label><select id="plot-${id}" data-plot="${id}">${values.map(item => `<option value="${escapeHtml(item.value ?? item.id)}" ${String(item.value ?? item.id) === String(value) ? "selected" : ""}>${escapeHtml(item.label.replace(/<[^>]*>/g,""))}</option>`).join("")}</select></div>`;
}
function profileSliderLabel() {
  const options=state.plotOptions,meta=state.metadata;
  if(options.xAxis==="time") return `<span>POSITION ${options.grid+1} / ${meta.positions.length}</span><span>${format(meta.positions[options.grid])} m</span>`;
  const units=options.xUnits||"min",factor=units==="s"?60:units==="hr"?1/60:1;
  return `<span>FRAME ${options.frame+1} / ${meta.frames}</span><span>${format(meta.times[options.frame]*factor)} ${units}</span>`;
}
function renderResultLayout() {
  const meta=state.metadata,options=state.plotOptions,profile=options.plot==="profile";
  const mixture=meta.kind==="MixturePrediction",paired=meta.energyBalance&&options.plot==="breakthrough";
  options.xAxis??="position";options.fps??=20;
  let controls=plotSelect("plot","Figure",meta.plots,options.plot);
  if(profile) controls+=plotSelect("xAxis","X-axis",[{value:"position",label:"Position"},{value:"time",label:"Time"}],options.xAxis);
  if(!mixture&&(!profile||options.xAxis==="time")) controls+=plotSelect("xUnits","Time units",[{value:"min",label:"Minutes"},{value:"s",label:"Seconds"},{value:"hr",label:"Hours"}],options.xUnits||"min");
  if(options.plot==="breakthrough") controls+=plotSelect("yUnits","Vertical axis",[{value:"normalized_concentration",label:"Normalized concentration"},{value:"concentration",label:"Concentration [mol/m³]"},...(!meta.liquid?[{value:"molefraction",label:"Mole fraction"},{value:"partial_pressure",label:"Partial pressure [Pa]"}]:[])],options.yUnits||"normalized_concentration");
  if(options.plot==="temperature"||paired) controls+=plotSelect("temperatureUnits","Temperature units",[{value:"kelvin",label:"Kelvin"},{value:"celsius",label:"Celsius"}],options.temperatureUnits||"kelvin");
  if(profile) controls+=plotSelect("metric","Column metric",meta.metrics,options.metric);
  if(!mixture&&!profile) controls+=plotSelect("grid","Column position",meta.positions.map((z,i)=>({value:i,label:`${format(z)} m${i===meta.positions.length-1?" · outlet":i===0?" · inlet":""}`})),options.grid);
  const columnMetric=profile&&meta.metrics.find(item=>item.id===options.metric)?.source==="column";
  if(!columnMetric&&options.plot!=="temperature") controls+=`<div class="field"><label>Components</label>${meta.components.map(item=>`<label class="component-choice"><input type="checkbox" data-plot-component="${item.index}" ${options.components.includes(item.index)?"checked":""}>${escapeHtml(item.name)}${item.carrier?" · carrier":""}</label>`).join("")}</div>`;
  controls+=`<label class="component-choice"><input type="checkbox" data-plot="markers" ${options.markers?"checked":""}>Show markers</label>`;
  if(profile) controls+=`<div class="video-export"><div class="field">${fieldHeader("video-fps","Frames per second","Set the playback and MP4 frame rate from 1 to 60 frames per second.")}<input id="video-fps" type="number" min="1" max="60" step="1" value="${options.fps}"></div><button class="button" data-action="export-video">↧ Download MP4</button><p class="helper">Exports ${options.xAxis==="time"?meta.positions.length:meta.frames} frames, sweeping ${options.xAxis==="time"?"position":"time"} with the current metric and components.</p></div>`;
  const position=options.xAxis==="time",sliderId=position?"position-slider":"frame-slider";
  const slider=profile?`<div class="range-row"><label for="${sliderId}">${profileSliderLabel()}</label><div class="playback"><button class="button" data-action="play" id="play-button" aria-label="Play column profiles">▶</button><input id="${sliderId}" aria-label="${position?"Column position":"Column time frame"}" type="range" min="0" max="${position?meta.positions.length-1:meta.frames-1}" value="${position?options.grid:options.frame}"></div></div>`:"";
  const second=paired?`<section class="card plot-card"><div class="card-header"><div><h2>Temperature profile</h2><p>Fluid, solid, and wall temperature at the selected column position.</p></div></div><div class="card-body"><div id="temperature-canvas" class="plot-canvas"></div></div></section>`:"";
  $("#content").innerHTML=runSelector()+`<div class="results-grid"><section class="card plot-options"><div class="card-header"><h2>Explore results</h2></div><div class="card-body">${controls}</div></section><div><section class="card plot-card"><div class="card-header"><div><h2>${escapeHtml(meta.plots.find(item=>item.id===options.plot).label)}</h2><p>Saved output · ${escapeHtml(activeRun().name||activeRun().id)}</p></div>${pill(activeRun().status)}</div>${slider}<div class="card-body"><div id="plot-canvas" class="plot-canvas"></div></div><div id="plot-error" class="plot-info"></div>${meta.sampled?`<div class="plot-info">Displaying ${meta.frames} sampled frames from ${format(meta.totalFrames)} saved frames.</div>`:""}</section>${second}${card("Keep exploring",`<div class="actions"><button class="button" data-action="report">↧ Interactive HTML report</button><button class="button" data-action="download" data-file="analysis.ipynb">↧ Analysis notebook</button><button class="button subtle" data-tab="log">View run log →</button></div><p class="helper">The notebook uses the same plotting library. Its editable cells remain available for custom analysis.</p>`)}${card("Output files",`<div class="file-list">${state.runInfo.files.map(file=>`<button data-action="download" data-file="${escapeHtml(file.name)}">${escapeHtml(file.name)} <span>${format(file.bytes/1024)} KB</span></button>`).join("")}</div>`)}</div></div>`;
}
async function updatePlot() {
  const canvas=$("#plot-canvas");if(!canvas) return;
  if(!window.Plotly) {$("#plot-error").textContent="Plotly is missing. Install ruptura[web] in the Python environment shown in Settings, then restart the Lab.";return;}
  const request=++state.chartRequest,id=state.selectedRun,options={...clone(state.plotOptions),theme:state.theme};
  const targets=[{canvas,options}];
  const temperature=$("#temperature-canvas");if(temperature) targets.push({canvas:temperature,options:{...options,plot:"temperature"}});
  $("#plot-error").textContent="Loading figure…";
  try {
    const figures=await Promise.all(targets.map(target=>api(`/api/runs/${id}/plot`,target.options)));
    if(request!==state.chartRequest||state.activeTab!=="results"||id!==state.selectedRun||!canvas.isConnected) return;
    await Promise.all(targets.map((target,i)=>Plotly.react(target.canvas,figures[i].data,figures[i].layout,{responsive:true,displaylogo:false,toImageButtonOptions:{format:"svg",filename:i?"ruptura-temperature":"ruptura-plot"}})));
    if(request===state.chartRequest&&$("#plot-error")) $("#plot-error").textContent="Hover to inspect · drag to zoom · click a legend entry to toggle a component.";
  } catch(error) {if(request===state.chartRequest&&$("#plot-error")) $("#plot-error").textContent=error.message;}
}
async function exportVideo(button) {
  const fps=Number($("#video-fps").value);
  if(!Number.isInteger(fps)||fps<1||fps>60) throw new Error("Frames per second must be a whole number from 1 to 60.");
  const id=state.selectedRun,options={...clone(state.plotOptions),fps,theme:state.theme};
  button.disabled=true;button.textContent="Rendering MP4…";stopPlayback();
  toast("Rendering the column animation locally. The download will start when it is ready.");
  try {const result=await api(`/api/runs/${id}/video`,options);await downloadFile(result.name,id);toast("MP4 ready.");}
  finally {button.disabled=false;button.textContent="↧ Download MP4";}
}
async function renderLog() {
  ensureRun(); logOffset = null;
  if (!state.selectedRun) { $("#content").innerHTML = emptyResults(); return; }
  $("#content").innerHTML = runSelector() + card("Simulation log", `<p class="card-description">${escapeHtml(activeRun()?.directory || "")}</p><pre id="run-log" class="log" aria-label="Simulation output">Loading…</pre>`);
  await updateLog();
}
async function updateLog() {
  const node = $("#run-log"), id = state.selectedRun; if (!node || !id) return;
  const data = await api(`/api/runs/${id}/log${logOffset == null ? "" : `?offset=${logOffset}`}`);
  if (!node.isConnected || id !== state.selectedRun) return;
  const nearBottom = node.scrollHeight-node.scrollTop-node.clientHeight < 50;
  node.textContent = logOffset == null ? data.text || "Waiting for output…" : (node.textContent+data.text).slice(-100000);
  logOffset = data.offset;
  if (nearBottom || logOffset == null) node.scrollTop = node.scrollHeight;
}

function modal(title, body) {
  stopPlayback();
  $("#modal-content").innerHTML = `<div class="modal-header"><h2>${escapeHtml(title)}</h2><button class="icon-button" data-action="close-modal" aria-label="Close dialog">✕</button></div><div class="modal-body">${body}</div>`;
  if (!$("#modal").open) $("#modal").showModal();
}
function newDialog() {
  const group = paper => `<div class="preset-grid">${Object.keys(state.presets).filter(key => (state.presetInfo[key]?.group === "paper") === paper).map(key => {
    const info = state.presetInfo[key] || {};
    return `<button class="preset-button" data-action="create" data-preset="${escapeHtml(key)}"><strong>${escapeHtml(info.title || state.presets[key].DisplayName)}</strong><span>${escapeHtml(info.description || "")}</span></button>`;
  }).join("")}</div>`;
  const pasteHelp = "Paste one complete simulation JSON object to validate it and create a new untitled simulation.";
  const clear = `<button type="button" class="reset" data-action="clear-pasted-json" title="Clear pasted JSON" aria-label="Clear pasted JSON">↺</button>`;
  modal("Quick start",`<p class="card-description">Choose an example or paste your own input below.</p><h3 class="preset-heading">From the paper</h3>${group(true)}<h3 class="preset-heading">Quick starts</h3>${group(false)}<section class="paste-json"><h3 class="preset-heading">From JSON</h3>${fieldHeader("new-simulation-json", "Simulation JSON", pasteHelp, clear)}<textarea id="new-simulation-json" class="json-editor compact" spellcheck="false" placeholder='{"SimulationType": "Breakthrough", "Components": […]}'></textarea><p id="new-json-error" class="inline-error" role="alert" hidden></p><div class="editor-actions"><span class="helper">New simulations receive an untitled name.</span><button class="button primary" data-action="create-json">Create from JSON</button></div></section>`);
}
async function browse(path) {
  const data = await api(`/api/files?path=${encodeURIComponent(path || state.root)}`);
  modal("Import from your workspace",`<p class="card-description">Open a simulation JSON or saved Qt/web Lab workspace. Local file selection preserves relative restart-file paths.</p>${fieldHeader("browse-path", "Directory", "Enter a local directory to browse for simulation inputs and saved workspaces.")}<div class="browser-path"><button class="button" data-action="browse" data-directory="${escapeHtml(data.parent)}" aria-label="Parent folder">↑</button><input class="text-input" id="browse-path" aria-label="Directory path" value="${escapeHtml(data.path)}"><button class="button" data-action="browse-path">Go</button></div><div class="file-browser">${data.entries.map(item => `<button data-action="${item.directory ? "browse" : "import-path"}" data-directory="${escapeHtml(item.path)}"><span>${item.directory ? "▱" : "{ }"}</span>${escapeHtml(item.name)}</button>`).join("") || `<p class="helper" style="padding:14px">No folders or JSON files in this directory.</p>`}</div><div class="editor-actions"><span class="helper">Or upload a JSON file from your browser.</span><button class="button" data-action="upload">Choose file…</button></div>${workspacePasteBox()}`);
}
function workspacePasteBox() {
  const clear = `<button type="button" class="reset" data-action="clear-import-json" title="Clear pasted JSON" aria-label="Clear workspace JSON">↺</button>`;
  return `<section class="paste-json"><h3 class="preset-heading">Paste JSON</h3>${fieldHeader("import-workspace-json", "Workspace or simulation JSON", "Paste a saved Qt or web Lab workspace, or a single simulation input, to add its simulations to this workspace.", clear)}<textarea id="import-workspace-json" class="json-editor compact" spellcheck="false" placeholder='{"Format": "RupturaWebLabState", …}'></textarea><p id="import-json-error" class="inline-error" role="alert" hidden></p><div class="editor-actions"><span class="helper">Adds simulations to the current workspace.</span><button class="button primary" data-action="import-json">Import JSON</button></div></section>`;
}
function nextUntitledName() {
  let number = Math.max(1, Number(state.workspace.nextUntitledIndex) || 1);
  for (const item of state.workspace.simulations) {
    const match = /^untitled_(\d+)$/i.exec(item.input.DisplayName || "");
    if (match) number = Math.max(number, Number(match[1]) + 1);
  }
  state.workspace.nextUntitledIndex = number + 1;
  return `untitled_${number}`;
}
function createSimulation(input,notes="") {
  const item={id:uid(),notes,input:{...clone(input),DisplayName:nextUntitledName()}};
  state.workspace.simulations.push(item);state.workspace.selectedId=item.id;
  state.selectedRun=null;state.metadata=null;state.plotOptions={};
  notice();$("#modal").close();selectTab("setup");render();
  return item;
}
async function importWorkspace(payload) {
  const imported = await api("/api/import",payload);
  if (imported.documentType === "simulation") for (const item of imported.simulations) item.input.DisplayName = nextUntitledName();
  const existing = state.workspace.simulations;
  for (const item of imported.simulations) { if (existing.some(old => old.id === item.id)) { if (imported.selectedId === item.id) imported.selectedId = item.id = uid(); else item.id = uid(); } }
  state.workspace.simulations.push(...imported.simulations);
  state.workspace.selectedId = imported.selectedId;
  if (imported.legacyLibrary) state.workspace.legacyLibrary = imported.legacyLibrary;
  state.selectedRun = null; state.metadata = null;
  $("#modal").close(); selectTab("setup"); render(); await saveNow(); toast(`Imported ${imported.simulations.length} simulation${imported.simulations.length === 1 ? "" : "s"}.`);
}
function downloadJson(value,name) {
  const blob = new Blob([JSON.stringify(value,null,2)+"\n"],{type:"application/json"}); downloadBlob(blob,name);
}
function downloadBlob(blob,name) {
  const link = document.createElement("a"), url = URL.createObjectURL(blob); link.href = url; link.download = name; link.click(); setTimeout(() => URL.revokeObjectURL(url),1000);
}
async function downloadFile(name,id=state.selectedRun) {
  const response = await fetch(`/api/runs/${id}/file/${encodeURIComponent(name)}`,{headers:{"X-Ruptura-Session":state.session}});
  if (!response.ok) { const error = await response.json(); throw new Error(error.error); }
  downloadBlob(await response.blob(),name);
}
async function run() {
  const button = $("#run-button"); button.disabled = true; notice();
  try {
    await saveNow();
    const record = await api("/api/runs",{input:config(),simulationId:simulation().id});
    state.runs.unshift(record); simulation().lastRunId = record.id; state.selectedRun = record.id; state.metadata = null; state.plotOptions = {};
    await saveNow(); header(); selectTab("log"); toast("Simulation submitted to the local engine.");
  } finally { button.disabled = false; }
}

function syncLengths() {
  if (config().ColumnSections?.length) config().ColumnLength = config().ColumnSections.reduce((sum,section) => sum+(section.Length ?? section.InterfaceLength ?? 0),0);
  else if (config().Adsorbents?.length && config().Adsorbents.every(bed => bed.AdsorbentLength != null)) config().ColumnLength = config().Adsorbents.reduce((sum,bed) => sum+bed.AdsorbentLength,0);
}
function arrangeBeds(mode) {
  const input = config(), length = columnLength();
  if (mode === "single") { delete input.Adsorbents; delete input.ColumnSections; return; }
  input.Adsorbents = input.Adsorbents?.length ? input.Adsorbents : [{Name:"Front bed"},{Name:"Polishing bed"}];
  delete input.ColumnSections;
  if (mode === "mixed") input.Adsorbents.forEach(bed => { delete bed.AdsorbentLength; bed.MixFraction = 1/input.Adsorbents.length; });
  else {
    input.ColumnSections = [];
    input.Adsorbents.forEach((bed,i) => { delete bed.MixFraction; delete bed.AdsorbentLength; input.ColumnSections.push({Adsorbent:bed.Name,Length:length/input.Adsorbents.length}); if (i+1<input.Adsorbents.length) input.ColumnSections.push({InterfaceLength:0}); });
  }
  input.ColumnLength = length;
}
function ensureSections() {
  const input=config();if(input.ColumnSections?.length) return;
  const length=columnLength();input.ColumnSections=[];
  input.Adsorbents.forEach((bed,i)=>{input.ColumnSections.push({Adsorbent:bed.Name,Length:bed.AdsorbentLength??length/input.Adsorbents.length});delete bed.AdsorbentLength;if(i+1<input.Adsorbents.length) input.ColumnSections.push({InterfaceLength:0});});
}
function renameComponent(oldName,newName) {
  for (const bed of config().Adsorbents || []) {
    if (bed.ComponentParameters && Object.hasOwn(bed.ComponentParameters,oldName)) { const value = bed.ComponentParameters[oldName]; delete bed.ComponentParameters[oldName]; bed.ComponentParameters[newName] = value; }
    for (const component of bed.Components || []) if (component.Name === oldName) component.Name = newName;
  }
  for (const reaction of config().Reactions || []) for (const side of ["Reactants","Products"]) reaction[side] = reaction[side].map(item => typeof item === "string" && item.toLowerCase() === oldName.toLowerCase() ? newName : item);
}
function initializeReactionArray(path) {
  if (path[0] !== "Reactions") return;
  if (path[2] === "Stoichiometry") {
    const reaction=config().Reactions[path[1]];
    if(!reaction.Stoichiometry) reaction.Stoichiometry={};
    if(!Array.isArray(reaction.Stoichiometry)) for(const side of ["Reactants","Products"]) if(!reaction.Stoichiometry[side]) reaction.Stoichiometry[side]=reaction[side].map(()=>1);
    return;
  }
  if (!["ForwardOrders","BackwardOrders"].includes(path[2]) || config().Reactions[path[1]][path[2]]) return;
  const reaction = config().Reactions[path[1]], side = path[2] === "ForwardOrders" ? "Reactants" : "Products";
  const offset = side === "Reactants" ? 0 : reaction.Reactants.length;
  reaction[path[2]] = reaction[side].map((_,i) => Array.isArray(reaction.Stoichiometry) ? reaction.Stoichiometry[offset+i] ?? 1 : reaction.Stoichiometry?.[side]?.[i] ?? 1);
}
function changeField(target) {
  const path = JSON.parse(target.dataset.path), previous = getPath(config(),path);
  let value = target.type === "checkbox" ? target.checked : target.value;
  if(target.tagName === "SELECT" && value === "" && target.dataset.default !== undefined) value=JSON.parse(target.dataset.default);
  if (target.type !== "checkbox" && value === "") {
    if (typeof path.at(-1) === "number") { target.setCustomValidity("Enter a value for this parameter."); target.reportValidity(); return; }
    value = undefined;
  } else if (target.type === "number" || (path.at(-1) === "NumberOfTimeSteps" && value !== "auto") || (path[0] === "Reactions" && ["Reactants","Products"].includes(path[2]))) {
    value = Number(value); if (!Number.isFinite(value)) { target.setCustomValidity("Enter a finite number or auto for time steps."); target.reportValidity(); return; }
  }
  target.setCustomValidity(""); initializeReactionArray(path); setPath(config(),path,value);
  let refresh = target.tagName === "SELECT" || target.type === "checkbox";
  if (path[0] === "Components" && path.at(-1) === "Name" && previous !== value) renameComponent(previous,value);
  if (path[0] === "Adsorbents" && path.at(-1) === "Name") for (const section of config().ColumnSections || []) if (section.Adsorbent === previous) section.Adsorbent = value;
  if (path.at(-1) === "Type" && ["PhysisorptionSites","ChemisorptionSites"].some(key => path.includes(key))) {
    const site = getPath(config(),path.slice(0,-1)); site.Parameters = state.schema.isotherms[value].map(def => def[2]);
  }
  if (path.at(-1) === "RateEquation") { const site = getPath(config(),path.slice(0,-1)); site.RateEquationParameters = Object.fromEntries(state.schema.kinetics[value].map(([key,,defaultValue]) => [key,defaultValue])); }
  if (path[0] === "Geometry" && path.at(-1) === "Type" && value === "Monolith") Object.assign(config().Geometry,{ChannelShape:config().Geometry.ChannelShape || "circular",InternalChannelDimension:config().Geometry.InternalChannelDimension || .001,NumberOfChannels:config().Geometry.NumberOfChannels || 10,OuterDiameter:config().Geometry.OuterDiameter || .01});
  if (path.length === 1 && ["ColumnVoidFraction","ParticleDiameter"].includes(path[0]) && config().Geometry?.[path[0]] === previous) setPath(config(),["Geometry",path[0]],value);
  if (path.at(-1) === "CarrierGas" && value) config().Components.forEach((item,i) => { if (i !== path[1]) delete item.CarrierGas; });
  if (path[0] === "FluidPhase" && value === "Liquid") { config().Components.forEach(item => delete item.CarrierGas); config().MixturePredictionMethod = "SPI"; }
  if (path[0] === "SimulationType" && value === "SwingAdsorption" && !config().SwingAdsorptionPhases?.length) config().SwingAdsorptionPhases = [{Name:"Adsorption",NumberOfTimeSteps:500},{Name:"Regeneration",InletPressure:20000,NumberOfTimeSteps:500}];
  if (path[0] === "SimulationType" && value === "MixturePrediction") Object.assign(config(),{PressureStart:config().PressureStart ?? 1,PressureEnd:config().PressureEnd ?? 1e6,NumberOfPressurePoints:config().NumberOfPressurePoints ?? 100});
  if (path[0] === "SimulationType" && ["Breakthrough","SwingAdsorption"].includes(value)) {
    for(const key of ["Geometry","BoundaryCondition","InletPressure","ColumnEntranceVelocity","ColumnLength","ColumnVoidFraction","ParticleDensity","ParticleDiameter","DynamicViscosity","BreakthroughIntegrator","NumberOfTimeSteps","TimeStep","PrintEvery","WriteEvery","NumberOfGridPoints"])
      if(config()[key]===undefined) config()[key]=clone(state.presets.breakthrough[key]);
  }
  if (path[0] === "Reactions" && path.at(-1) === "Style" && value !== "GeneralPowerLaw") config().Reactions[path[1]].Phase = "PoreConcentration";
  if (path.includes("ColumnSections") || path.at(-1) === "AdsorbentLength") syncLengths();
  scheduleSave(); header(); if (refresh) renderContent();
}

async function action(name,button) {
  const index = Number(button?.dataset.index);
  switch(name) {
    case "new": newDialog(); break;
    case "close-modal": $("#modal").close(); break;
    case "create": {
      createSimulation(state.presets[button.dataset.preset],state.presetNotes[button.dataset.preset]||""); await saveNow(); break;
    }
    case "clear-pasted-json": $("#new-simulation-json").value=""; $("#new-json-error").hidden=true; break;
    case "create-json": {
      const errorNode = $("#new-json-error"); errorNode.hidden = true; button.disabled = true;
      try {
        const text = $("#new-simulation-json").value.trim();
        if (!text) throw new Error("Paste a simulation JSON object first.");
        let document;
        try { document = JSON.parse(text); } catch { throw new Error("The pasted text is not valid JSON."); }
        if (!document || Array.isArray(document) || typeof document !== "object" || document.Format || document.Simulations)
          throw new Error("Paste a single simulation input; use Import for a saved workspace.");
        const checked = await api("/api/validate", {input:document});
        createSimulation(checked.input); await saveNow(); toast("Simulation created from JSON.");
      } catch (error) { errorNode.textContent = error.message; errorNode.hidden = false; }
      finally { button.disabled = false; }
      break;
    }
    case "select-simulation": state.workspace.selectedId = button.dataset.id; state.selectedRun = simulation().lastRunId || null; state.metadata = null; state.plotOptions = {}; notice(); render(); window.scrollTo({top:0,behavior:"instant"}); scheduleSave(); break;
    case "duplicate": {
      const source = button.dataset.id ? state.workspace.simulations.find(item=>item.id===button.dataset.id) : simulation();
      if (!source) throw new Error("The source simulation is no longer in this workspace.");
      const item=clone(source);item.id=uid();item.input.DisplayName=nextUntitledName();delete item.lastRunId;delete item.presetId;
      state.workspace.simulations.push(item);state.workspace.selectedId=item.id;
      state.selectedRun=null;state.metadata=null;state.plotOptions={};notice();selectTab("setup");render();await saveNow();break;
    }
    case "delete-simulation": {
      if(state.workspace.simulations.length === 1) throw new Error("Keep at least one simulation in the workspace.");
      const id=button.dataset.id || simulation().id;
      if(!state.workspace.simulations.some(item=>item.id===id)) throw new Error("The simulation is no longer in this workspace.");
      state.workspace.simulations=state.workspace.simulations.filter(item=>item.id!==id);
      if(state.workspace.selectedId===id) {
        state.workspace.selectedId=state.workspace.simulations[0].id; state.selectedRun=null; state.metadata=null; state.plotOptions={}; if(state.page==="simulation") render(); else header();
      } else header();
      await saveNow(); toast("Simulation deleted."); break;
    }
    case "copy-simulation-json": {
      const source=state.workspace.simulations.find(item=>item.id===button.dataset.id);
      if(!source) throw new Error("The simulation is no longer in this workspace.");
      await copyToClipboard(JSON.stringify(source.input,null,2)+"\n"); toast("Simulation JSON copied to clipboard."); break;
    }
    case "import": await browse(state.root); break;
    case "browse": await browse(button.dataset.directory); break;
    case "browse-path": await browse($("#browse-path").value); break;
    case "import-path": await importWorkspace({path:button.dataset.directory}); break;
    case "clear-import-json": $("#import-workspace-json").value=""; $("#import-json-error").hidden=true; break;
    case "import-json": {
      const errorNode=$("#import-json-error"); errorNode.hidden=true; button.disabled=true;
      try {
        const text=$("#import-workspace-json").value.trim();
        if (!text) throw new Error("Paste workspace or simulation JSON first.");
        let document;
        try {document=JSON.parse(text);} catch {throw new Error("The pasted text is not valid JSON.");}
        if (!document || Array.isArray(document) || typeof document!=="object") throw new Error("Paste a workspace or simulation JSON object.");
        await importWorkspace({document});
      } catch(error) {errorNode.textContent=error.message; errorNode.hidden=false;}
      finally {button.disabled=false;}
      break;
    }
    case "upload": $("#upload").click(); break;
    case "export-state": await saveNow(); downloadJson(state.workspace,"ruptura-workspace.json"); break;
    case "export-input": downloadJson(config(),"simulation.json"); break;
    case "export-input-run": { const info=await api(`/api/runs/${state.selectedRun}`); downloadJson(info.input,"simulation.json"); break; }
    case "validate": { const node=$("#validate-button"); node.disabled=true; notice(); try { await api("/api/validate",{input:config()}); notice("Input validated by the Ruptura C++ engine. Ready to run."); } finally { node.disabled=false; } break; }
    case "run": await run(); break;
    case "workspace": selectTab(state.activeTab);break;
    case "license": modal("MIT license",`<pre class="license-text">${escapeHtml(state.about.licenseText)}</pre>`);break;
    case "about": renderAbout(); break;
    case "back-to-simulation": selectTab(state.activeTab); break;
    case "download-citation": downloadBlob(new Blob([state.about.citation.bibtex],{type:"application/x-bibtex"}),"ruptura.bib"); break;
    case "cancel": { const record=await api(`/api/runs/${state.selectedRun}/cancel`,{}); state.runs=state.runs.map(item=>item.id===record.id?record:item); header(); renderContent(); break; }
    case "settings": modal("Settings",`<div class="field">${fieldHeader("theme-select","Appearance","Choose a dark or light workspace; your preference is remembered in this browser.")}<select id="theme-select"><option value="dark" ${state.theme==="dark"?"selected":""}>Dark mode</option><option value="light" ${state.theme==="light"?"selected":""}>Light mode</option></select></div><h3 class="subheading">Local engine</h3><div class="field">${fieldHeader("cli-path", "Ruptura executable", "Enter the path to the local ruptura binary that validates and runs your simulations.")}<input id="cli-path" class="text-input" value="${escapeHtml(state.capabilities.cli || "")}" placeholder="/path/to/ruptura"></div><div class="subtle-box" style="margin:18px 0"><strong>Python</strong><br>${escapeHtml(state.python)}<br><br><strong>Workspace</strong><br>${escapeHtml(state.root)}<br><br>CVODE: ${state.capabilities.sundials ? "available" : "unavailable"}</div><button class="button primary" data-action="save-settings">Save engine path</button>`); break;
    case "save-settings": state.capabilities=await api("/api/settings",{cli:$("#cli-path").value}); $("#modal").close(); header(); notice(); toast("Engine connected."); break;
    case "add-component": {
      let number=config().Components.length+1;while(config().Components.some(item=>item.Name===`Component ${number}`)) number++;
      config().Components.push({Name:`Component ${number}`,...(config().FluidPhase==="Liquid"?{LiquidPhaseConcentration:0}:{GasPhaseMolFraction:0}),MassTransferCoefficient:.08,PhysisorptionSites:[{Type:"Langmuir",Parameters:[1,1e-5]}]});render();scheduleSave();break;
    }
    case "normalize-feed": {const sum=config().Components.reduce((total,item)=>total+(item.GasPhaseMolFraction||0),0);if(sum<=0)throw new Error("Set a positive feed fraction before normalizing.");config().Components.forEach(item=>item.GasPhaseMolFraction=(item.GasPhaseMolFraction||0)/sum);renderContent();scheduleSave();break;}
    case "remove-component": {
      const old=config().Components[index].Name;
      if((config().Reactions||[]).some(reaction=>[...(reaction.Reactants||[]),...(reaction.Products||[])].some(item=>item===index||typeof item==="string"&&item.toLowerCase()===old.toLowerCase()))) throw new Error("Remove this component from the reaction network first.");
      if(config().Components.length===1) throw new Error("Keep at least one component.");
      config().Components.splice(index,1);
      for(const reaction of config().Reactions||[]) for(const side of ["Reactants","Products"]) reaction[side]=(reaction[side]||[]).map(item=>typeof item==="number"&&item>index?item-1:item);
      for(const bed of config().Adsorbents||[]) {if(bed.ComponentParameters) for(const name of Object.keys(bed.ComponentParameters)) if(name.toLowerCase()===old.toLowerCase()) delete bed.ComponentParameters[name];if(bed.Components) bed.Components=bed.Components.filter(item=>item.Name.toLowerCase()!==old.toLowerCase());}
      if(config().pHComponent===index) delete config().pHComponent;else if(config().pHComponent>index) config().pHComponent--;
      render();scheduleSave();break;
    }
    case "add-site": { const key=button.dataset.kind; const site={Type:"Langmuir",Parameters:[1,1e-5]}; if(key==="ChemisorptionSites") Object.assign(site,{RateEquation:"FirstOrder",RateEquationParameters:{MassTransferCoefficient:.04},MaximumLoading:.8,HeatOfAdsorption:42000}); (config().Components[index][key]??=[]).push(site); renderContent(); scheduleSave(); break; }
    case "remove-site": { const path=JSON.parse(button.dataset.path); getPath(config(),path.slice(0,-1)).splice(path.at(-1),1); renderContent(); scheduleSave(); break; }
    case "add-phase": (config().SwingAdsorptionPhases??=[]).push({Name:`Phase ${(config().SwingAdsorptionPhases||[]).length+1}`,NumberOfTimeSteps:500}); renderContent(); scheduleSave(); break;
    case "remove-phase": config().SwingAdsorptionPhases.splice(index,1); renderContent(); scheduleSave(); break;
    case "add-bed": {
      const input=config(),mixed=input.Adsorbents[0].MixFraction!==undefined;
      let number=input.Adsorbents.length+1; while(input.Adsorbents.some(bed=>bed.Name===`Bed ${number}`)) number++;
      if(!mixed) ensureSections();
      const bed={Name:`Bed ${number}`};
      if(mixed) bed.MixFraction=0;
      else input.ColumnSections.push({InterfaceLength:0},{Adsorbent:bed.Name,Length:.05});
      input.Adsorbents.push(bed);syncLengths();renderContent();scheduleSave();break;
    }
    case "remove-bed": {
      const input=config(),name=input.Adsorbents[index].Name,mixed=input.Adsorbents[0].MixFraction!==undefined;
      if(!mixed) {
        ensureSections(); const section=input.ColumnSections.findIndex(item=>item.Adsorbent===name);
        if(section>=0) {if(input.ColumnSections[section+1]?.InterfaceLength!==undefined) input.ColumnSections.splice(section,2);else input.ColumnSections.splice(Math.max(0,section-1),section===0?1:2);}
      }
      input.Adsorbents.splice(index,1);
      if(!input.Adsorbents.length) {delete input.Adsorbents;delete input.ColumnSections;}else syncLengths();
      renderContent();scheduleSave();break;
    }
    case "add-reaction": (config().Reactions??=[]).push({Phase:"PoreConcentration",Style:"GeneralPowerLaw",Site:0,Reactants:[0],Products:[Math.min(1,config().Components.length-1)],Kinetics:{forwardRateCoefficient:1,forwardActivationEnergy:0,equilibriumConstant:1,gibbsFreeEnergy:0},RateLimitTime:.0001});renderContent();scheduleSave();break;
    case "remove-reaction": config().Reactions.splice(index,1);renderContent();scheduleSave();break;
    case "add-participant": {
      const reaction=config().Reactions[index],side=button.dataset.side,used=[...reaction.Reactants,...reaction.Products].map(item=>typeof item==="number"?item:config().Components.findIndex(component=>component.Name.toLowerCase()===item.toLowerCase()));
      const available=config().Components.findIndex((_,i)=>!used.includes(i));if(available<0) throw new Error("Add another component before adding a reaction participant.");
      const offset=side==="Reactants"?reaction.Reactants.length:reaction.Reactants.length+reaction.Products.length;
      if(Array.isArray(reaction.Stoichiometry)) reaction.Stoichiometry.splice(offset,0,1);else if(reaction.Stoichiometry?.[side]) reaction.Stoichiometry[side].push(1);
      reaction[side==="Reactants"?"ForwardOrders":"BackwardOrders"]?.push(1);
      reaction[side].push(available);renderContent();scheduleSave();break;
    }
    case "remove-participant": { const reaction=config().Reactions[index],side=button.dataset.side,row=Number(button.dataset.row); if(Array.isArray(reaction.Stoichiometry)) reaction.Stoichiometry.splice(side==="Reactants"?row:reaction.Reactants.length+row,1); else reaction.Stoichiometry?.[side]?.splice(row,1); reaction[side].splice(row,1); reaction[side==="Reactants"?"ForwardOrders":"BackwardOrders"]?.splice(row,1); renderContent(); scheduleSave(); break; }
    case "reset": {
      const path=JSON.parse(button.dataset.path),control=document.getElementById(`field-${path.join("-")}`);
      const value=button.dataset.default === undefined ? undefined : JSON.parse(button.dataset.default);
      if(control?.tagName === "SELECT") {control.value=value??"";changeField(control);}
      else {setPath(config(),path,value);render();scheduleSave();}
      break;
    }
    case "reset-json": {
      const name=button.dataset.editor;
      const value=name==="simulation"?config():name==="grid"?(config().ColumnDistances||[]):Object.fromEntries(["ComponentParameters","Components"].filter(key=>config().Adsorbents[Number(name.split("-")[1])][key]!==undefined).map(key=>[key,config().Adsorbents[Number(name.split("-")[1])][key]]));
      $(`#json-${name}`).value=JSON.stringify(value,null,2);break;
    }
    case "apply-json": {
      const name=button.dataset.editor,value=JSON.parse($(`#json-${name}`).value);
      if(name==="simulation") { const checked=await api("/api/validate",{input:value}); simulation().input=checked.input; }
      else if(name==="grid") { if(!Array.isArray(value)||value.some(item=>typeof item!=="number"||!Number.isFinite(item))) throw new Error("Grid positions must be a JSON array of finite numbers."); if(value.length) config().ColumnDistances=value; else delete config().ColumnDistances; }
      else if(name.startsWith("bed-")) { if(!value||Array.isArray(value)||typeof value!=="object"||Object.keys(value).some(key=>!["ComponentParameters","Components"].includes(key))) throw new Error("Enter an object containing ComponentParameters and/or Components."); const bed=config().Adsorbents[Number(name.split("-")[1])]; delete bed.ComponentParameters; delete bed.Components; Object.assign(bed,value); }
      renderContent(); scheduleSave(); toast("Changes applied."); break;
    }
    case "download": await downloadFile(button.dataset.file); break;
    case "report": button.disabled=true; try { const result=await api(`/api/runs/${state.selectedRun}/report`,{}); await downloadFile(result.name); } finally { button.disabled=false; } break;
    case "export-video": await exportVideo(button);break;
    case "play": if(state.playing) {stopPlayback(); button.textContent="▶";} else {button.textContent="Ⅱ"; state.playing=setInterval(()=>{ const slider=$("#frame-slider")||$("#position-slider"); if(!slider) {stopPlayback();return;} slider.value=(Number(slider.value)+1)%(Number(slider.max)+1); slider.dispatchEvent(new Event("change",{bubbles:true})); },1000/(state.plotOptions.fps||20));} break;
  }
}

function renderNotes() {
  const help="Record experiment assumptions, provenance, and scientific references; notes are saved in the workspace alongside the simulation.";
  $("#content").innerHTML=card("Experiment notes",`<div class="field">${fieldHeader("experiment-notes","Notes",help)}<textarea id="experiment-notes" class="notes-editor" placeholder="Add assumptions, observations, or references…">${escapeHtml(simulation().notes||"")}</textarea></div><p class="helper">Notes are included when you save, import, or duplicate a workspace experiment.</p>`);
}
function applyTheme(theme,persist=true) {
  theme=theme==="light"?"light":"dark";state.theme=theme;
  document.documentElement.dataset.theme=theme;
  document.querySelector('meta[name="theme-color"]').content=theme==="light"?"#f4f7f9":"#10171e";
  if(persist) try {localStorage.setItem("ruptura-theme",theme);} catch {}
  if($("#plot-canvas")) updatePlot();
}

async function copyToClipboard(text) {
  if(navigator.clipboard?.writeText) {
    try {await navigator.clipboard.writeText(text);return;} catch { /* Try the user-gesture fallback below. */ }
  }
  const previous=document.activeElement, textarea=document.createElement("textarea");
  textarea.value=text; textarea.readOnly=true; textarea.tabIndex=-1;
  textarea.style.cssText="position:fixed;left:-9999px;top:0";
  document.body.append(textarea);
  try {
    textarea.focus({preventScroll:true}); textarea.select();
    if(!document.execCommand("copy")) throw new Error("Clipboard access was denied; allow clipboard access and try again.");
  } finally {textarea.remove();previous?.focus({preventScroll:true});}
}

const simulationMenu=document.createElement("div");
simulationMenu.id="simulation-menu"; simulationMenu.className="simulation-menu";
simulationMenu.setAttribute("role","menu"); simulationMenu.setAttribute("aria-label","Simulation actions");
simulationMenu.hidden=true; document.body.append(simulationMenu);
let simulationMenuTrigger=null;
function closeSimulationMenu(restoreFocus=false) {
  const trigger=simulationMenuTrigger;
  trigger?.setAttribute("aria-expanded","false"); simulationMenuTrigger=null;
  simulationMenu.hidden=true; simulationMenu.replaceChildren();
  if(restoreFocus&&trigger?.isConnected) trigger.focus({preventScroll:true});
}
function positionSimulationMenu() {
  if(simulationMenu.hidden||!simulationMenuTrigger) return;
  if(!simulationMenuTrigger.isConnected) {closeSimulationMenu();return;}
  const rect=simulationMenuTrigger.getBoundingClientRect(),list=$("#experiments").getBoundingClientRect();
  if(rect.bottom<=Math.max(0,list.top)||rect.top>=Math.min(innerHeight,list.bottom)||rect.right<=Math.max(0,list.left)||rect.left>=Math.min(innerWidth,list.right)) {closeSimulationMenu();return;}
  const width=simulationMenu.offsetWidth,height=simulationMenu.offsetHeight;
  simulationMenu.style.left=`${Math.max(8,Math.min(rect.right-width,innerWidth-width-8))}px`;
  simulationMenu.style.top=`${rect.bottom+height+8<=innerHeight?rect.bottom+5:Math.max(8,rect.top-height-5)}px`;
}
function openSimulationMenu(button,last=false) {
  if(simulationMenuTrigger===button&&!simulationMenu.hidden) {closeSimulationMenu(true);return;}
  closeSimulationMenu(); hideHelp();
  const source=state.workspace.simulations.find(item=>item.id===button.dataset.id);
  if(!source) return;
  const id=escapeHtml(source.id);
  simulationMenu.innerHTML=`<button role="menuitem" tabindex="-1" data-action="duplicate" data-id="${id}">Duplicate</button><button role="menuitem" tabindex="-1" data-action="copy-simulation-json" data-id="${id}">Copy simulation JSON</button><div class="menu-separator" role="separator"></div><button class="menu-delete" role="menuitem" tabindex="-1" data-action="delete-simulation" data-id="${id}" ${state.workspace.simulations.length===1?'disabled title="Keep at least one simulation"':''}>Delete</button>`;
  simulationMenu.hidden=false; simulationMenuTrigger=button; button.setAttribute("aria-expanded","true");
  positionSimulationMenu();
  const items=simulationMenu.querySelectorAll('button:not(:disabled)');
  items[last?items.length-1:0]?.focus({preventScroll:true});
}
document.addEventListener("pointerdown",event=>{
  if(!simulationMenu.contains(event.target)&&!event.target.closest('[data-action="simulation-menu"]')) closeSimulationMenu();
});
document.addEventListener("scroll",event=>{if(!simulationMenu.contains(event.target)) positionSimulationMenu();},true);
window.addEventListener("resize",positionSimulationMenu);
document.addEventListener("keydown",event=>{
  const trigger=event.target.closest('[data-action="simulation-menu"]');
  if(trigger&&["ArrowDown","ArrowUp"].includes(event.key)) {event.preventDefault();openSimulationMenu(trigger,event.key==="ArrowUp");return;}
  if(simulationMenu.hidden) return;
  if(event.key==="Escape") {event.preventDefault();closeSimulationMenu(true);return;}
  if(event.key==="Tab") {closeSimulationMenu(true);return;}
  if(!simulationMenu.contains(event.target)||!["ArrowDown","ArrowUp","Home","End"].includes(event.key)) return;
  event.preventDefault();
  const items=[...simulationMenu.querySelectorAll('button:not(:disabled)')],index=items.indexOf(document.activeElement);
  const next=event.key==="Home"?0:event.key==="End"?items.length-1:(index+(event.key==="ArrowUp"?-1:1)+items.length)%items.length;
  items[next]?.focus();
});

const helpTooltip = document.createElement("div");
helpTooltip.id = "field-tooltip"; helpTooltip.className = "field-tooltip";
helpTooltip.setAttribute("role", "tooltip"); helpTooltip.hidden = true;
document.body.append(helpTooltip);
let activeHelp = null, helpTimer;
function hideHelp() {
  clearTimeout(helpTimer);
  activeHelp?.removeAttribute("aria-describedby"); activeHelp = null; helpTooltip.hidden = true;
}
function showHelp(button) {
  clearTimeout(helpTimer);
  activeHelp?.removeAttribute("aria-describedby"); activeHelp = button;
  (button.closest("dialog") || document.body).append(helpTooltip);
  helpTooltip.textContent = button.dataset.help;
  helpTooltip.hidden = false; button.setAttribute("aria-describedby", helpTooltip.id);
  const rect = button.getBoundingClientRect(), width = helpTooltip.offsetWidth, height = helpTooltip.offsetHeight;
  helpTooltip.style.left = `${Math.max(12, Math.min(rect.right-width, innerWidth-width-12))}px`;
  helpTooltip.style.top = `${rect.bottom+height+10 <= innerHeight ? rect.bottom+8 : Math.max(8,rect.top-height-8)}px`;
}
function deferHideHelp() {
  clearTimeout(helpTimer);
  helpTimer=setTimeout(()=>{if(document.activeElement!==activeHelp) hideHelp();},150);
}
document.addEventListener("pointerover",event=>{const button=event.target.closest("[data-help]");if(button) showHelp(button);});
document.addEventListener("pointerout",event=>{if(event.target.closest("[data-help]")) deferHideHelp();});
document.addEventListener("focusin",event=>{const button=event.target.closest("[data-help]");if(button) showHelp(button);});
document.addEventListener("focusout",event=>{if(event.target.closest("[data-help]")) hideHelp();});
document.addEventListener("scroll",hideHelp,true);
document.addEventListener("keydown",event=>{if(event.key==="Escape"&&!helpTooltip.hidden) {event.preventDefault();hideHelp();}});
helpTooltip.addEventListener("pointerenter",()=>clearTimeout(helpTimer));
helpTooltip.addEventListener("pointerleave",deferHideHelp);
document.addEventListener("click",event=>{
  const button=event.target.closest("[data-action]");
  if(button?.tagName==="A") event.preventDefault();
  if(button?.dataset.action==="simulation-menu") {openSimulationMenu(button);return;}
  closeSimulationMenu(simulationMenu.contains(event.target));
  const help=event.target.closest("[data-help]"); if(help) {showHelp(help);return;}
  hideHelp();
  const tab=event.target.closest("[data-tab]"); if(tab) {selectTab(tab.dataset.tab);return;}
  if(button&&!button.disabled) action(button.dataset.action,button).catch(fail);
});
document.addEventListener("change",event=>{
  const target=event.target;
  if(target.id==="theme-select") {applyTheme(target.value);return;}
  if(target.id==="experiment-notes") {simulation().notes=target.value;scheduleSave();return;}
  if(target.dataset.path) {try{changeField(target);}catch(error){fail(error);} return;}
  if(target.id==="bed-mode") {arrangeBeds(target.value); render(); scheduleSave();return;}
  if(target.id==="run-selector") {stopPlayback();state.selectedRun=target.value;state.metadata=null;state.plotOptions={};renderContent();return;}
  if(target.id==="upload") {const file=target.files[0];if(file) file.text().then(text=>importWorkspace({document:JSON.parse(text)})).catch(fail);target.value="";return;}
  if(target.dataset.plot) {const key=target.dataset.plot;state.plotOptions[key]=target.type==="checkbox"?target.checked:key==="grid"?Number(target.value):target.value; if(["plot","metric","xAxis"].includes(key)) {stopPlayback();renderResultLayout();} updatePlot();return;}
  if(target.dataset.plotComponent!==undefined) {const i=Number(target.dataset.plotComponent);state.plotOptions.components=target.checked?[...state.plotOptions.components,i]:state.plotOptions.components.filter(value=>value!==i);updatePlot();return;}
  if(target.id==="video-fps") {const fps=Number(target.value);if(!Number.isInteger(fps)||fps<1||fps>60) {target.setCustomValidity("Enter a whole number from 1 to 60.");target.reportValidity();return;}target.setCustomValidity("");state.plotOptions.fps=fps;stopPlayback();if($("#play-button")) $("#play-button").textContent="▶";return;}
  if(["frame-slider","position-slider"].includes(target.id)) {state.plotOptions[target.id==="frame-slider"?"frame":"grid"]=Number(target.value);target.closest(".range-row").querySelector("label").innerHTML=profileSliderLabel();updatePlot();}
});
document.addEventListener("keydown",event=>{if(event.key==="Enter"&&event.target.id==="browse-path") action("browse-path").catch(fail);});
window.addEventListener("resize",()=>{if(window.Plotly) $$(".plot-canvas").forEach(canvas=>Plotly.Plots.resize(canvas));});

async function poll() {
  if(polling||!state.session) return; polling=true;
  try {
    const before=activeRun()?.status;
    state.runs=await api("/api/runs");
    const after=activeRun()?.status;
    if(before!==after) {header();if(["results","log"].includes(state.activeTab)) renderContent();}
    if(state.activeTab==="log") await updateLog();
  } catch(error) {$("#save-status").textContent="Server disconnected";}
  finally {polling=false;}
}
async function boot() {
  const data=await api("/api/bootstrap");
  Object.assign(state,{workspace:data.state,schema:data.schema,presets:data.presets,presetInfo:data.presetInfo,presetNotes:data.presetNotes,about:data.about,runs:data.runs,session:data.session,capabilities:data.capabilities,root:data.workspace,python:data.python});
  $("#workspace-path").textContent=data.workspace; $("#workspace-path").title=data.workspace;
  $("#save-status").textContent="All changes saved";
  let theme="dark";try {theme=localStorage.getItem("ruptura-theme")||"dark";} catch {}
  applyTheme(theme,false);render();renderAbout();
  if(data.capabilities.error) notice(data.capabilities.error,true);
  else if(data.missingDependencies.length) notice(`Plotting dependencies missing: ${data.missingDependencies.join(", ")}. Install ruptura[web] in ${data.python}.`,true);
  setInterval(poll,1000);
}
boot().catch(error=>{$("#content").innerHTML=`<section class="card"><div class="empty"><h2>Could not open the workspace</h2><p>${escapeHtml(error.message)}</p></div></section>`;fail(error);});
