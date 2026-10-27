"""Optional browser smoke test. Requires Playwright and its Chromium headless shell.

Run against a fresh server started with a disposable --workspace directory.
"""
from pathlib import Path
import json
import sys
from playwright.sync_api import sync_playwright, expect

import argparse
import atexit
import os
import subprocess
import tempfile
parser=argparse.ArgumentParser(description='Browser smoke check; launches a disposable Lab by default.')
parser.add_argument('--url', help='An existing server with a disposable workspace')
parser.add_argument('--launcher', type=Path, default=Path(__file__).resolve().parents[1]/'build/app/ruptura_lab')
parser.add_argument('--screenshots', type=Path, default=Path(tempfile.gettempdir())/'ruptura-web-screenshots')
args=parser.parse_args()
out=args.screenshots
out.mkdir(parents=True,exist_ok=True)
if not args.url:
    root=Path(__file__).resolve().parents[1]
    workspace_directory=tempfile.TemporaryDirectory(prefix='ruptura-browser-')
    process=subprocess.Popen([str(args.launcher),'--workspace',workspace_directory.name,'--port','0','--no-browser','--cli',str(root/'build/ruptura')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    def cleanup():
        process.terminate()
        try: process.wait(timeout=5)
        except subprocess.TimeoutExpired: process.kill(); process.wait()
        workspace_directory.cleanup()
    atexit.register(cleanup)
    line=process.stdout.readline()
    if not line.startswith('Ruptura Lab: '): raise RuntimeError(line+process.stdout.read())
    args.url=line.removeprefix('Ruptura Lab: ').strip()
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True)
    page=browser.new_page(viewport={'width':1512,'height':1100},device_scale_factor=1)
    errors=[]
    page.on('pageerror',lambda error:errors.append(str(error)))
    page.goto(args.url,wait_until='networkidle')
    page.context.grant_permissions(['clipboard-read','clipboard-write'],origin=args.url)
    expect(page.locator('#about-page')).to_be_visible()
    expect(page.locator('.sidebar')).to_be_visible()
    assert page.locator('.library-column, #library-editor-page, #workspace-page').count()==0
    assert page.locator('.brand-sub').count()==0
    assert page.locator('.brand-logo').evaluate('(img)=>img.complete&&img.naturalWidth===200')
    page.locator('[data-action="license"]').click()
    expect(page.locator('.license-text')).to_contain_text('Permission is hereby granted')
    page.locator('[data-action="close-modal"]').click()
    page.locator('#workspace-button').click()
    expect(page.locator('#simulation-page')).to_be_visible()
    assert page.locator('.workspace').evaluate('(node)=>getComputedStyle(node).marginLeft')=='250px'
    expect(page.locator('#experiments .experiment')).to_have_count(4)
    assert page.locator('.experiment-row').count()==page.locator('.experiment-row .simulation-menu-button').count()
    page.screenshot(path=str(out/'30-simulation-sidebar-dark.png'),full_page=True)
    boot=page.request.get(args.url+'/api/bootstrap').json()
    source=boot['state']['simulations'][0]
    gas_id=source['id']
    assert 'library' not in boot['state'] and 'library' not in boot['schema']

    def select_simulation(identifier,tab='setup'):
        page.locator(f'#experiments [data-action="select-simulation"][data-id="{identifier}"]').click()
        expect(page.locator('#simulation-page')).to_be_visible()
        page.locator(f'[data-tab="{tab}"]').click()

    def open_menu(identifier):
        page.locator(f'[data-action="simulation-menu"][data-id="{identifier}"]').click()
        menu=page.locator('#simulation-menu')
        expect(menu).to_be_visible()
        return menu

    def scientific_input():
        page.locator('[data-tab="json"]').click()
        return json.loads(page.locator('#json-simulation').input_value())

    def export_workspace():
        with page.expect_download() as saved:
            page.locator('[data-action="export-state"]').click()
        return json.loads(Path(saved.value.path()).read_text())

    def edit(selector,value):
        page.locator(selector).fill(str(value))
        page.locator(selector).blur()

    # The examples, About, help, and notes remain available beside the list.
    select_simulation('default-paper-t15')
    expect(page.locator('#field-Temperature')).to_have_value('298')
    assert page.locator('.tabs').bounding_box()['height']>=60
    assert page.locator('.tabs [data-tab="setup"]').evaluate('(node)=>getComputedStyle(node).fontSize')=='16px'
    page.locator('[data-tab="json"]').click()
    original_example=json.loads(page.locator('#json-simulation').input_value())
    assert original_example['Geometry']['Type']=='PackedBed'
    page.locator('[data-action="apply-json"][data-editor="simulation"]').click()
    expect(page.locator('#toast')).to_have_text('Changes applied.')
    # Both reset and selecting Default must serialize the C++-required type.
    for use_reset in (True,False):
        page.locator('[data-tab="column"]').click()
        if use_reset:
            page.locator('[data-action="reset"][data-path=\'["Geometry","Type"]\']').click()
        else:
            page.locator('#field-Geometry-Type').select_option('')
        expect(page.locator('#field-Geometry-Type')).to_have_value('PackedBed')
        assert scientific_input()==original_example
        page.locator('[data-action="apply-json"][data-editor="simulation"]').click()
        expect(page.locator('#toast')).to_have_text('Changes applied.')
        page.locator('#validate-button').click()
        expect(page.locator('#notice')).to_contain_text('Input validated')
    page.screenshot(path=str(out/'50-expanded-tabs-json-validated.png'),full_page=True)
    page.locator('[data-tab="setup"]').click()
    help_button=page.locator('.field').filter(has=page.locator('#field-Temperature')).locator('[data-help]')
    help_button.hover()
    expect(page.locator('#field-tooltip')).to_contain_text('temperature')
    page.locator('[data-tab="notes"]').click()
    assert '10.1007/s10450-013-9593-5' in page.locator('#experiment-notes').input_value()
    assert page.locator('#experiment-notes').bounding_box()['height']>=480
    assert page.locator('.card').first.evaluate('(node)=>getComputedStyle(node).borderTopWidth')=='2px'
    page.screenshot(path=str(out/'40-notes-expanded.png'),full_page=True)
    page.locator('.brand').click()
    expect(page.locator('#about-page')).to_be_visible()
    select_simulation(gas_id)

    # Three-dot actions operate on their row and support the keyboard.
    menu=open_menu(gas_id)
    expect(menu.get_by_role('menuitem')).to_have_text(['Duplicate','Copy simulation JSON','Delete'])
    page.keyboard.press('End')
    expect(menu.get_by_role('menuitem',name='Delete',exact=True)).to_be_focused()
    page.keyboard.press('Escape')
    expect(menu).to_be_hidden()
    expect(page.locator(f'[data-action="simulation-menu"][data-id="{gas_id}"]')).to_be_focused()
    open_menu(gas_id).get_by_role('menuitem',name='Duplicate',exact=True).click()
    expect(page.locator('#project-title')).to_have_text('untitled_1')
    clone_id=page.locator('#experiments .experiment.active').get_attribute('data-id')
    page.locator('[data-tab="notes"]').click()
    edit('#experiment-notes','An independent model, with my observations.')
    page.locator('[data-tab="setup"]').click()
    edit('#field-Temperature',305)
    page.locator('[data-tab="components"]').click()
    page.locator('[data-action="add-site"][data-index="0"][data-kind="PhysisorptionSites"]').click()
    page.locator('[data-action="add-site"][data-index="0"][data-kind="ChemisorptionSites"]').click()
    expect(page.locator('#field-Components-0-CarrierGas')).to_be_checked()
    expect(page.locator('#field-Components-0-PhysisorptionSites-0-Parameters-0')).to_have_value('1')
    expect(page.locator('#field-Components-0-ChemisorptionSites-0-MaximumLoading')).to_have_value('0.8')
    edit('#field-Components-1-MolecularWeight',0.049)
    edit('#field-Components-1-HeatOfAdsorption',42000)
    edit('#field-Components-1-PhysisorptionSites-0-Parameters-0',3)
    page.screenshot(path=str(out/'31-independent-components.png'),full_page=True)
    page.locator('[data-action="add-component"]').click()
    expect(page.locator('#modal')).to_be_hidden()
    expect(page.locator('#component-count')).to_have_text('3')
    page.locator('[data-action="add-site"][data-index="2"][data-kind="ChemisorptionSites"]').click()
    expect(page.locator('#field-Components-2-ChemisorptionSites-0-MaximumLoading')).to_have_value('0.8')
    page.locator('[data-action="remove-component"][data-index="2"]').click()
    expect(page.locator('#component-count')).to_have_text('2')
    page.locator('[data-tab="column"]').click()
    edit('#field-ColumnLength',0.2)
    configured=scientific_input()
    assert configured['Components'][1]['MolecularWeight']==0.049
    assert configured['Components'][1]['HeatOfAdsorption']==42000
    assert configured['Components'][1]['PhysisorptionSites'][0]['Parameters'][0]==3
    assert configured['Temperature']==305 and configured['ColumnLength']==0.2
    select_simulation(gas_id)
    assert scientific_input()==source['input']
    # Copying an unselected row leaves selection unchanged and copies its full input.
    open_menu(clone_id).get_by_role('menuitem',name='Copy simulation JSON',exact=True).click()
    expect(page.locator('#toast')).to_contain_text('copied to clipboard')
    assert json.loads(page.evaluate('navigator.clipboard.readText()'))==configured
    assert page.locator('#experiments .experiment.active').get_attribute('data-id')==gas_id
    expect(page.locator('#save-status')).to_have_text('All changes saved')
    page.reload(wait_until='networkidle')
    expect(page.locator('#about-page')).to_be_visible()
    select_simulation(clone_id)
    assert scientific_input()==configured
    page.locator('[data-tab="notes"]').click()
    expect(page.locator('#experiment-notes')).to_have_value('An independent model, with my observations.')

    # Layered columns and reactions are edited locally, without a library picker.
    page.locator('#new-simulation-button').click()
    quick_names=page.locator('.preset-button strong').all_text_contents()
    assert 'Gas breakthrough' in quick_names and 'Liquid & pH' in quick_names
    assert sum('Monolith' in name for name in quick_names)==1
    page.locator('[data-action="create"][data-preset="multibed"]').click()
    expect(page.locator('#project-title')).to_have_text('untitled_2')
    layered_id=page.locator('#experiments .experiment.active').get_attribute('data-id')
    page.locator('[data-tab="column"]').click()
    expect(page.locator('#bed-mode')).to_have_value('layered')
    edit('#field-ColumnSections-0-Length',0.07)
    page.locator('[data-action="add-bed"]').click()
    expect(page.locator('#modal')).to_be_hidden()
    layered=scientific_input()
    assert [section['Length'] for section in layered['ColumnSections'] if 'Length' in section]==[0.07,0.05,0.05]
    assert abs(layered['ColumnLength']-0.17)<1e-12
    page.locator('[data-tab="column"]').click()
    page.locator('[data-action="remove-bed"][data-index="2"]').click()
    page.locator('#json-bed-0').fill(json.dumps({'ComponentParameters':{'CO2':{'HeatOfAdsorption':25000,'PhysisorptionSites':[{'Type':'Langmuir','Parameters':[4,0.00002]}]}}}))
    page.locator('[data-action="apply-json"][data-editor="bed-0"]').click()
    expect(page.locator('#toast')).to_have_text('Changes applied.')
    page.locator('[data-tab="components"]').click()
    page.locator('[data-action="add-component"]').click()
    page.locator('[data-tab="reactions"]').click()
    page.locator('[data-action="add-reaction"]').click()
    expect(page.locator('#modal')).to_be_hidden()
    edit('#field-Reactions-0-Stoichiometry-Reactants-0',2)
    edit('#field-Reactions-0-ForwardOrders-0',0.5)
    page.locator('[data-action="add-participant"][data-side="Reactants"]').click()
    reacted=scientific_input()
    assert reacted['Reactions'][0]['Stoichiometry']['Reactants']==[2,1]
    assert reacted['Reactions'][0]['ForwardOrders']==[0.5,1]
    assert reacted['Adsorbents'][0]['ComponentParameters']['CO2']['HeatOfAdsorption']==25000
    page.locator('[data-tab="reactions"]').click()
    page.locator('[data-action="remove-reaction"]').click()
    select_simulation(gas_id)
    assert scientific_input()==source['input']
    open_menu(layered_id).get_by_role('menuitem',name='Delete',exact=True).click()
    expect(page.locator(f'#experiments .experiment[data-id="{layered_id}"]')).to_have_count(0)
    assert page.locator('#experiments .experiment.active').get_attribute('data-id')==gas_id

    # Pasted inputs retain native validation and the local run/plot pipeline.
    page.locator('#new-simulation-button').click()
    page.locator('#new-simulation-json').fill('not JSON')
    page.locator('[data-action="create-json"]').click()
    expect(page.locator('#new-json-error')).to_contain_text('not valid JSON')
    short=dict(source['input'],SimulationTime=1e-5,TimeStep=1e-6,WriteInterval=1e-6,NumberOfGridPoints=3,energyBalance=True,
        wallDensity=7800,gasThermalConductivity=0.03,wallThermalConductivity=16,heatTransferGasSolid=20,
        heatTransferGasWall=20,heatTransferWallExternal=5,heatCapacityGas=1000,heatCapacitySolid=1000,heatCapacityWall=500)
    if short.get('BreakthroughIntegrator')=='CVODE': short.pop('TimeStep')
    short['Geometry']={**short['Geometry'],'InternalDiameter':0.01,'OuterDiameter':0.012}
    page.locator('#new-simulation-json').fill(json.dumps(short))
    page.locator('[data-action="create-json"]').click()
    expect(page.locator('#project-title')).to_have_text('untitled_3')
    new_id=page.locator('#experiments .experiment.active').get_attribute('data-id')
    page.locator('#run-button').click()
    expect(page.locator('[data-tab="log"]')).to_have_class('active')
    expect(page.locator('#run-log')).to_contain_text('Working directory',timeout=10000)
    expect(page.locator('#metrics')).to_contain_text('completed',timeout=30000)
    page.locator('[data-tab="results"]').click()
    expect(page.locator('#plot-canvas .main-svg').first).to_be_visible(timeout=30000)
    expect(page.locator('#temperature-canvas .main-svg').first).to_be_visible(timeout=30000)
    assert page.locator('#temperature-canvas').bounding_box()['y']>page.locator('#plot-canvas').bounding_box()['y']
    assert page.locator('[data-plot="showCarrier"] ').count()==0
    page.locator('[data-plot-component="0"]').check()
    page.wait_for_function('document.querySelector("#plot-canvas").data.some(trace=>trace.name.startsWith("N2"))')
    page.locator('[data-action="settings"]:visible').click()
    page.locator('#theme-select').select_option('light')
    page.locator('[data-action="close-modal"]').click()
    page.wait_for_function('document.querySelector("#plot-canvas").layout.paper_bgcolor==="#ffffff"')
    page.wait_for_function('document.querySelector("#temperature-canvas").layout.paper_bgcolor==="#ffffff"')
    page.screenshot(path=str(out/'41-energy-results-stacked.png'),full_page=True)
    page.locator('#plot-plot').select_option('profile')
    expect(page.locator('#frame-slider')).to_be_visible()
    assert page.locator('#position-slider').count()==0
    page.wait_for_function('document.querySelector("#plot-canvas").layout?.xaxis.title.text==="Position [m]"')
    page.locator('#plot-xAxis').select_option('time')
    expect(page.locator('#position-slider')).to_be_visible()
    assert page.locator('#frame-slider').count()==0
    page.locator('#position-slider').fill('0')
    page.locator('#position-slider').dispatch_event('change')
    page.wait_for_function('document.querySelector("#plot-canvas").layout?.xaxis.title.text==="Time [min]"')
    expect(page.locator('.range-row label')).to_contain_text('POSITION 1 / 4')
    edit('#video-fps',12)
    with page.expect_download(timeout=60000) as movie:
        page.locator('[data-action="export-video"]').click()
    assert b'ftyp' in Path(movie.value.path()).read_bytes()[:32]
    inspected=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-of','json',str(movie.value.path())]))
    assert inspected['streams'][0]['r_frame_rate']=='12/1'
    assert inspected['streams'][0]['nb_frames']=='4'
    page.screenshot(path=str(out/'42-column-time-viewer.png'),full_page=True)
    page.locator('#plot-xAxis').select_option('position')
    expect(page.locator('#frame-slider')).to_be_visible()
    page.locator('#frame-slider').fill('0')
    page.locator('#frame-slider').dispatch_event('change')
    expect(page.locator('.range-row label')).to_contain_text('FRAME 1 /')
    with page.expect_download() as notebook:
        page.locator('[data-action="download"][data-file="analysis.ipynb"]').first.click()
    assert json.loads(Path(notebook.value.path()).read_text())['nbformat']==4
    page.reload(wait_until='networkidle')
    expect(page.locator('#about-page')).to_be_visible()
    assert page.evaluate('getComputedStyle(document.documentElement).colorScheme')=='light'

    # Workspace paste/export contains independent inputs, notes, and unique IDs.
    saved=export_workspace()
    assert saved['Version']==3 and 'library' not in saved
    assert all(not {'bindings','libraryDefaults','templateId'}.intersection(item) for item in saved['simulations'])
    page.locator('[data-action="import"]').click()
    page.locator('#import-workspace-json').fill('not JSON')
    page.locator('[data-action="import-json"]').click()
    expect(page.locator('#import-json-error')).to_contain_text('not valid JSON')
    page.locator('#import-workspace-json').fill(json.dumps(saved))
    page.locator('[data-action="import-json"]').click()
    expect(page.locator('#modal')).to_be_hidden()
    imported=export_workspace()
    assert len(imported['simulations'])==2*len(saved['simulations'])
    assert len({item['id'] for item in imported['simulations']})==len(imported['simulations'])
    assert [item['input'] for item in imported['simulations'][len(saved['simulations']):]]==[item['input'] for item in saved['simulations']]
    assert [item['notes'] for item in imported['simulations'][len(saved['simulations']):]]==[item['notes'] for item in saved['simulations']]
    assert next(item for item in imported['simulations'] if item['id']==new_id).get('lastRunId')

    # Narrow layouts keep the simulation list and menu usable without page overflow.
    select_simulation(clone_id)
    page.set_viewport_size({'width':390,'height':900})
    page.locator('[data-tab="setup"]').click()
    assert page.evaluate('document.body.scrollWidth<=innerWidth'),page.evaluate('document.body.scrollWidth')
    expect(page.locator('#experiments')).to_be_visible()
    menu=open_menu(clone_id)
    bounds=menu.bounding_box()
    assert bounds['x']>=0 and bounds['x']+bounds['width']<=390
    page.screenshot(path=str(out/'33-simulation-sidebar-mobile.png'))
    page.keyboard.press('Escape')
    # Removing the selected simulation chooses another one and keeps settings visible.
    open_menu(clone_id).get_by_role('menuitem',name='Delete',exact=True).click()
    expect(page.locator('#simulation-page')).to_be_visible()
    expect(page.locator('#project-title')).not_to_have_text('untitled_1')
    assert not errors,errors
    browser.close()
print('Browser smoke passed: independent sidebar simulations, forms, migration-compatible import/export, menus, notes, local run/plots, themes, and mobile layout.')
