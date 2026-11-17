"""Regression tests for PR revision groups, independent physics selection, and overview statistics."""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.skipif(shutil.which('node') is None, reason='Node.js not available for JS policy test')
def test_pr_overview_policy():
    template = (Path(__file__).resolve().parents[1] / 'perf/templates/dashboard.html').read_text()
    # Run the actual dashboard policy functions, without browser dependencies.
    script = template[template.index('function prRunKey('):template.index('function fillPRs(')]
    checks = r"""
const assert = require('node:assert/strict');
const base={pr:240,head_sha:'head',base_sha:'base',candidate_sha:'merge',attempt:1};
const cpu={...base,run_id:'1',rows:[{name:'slow',delta:24.9,pair:{status:'OK'}},
 {name:'fast',delta:-12,pair:{status:'OK'},scheduled:{status:'CHANGED'}},
 {name:'unchanged',delta:0,pair:{status:'UNVERIFIED'}},
 {name:'missing',delta:null,pair:{status:'NOT_CHECKED'}}]};
const oldPhysics={...base,run_id:'2',kind:'physics',conclusion:'success',rows:[{status:'PASSED'}]};
const newPhysics={...oldPhysics,run_id:'3',base_sha:'new-base',candidate_sha:'new-merge'};
const fluka={...cpu,run_id:'4',profile:'fluka',machine:'CERN',rows:[{name:'slow',delta:15,pair:{status:'OK'}}]};
const geant4={...cpu,run_id:'5',profile:'geant4',machine:'CERN',rows:[{name:'fast',delta:-20,pair:{status:'OK'}}]};
const oldFluka={...fluka,run_id:'6',rows:[{name:'old',delta:99}]};
const other={...fluka,run_id:'7',base_sha:'new-base',candidate_sha:'new-merge'};
const groups=benchmarkPRGroups([geant4,fluka,other,oldFluka,cpu]);
assert.equal(groups.length,2);assert.deepEqual(groups[0].components,[cpu,fluka,geant4]);
assert.equal(groups[0].rows.length,6);assert.equal(groups[0].rows.filter(x=>x.name==='slow').length,2);
assert.equal(groups[0].rows.find(x=>x.engine==='FLUKA').sourceRun,'4');
assert.equal(groups[0].rows.find(x=>x.engine==='FLUKA').sourceMachine,'CERN');
assert.equal(groups[1].rows.length,1);assert.equal(groups[1].candidate_sha,'new-merge');
assert.equal(benchmarkPROverview(groups[0]).slow,2);assert.equal(benchmarkPROverview(groups[0]).fast,2);
assert.equal(benchmarkPROverview(groups[0]).comparable,5);
assert.equal(benchmarkPROverview(groups[0]).worst.delta,24.9);
let chosen=defaultPRSelection([geant4,fluka,newPhysics,oldPhysics,cpu]);
assert.deepEqual(chosen.benchmark.components,[cpu,fluka,geant4]);assert.equal(chosen.validation,oldPhysics);
chosen=defaultPRSelection([other,geant4,fluka,newPhysics,oldPhysics,cpu],prRevisionKey(cpu),'3-1');
assert.deepEqual(chosen.benchmark.components,[cpu,fluka,geant4]);assert.equal(chosen.validation,newPhysics);
assert.equal(defaultPRSelection([other,cpu],prRevisionKey(other)).benchmark.candidate_sha,'new-merge');
assert(!samePRRevision(chosen.benchmark,chosen.validation));
assert.equal(defaultPRSelection([newPhysics,cpu]).validation,newPhysics);
assert.equal(defaultPRSelection([oldPhysics]).benchmark,undefined);
assert.equal(defaultPRSelection([cpu]).validation,undefined);
assert.equal(defaultPRSelection([]).validation,undefined);
assert.equal(prRunLabel(geant4),'Geant4');assert.equal(prRunLabel(cpu),'CPU');
assert.deepEqual(physicsPROverview(oldPhysics).status,['All passed','good']);
const bad={...oldPhysics,conclusion:'failure',rows:[{status:'PASSED'},{status:'FAILED'},
 {status:'ERROR'},{status:'NOT_RUN'},{status:'SKIPPED'}]};
const p=physicsPROverview(bad);assert.equal(p.passed,1);assert.equal(p.failed,2);
assert.equal(p.notRun,1);assert.equal(p.skipped,1);assert.equal(p.status[1],'bad');
assert.equal(physicsPROverview({...oldPhysics,conclusion:'failure'}).status[0],'Workflow failed');
assert.equal(physicsPROverview({...oldPhysics,rows:[]}).status[0],'Incomplete');
assert.equal(physicsPROverview({...oldPhysics,conclusion:'cancelled'}).status[0],'Incomplete');
const b=benchmarkPROverview(cpu);assert.equal(b.comparable,3);assert.equal(b.slow,1);
assert.equal(b.fast,1);assert.equal(b.flags,1);assert.equal(b.unverified,2);
assert.equal(b.worst.name,'slow');assert.equal(b.best.name,'fast');
const noSlow=benchmarkPROverview({...cpu,rows:[{name:'flat',delta:0,pair:{status:'OK'}}]});
assert.equal(noSlow.worst,null);assert.equal(noSlow.best,null);
assert.equal(noSlow.status[0],'No slowdowns ≥10%');
assert.equal(benchmarkPROverview({...cpu,rows:[{delta:null}]}).status[1],'warn');
"""
    subprocess.run(['node', '-e', script + checks], check=True, capture_output=True, text=True)
