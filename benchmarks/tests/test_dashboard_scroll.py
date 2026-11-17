"""Regression checks for the benchmark history's rightmost initial viewport."""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.skipif(shutil.which('node') is None, reason='Node.js not available')
def test_timeline_starts_at_latest_and_respects_manual_scrolling():
    template = (Path(__file__).resolve().parents[1] /
                'perf/templates/dashboard.html').read_text()
    source = template.split('function initialiseTimelineScroll() {', 1)[1]
    source = 'function initialiseTimelineScroll() {' + source.split('\ninit();', 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const area={scrollWidth:1250,scrollLeft:0};
const $=id=>{assert.equal(id,'timeline');return area;};
let observer, observed;
class MutationObserver {
  constructor(cb){this.cb=cb;observer=this;}
  observe(element,opts){observed=element;assert.deepEqual(opts,{childList:true});}
}
const requestAnimationFrame=fn=>fn();
initialiseTimelineScroll();
assert.equal(observed,area);
assert.equal(area.scrollLeft,1250,'Initial heatmap aligns to right edge');
area.scrollLeft=410;
assert.equal(area.scrollLeft,410,'Manual scroll is not overridden');
area.scrollWidth=2350;
area.scrollLeft=0;
observer.cb();
assert.equal(area.scrollLeft,2350,'Newly rendered trend/heatmap aligns to right edge');
"""
    subprocess.run(['node', '-e', source + '\n' + checks],
                   check=True, text=True, capture_output=True)
