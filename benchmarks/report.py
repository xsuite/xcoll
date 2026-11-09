"""Create an offline interactive dashboard for GitHub Pages or local viewing.

Reads any number of pytest-benchmark JSON files; no network or JS libraries.

No application server or third-party Python dependencies required.

    python benchmarks/report.py .benchmarks/release.json \
        .benchmarks/gpu.json -o .benchmarks/dashboard.html

Open the resulting HTML in any browser. Inputs may come from different Xsuite
projects, but comparisons are restricted to identical scenario/protocol/testbed.
"""

import argparse
import json
from pathlib import Path

COMPARE_FIELDS = (
    "protocol", "particles", "turns", "crystals", "beam", "plane",
    "interpolation_m", "physics_case", "geometry", "tilt_radians",
    "sequence_sha256", "collimator_database_sha256",
)


def read_runs(paths):
    runs = []
    for path in paths:
        path = Path(path)
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data.get("benchmarks"), list):
            raise ValueError(f"Not a pytest-benchmark JSON result: {path}")
        machine = data.get("machine_info") or {}
        cpu = machine.get("cpu") or {}
        processor = machine.get("processor") or cpu.get("brand_raw") or cpu.get("brand") or "unknown CPU"
        host = machine.get("node", "unknown-host")
        testbed = data.get("testbed", "unlabelled-local")
        # CI hostnames are disposable: comparing by them destroys history.
        # Compare on labelled testbed + hardware. Add hostname only for local
        # runs (where a developer's machine is assumed to be persistent).
        machine_key = f"{testbed} | {machine.get('system', '?')} | {processor}"
        if not str(testbed).startswith("github-"):
            machine_key += f" | {host}"
        ci = data.get("commit_info") or {}
        entries = []
        for b in data["benchmarks"]:
            stats = b.get("stats") or {}
            if stats.get("median") is None:
                continue
            extra = b.get("extra_info") or {}
            entries.append({
                "name": b.get("fullname") or b.get("name") or "unnamed",
                "group": b.get("group") or "Ungrouped",
                "median": float(stats["median"]),
                "stddev": (float(stats["stddev"]) if stats.get("stddev") is not None else None),
                "rounds": (int(stats["rounds"]) if stats.get("rounds") is not None else None),
                "extra": extra,
                "signature": [extra.get(k) for k in COMPARE_FIELDS],
            })
        runs.append({
            "id": path.name,
            "datetime": data.get("datetime") or "",
            "sha": (data.get("ci") or {}).get("source_sha") or ci.get("id") or "unknown",
            "branch": (data.get("ci") or {}).get("source_branch") or ci.get("branch") or "unknown",
            "machine": machine_key,
            "versions": data.get("xsuite_versions") or {},
            "entries": entries,
        })
    return sorted(runs, key=lambda r: (r["datetime"], r["id"]))


PAGE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Xsuite performance</title>
<style>
:root{color-scheme:light;font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
 --ink:#12223b;--sub:#53637b;--border:#dce4ed;--surface:#f4f7fb;--blue:#2663c7}
body{margin:0;background:white;color:var(--ink)}
main{max-width:1280px;margin:auto;padding:28px 24px 60px}
header{display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:15px}
h1{font-size:1.65rem;margin:0}.sub{color:var(--sub);font-size:.87rem}
.controls{display:flex;gap:12px;flex-wrap:wrap;margin:22px 0 16px}
label{display:grid;gap:5px;color:var(--sub);font-size:.85rem}
select,input{min-width:180px;max-width:360px;padding:8px;border:1px solid #b8c5d5;border-radius:7px;background:white;color:var(--ink)}
.stats{display:flex;gap:10px;margin:10px 0 17px;flex-wrap:wrap}
.stat{background:var(--surface);border:1px solid var(--border);padding:12px 18px;border-radius:9px;min-width:130px}
.stat strong{display:block;font-size:1.55rem}.stat span{font-size:.78rem;color:var(--sub)}
.grid{display:grid;grid-template-columns:minmax(0,1.3fr) minmax(300px,.85fr);gap:16px}
.panel{border:1px solid var(--border);border-radius:10px;overflow:hidden;background:white}
.panel h2{font-size:1.05rem;margin:0;padding:15px 16px;background:var(--surface);border-bottom:1px solid var(--border)}
.scroll{max-height:690px;overflow:auto}
table{border-collapse:collapse;width:100%;font-size:.84rem}th,td{padding:9px 11px;border-bottom:1px solid #e7edf4;text-align:left}th{background:var(--surface);position:sticky;top:0;color:var(--sub)}
tr.selectable{cursor:pointer}tr.selectable:hover,tr.selected{background:#eaf1ff}
td.name{overflow-wrap:anywhere;max-width:370px}td.num{font-variant-numeric:tabular-nums;white-space:nowrap}
.warn{color:#ae4318;font-weight:650}.good{color:#16754b;font-weight:650}.muted{color:var(--sub)}
.detail{padding:16px}.detail h3{margin:0 0 4px;font-size:.98rem;overflow-wrap:anywhere}
svg{max-width:100%;height:auto}.note{font-size:.8rem;color:var(--sub);margin:11px 0}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:var(--surface);border:1px solid var(--border);border-radius:7px;padding:12px;font-size:.73rem}
@media(max-width:850px){.grid{grid-template-columns:1fr}main{padding:18px 12px}select,input{min-width:130px}}
</style></head><body><main>
<header><div><h1>Xsuite performance dashboard</h1><div class="sub">Offline viewer · pytest-benchmark results · identical-machine comparisons only</div></div>
<div class="sub" id="subtitle"></div></header>
<div class="controls">
<label>Machine / testbed<select id="machine"></select></label>
<label>Latest run<select id="latest"></select></label>
<label>Filter benchmarks<input id="search" placeholder="Everest, geometry, LHC…"></label>
</div>
<div class="stats" id="stats"></div>
<div class="grid"><section class="panel"><h2>Latest run vs previous comparable run</h2>
<div class="scroll"><table><thead><tr><th>Benchmark</th><th>Median</th><th>Change</th><th>Physics</th><th>Rounds</th></tr></thead>
<tbody id="rows"></tbody></table></div></section>
<section class="panel"><h2>Benchmark history</h2><div class="detail"><h3 id="detailName">Select a benchmark</h3>
<div class="sub" id="detailSubtitle"></div><div id="plot"></div>
<div class="note" id="detailNote"></div><pre id="details"></pre></div></section></div>
<p class="note">Physical checks are indicative comparisons of recorded observables, not full physics validation. Missing observables are labelled unverified. Changes are informational, not statistical release gates. Warning: ≥10% slowdown; improvement: ≥10% speed-up. Historical points must have identical workload signatures and machine identities. Compare repeated runs before concluding a regression.</p>
</main><script>
const runs = __RUN_DATA__;
const machines = [...new Set(runs.map(x=>x.machine))].sort();
const $=id=>document.getElementById(id);
let selected=null;
function el(tag,txt,cls){let e=document.createElement(tag);if(txt!==undefined)e.textContent=String(txt);if(cls)e.className=cls;return e}
function fill(select,values,label){select.replaceChildren();values.forEach(x=>{let o=el('option',label(x));o.value=x;select.append(o)})}
function fmt(s){if(!Number.isFinite(s))return '—';return s<.001?(s*1e6).toFixed(1)+' µs':s<1?(s*1e3).toFixed(2)+' ms':s.toFixed(3)+' s'}
function pct(a,b){return (a/b-1)*100}
function diffClass(d){return d>=10?'warn':d<=-10?'good':''}
// Physical results are screening checks. Use compare_physics.py for the
// authoritative per-observable comparison and exit code.
const physLimits={'lost_fraction':.01,'survivors_with_kick_gt_0.1urad':.02,
'channeling_proxy_fraction_of_survivors':.02,
'reference_collimator_losses_per_particle':.02,'reference_aperture_losses_per_particle':.02};
function physicals(e){let v=e.extra||{},p=Number(v.particles),out={};
for(let k of ['lost_fraction','survivors_with_kick_gt_0.1urad','channeling_proxy_fraction_of_survivors'])
if(typeof v[k]==='number')out[k]=v[k];
if(p>0)for(let k of ['reference_collimator_losses','reference_aperture_losses'])
if(typeof v[k]==='number')out[k+'_per_particle']=v[k]/p;
return out}
function physicalStatus(entry,baseline){if(!baseline)return '—';
let x=physicals(entry),y=physicals(baseline.entry);
let keys=Object.keys(x).sort(),other=Object.keys(y).sort();
if(!keys.length||JSON.stringify(keys)!==JSON.stringify(other))return 'unverified';
return keys.some(k=>!Number.isFinite(x[k])||!Number.isFinite(y[k])||
Math.abs(x[k]-y[k])>physLimits[k])?'changed':'ok'}
function sig(e){return JSON.stringify(e.signature)}
function getBaseline(run,e){let i=runs.indexOf(run);for(let k=i-1;k>=0;k--){let r=runs[k];if(r.machine!==run.machine)continue;let m=r.entries.find(x=>x.name===e.name && sig(x)===sig(e));if(m)return {run:r,entry:m}}return null}
function currentRun(){return runs.find(x=>x.id===$('latest').value && x.machine===$('machine').value)}
function machineChange(){let available=runs.filter(x=>x.machine===$('machine').value);fill($('latest'),available.map(x=>x.id).reverse(),x=>{let r=available.find(y=>y.id===x);return (r.datetime||r.id).slice(0,19)+' · '+r.branch+' · '+r.sha.slice(0,8)});render()}
function render(){let run=currentRun();if(!run)return;let query=$('search').value.toLowerCase();let entries=run.entries.filter(x=>x.name.toLowerCase().includes(query)||x.group.toLowerCase().includes(query));entries.sort((a,b)=>{let aa=getBaseline(run,a),bb=getBaseline(run,b);let da=aa?pct(a.median,aa.entry.median):-9999,db=bb?pct(b.median,bb.entry.median):-9999;return db-da});
$('rows').replaceChildren();let slow=0,fast=0,missing=0;
entries.forEach(e=>{let base=getBaseline(run,e);let d=base?pct(e.median,base.entry.median):null;if(d===null)missing++;else if(d>=10)slow++;else if(d<=-10)fast++;
let tr=el('tr');tr.className='selectable'+(selected===e.name?' selected':'');let n=el('td',e.name,'name');let t=el('td',fmt(e.median),'num');let c=el('td',d===null?'no baseline':(d>0?'+':'')+d.toFixed(1)+'%',d===null?'muted':diffClass(d));let phy=physicalStatus(e,base);let phytd=el('td',phy,phy==='changed'?'warn':phy==='ok'?'good':'muted');let rounds=el('td',e.rounds===null?'—':e.rounds,'num');tr.append(n,t,c,phytd,rounds);tr.addEventListener('click',()=>{selected=e.name;render()});$('rows').append(tr)});
$('stats').replaceChildren();for(const [label,value] of [['Results',entries.length],['Slowdowns ≥10%',slow],['Speed-ups ≥10%',fast],['No baseline',missing]]){let box=el('div',undefined,'stat');box.append(el('strong',value),el('span',label));$('stats').append(box)}
$('subtitle').textContent='Branch '+run.branch+' · '+run.sha.slice(0,12)+' · '+run.datetime;
if(!entries.some(e=>e.name===selected))selected=entries.length?entries[0].name:null;
if(selected){let e=entries.find(x=>x.name===selected);draw(run,e)}else{$('plot').replaceChildren();$('details').textContent=''}
}
function draw(run,e){$('detailName').textContent=e.name;$('detailSubtitle').textContent=e.group;
let candidates=runs.filter(x=>x.machine===run.machine && runs.indexOf(x)<=runs.indexOf(run)).map(r=>{let b=r.entries.find(x=>x.name===e.name&&sig(x)===sig(e));return b?{r,b}:null}).filter(Boolean);
const w=480,h=240,pad=50;let mn=Math.min(...candidates.map(x=>x.b.median)),mx=Math.max(...candidates.map(x=>x.b.median));if(mx===mn){mn*=.9;mx*=1.1}let n=candidates.length;
const svg=document.createElementNS('http://www.w3.org/2000/svg','svg');svg.setAttribute('viewBox',`0 0 ${w} ${h}`);svg.setAttribute('role','img');svg.setAttribute('aria-label','Historical median runtime by run');
function node(tag,attrs){let x=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))x.setAttribute(k,v);svg.append(x);return x}
for(let i=0;i<5;i++){let val=mn+(mx-mn)*i/4;let y=h-pad-(h-2*pad)*i/4;node('line',{x1:pad,y1:y,x2:w-14,y2:y,stroke:'#dce4ed'});let tx=node('text',{x:pad-7,y:y+4,'text-anchor':'end',fill:'#53637b','font-size':'11'});tx.textContent=fmt(val)}
const pts=candidates.map((x,i)=>{let xx=pad+(w-pad-24)*(n===1?0.5:i/(n-1));let yy=h-pad-(x.b.median-mn)/(mx-mn)*(h-2*pad);return [xx,yy]});
node('polyline',{points:pts.map(x=>x.join(',')).join(' '),fill:'none',stroke:'#2663c7','stroke-width':'2.5'});
pts.forEach(([x,y],i)=>{let circle=node('circle',{cx:x,cy:y,r:4,fill:'#2663c7'});let t=document.createElementNS('http://www.w3.org/2000/svg','title');t.textContent=candidates[i].r.sha.slice(0,8)+' · '+fmt(candidates[i].b.median);circle.append(t)});
for(const [x,txt] of [[pad,candidates[0].r.sha.slice(0,7)],[w-14,candidates[n-1].r.sha.slice(0,7)]]){let t=node('text',{x,y:h-16,'text-anchor':x===pad?'start':'end',fill:'#53637b','font-size':'11'});t.textContent=txt}
$('plot').replaceChildren(svg);$('detailNote').textContent=candidates.length+' compatible run(s) on '+run.machine;
let base=getBaseline(run,e);$('details').textContent=JSON.stringify({latest_median:fmt(e.median),previous_median:base?fmt(base.entry.median):null,change_percent:base?Number(pct(e.median,base.entry.median).toFixed(2)):null,rounds:e.rounds,standard_deviation:e.stddev,physics_status:physicalStatus(e,base),physics_current:physicals(e),physics_previous:base?physicals(base.entry):null,extra_info:e.extra,versions:run.versions},null,2)
}
fill($('machine'),machines,x=>x);if(machines.length){$('machine').value=runs[runs.length-1].machine;machineChange()}
$('machine').addEventListener('change',machineChange);$('latest').addEventListener('change',()=>{selected=null;render()});$('search').addEventListener('input',render);
</script></body></html>'''


def create_dashboard(paths, output):
    runs = read_runs(paths)
    if not runs:
        raise ValueError("No benchmark results supplied")
    # Prevent untrusted JSON from terminating the script tag in the HTML.
    serialized = json.dumps(runs, ensure_ascii=True).replace("<", "\\u003c")
    html = PAGE.replace("__RUN_DATA__", serialized)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", help="pytest-benchmark JSON result paths")
    parser.add_argument("-o", "--output", default=".benchmarks/dashboard.html")
    args = parser.parse_args()
    print(create_dashboard(args.results, args.output))


if __name__ == "__main__":
    main()
