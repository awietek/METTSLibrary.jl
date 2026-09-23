# The catalogue page, embedded in the package so `serve()` has no asset paths
# to resolve and works from a read-only depot. Deliberately self-contained: no
# CDN, no fonts, no charting library, because this is often opened on a cluster
# node with no outbound network. Drawing is done with SVG, which is enough for
# a lattice and a time series and costs nothing.
#
# Navigation mirrors the library layout --
#   <model>/<project>/<lattice>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5
# -- one level at a time, behind a breadcrumb. A sidebar listing every project
# and its lattices does not survive the archive growing to dozens of projects.
#
# raw"" so the JavaScript's $ and \ are literal.

const CATALOGUE_HTML = raw"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>METTS library</title>
<style>
:root{
  --bg:#f7f7f5; --panel:#fff; --ink:#1b1c1a; --dim:#6a6c66; --line:#dcdcd6;
  --accent:#3a6ea5; --accent2:#a8562a; --good:#2f7d52; --warn:#b06a1c;
  --mono:ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --bg:#17181a; --panel:#1f2023; --ink:#e8e8e4; --dim:#9a9c96; --line:#33353a;
    --accent:#7fb0e0; --accent2:#e0956a; --good:#6fcb96; --warn:#e0b070;
  }
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:14px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
header{padding:12px 22px;border-bottom:1px solid var(--line);background:var(--panel);
       display:flex;align-items:baseline;gap:18px;flex-wrap:wrap}
header h1{margin:0;font-size:15px;font-weight:650;letter-spacing:-.01em;cursor:pointer}
header h1:hover{color:var(--accent)}
header .stat{color:var(--dim);font-size:12.5px}
header .stat b{color:var(--ink);font-variant-numeric:tabular-nums}
header .root{margin-left:auto;color:var(--dim);font-family:var(--mono);font-size:11.5px}
#crumbs{padding:9px 22px;border-bottom:1px solid var(--line);font-size:13px;
        display:flex;gap:7px;flex-wrap:wrap;align-items:center;background:var(--panel)}
#crumbs a{color:var(--accent);text-decoration:none;cursor:pointer;font-family:var(--mono)}
#crumbs a:hover{text-decoration:underline}
#crumbs .sep{color:var(--dim)}
#crumbs .lvl{margin-left:8px;color:var(--dim);font-size:12px}
main{padding:20px 22px 60px;max-width:1500px}
h2.sec{margin:26px 0 12px;font-size:13px;font-weight:650;text-transform:uppercase;
       letter-spacing:.06em;color:var(--dim);border-bottom:1px solid var(--line);
       padding-bottom:6px}
h2.sec:first-child{margin-top:0}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(270px,1fr));gap:12px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:13px 15px;
      cursor:pointer;text-align:left;font:inherit;color:inherit;transition:border-color .1s}
.card:hover{border-color:var(--accent)}
.card .t{font-weight:600;font-size:13.5px;font-family:var(--mono);word-break:break-all}
.card .s{color:var(--dim);font-size:12px;margin-top:5px;font-variant-numeric:tabular-nums}
.card .c{margin-top:7px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:7px;
       padding:15px;margin-bottom:16px}
.row{display:flex;gap:16px;flex-wrap:wrap;align-items:flex-start}
.grow{flex:1 1 420px;min-width:0}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{text-align:left;padding:5px 9px;border-bottom:1px solid var(--line);white-space:nowrap}
th{color:var(--dim);font-weight:600;position:sticky;top:0;background:var(--panel)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono)}
tbody tr{cursor:pointer}
tbody tr:hover{background:color-mix(in srgb,var(--accent) 9%,transparent)}
tbody tr.sel{background:color-mix(in srgb,var(--accent) 18%,transparent)}
.scroll{max-height:400px;overflow:auto;border:1px solid var(--line);border-radius:6px}
.tag{display:inline-block;padding:1px 7px;border-radius:99px;font-size:11px;
     border:1px solid var(--line);color:var(--dim);margin-right:4px}
.tag.x{border-color:var(--accent);color:var(--accent)}
.tag.z{border-color:var(--accent2);color:var(--accent2)}
select,button.ctl{font:inherit;font-size:12.5px;padding:3px 8px;border:1px solid var(--line);
                  border-radius:5px;background:var(--panel);color:var(--ink)}
.dim{color:var(--dim)}
.mono{font-family:var(--mono);font-size:12px}
.empty{color:var(--dim);padding:26px 4px}
svg{display:block;max-width:100%}
.legend{display:flex;gap:13px;flex-wrap:wrap;font-size:12px;color:var(--dim);margin-top:7px}
.legend i{display:inline-block;width:16px;height:3px;vertical-align:middle;margin-right:5px}
.swatch{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle}
.md{font-size:13.5px;max-width:78ch}
.md h3,.md h4,.md h5{margin:18px 0 7px;font-size:14px;font-weight:650}
.md h3:first-child{margin-top:0}
.md p{margin:0 0 10px}
.md ul,.md ol{margin:0 0 10px;padding-left:22px}
.md li{margin:2px 0}
.md code{font-family:var(--mono);font-size:12px;background:var(--bg);
         border:1px solid var(--line);border-radius:3px;padding:0 4px}
.md table{margin:0 0 12px;font-size:12.5px;width:auto}
.md th,.md td{white-space:normal}
.md hr{border:0;border-top:1px solid var(--line);margin:16px 0}
.md a{color:var(--accent)}
tbody tr:hover{background:color-mix(in srgb,var(--accent) 7%,transparent)}
td input[type=checkbox],th input[type=checkbox]{cursor:pointer;margin:0}
</style>
</head>
<body>
<header>
  <h1 onclick="go('')">METTS product-state library</h1>
  <span class="stat"><b id="hFiles">–</b> runs</span>
  <span class="stat"><b id="hSamples">–</b> samples</span>
  <span class="root" id="hRoot"></span>
</header>
<div id="crumbs"></div>
<main id="main"><div class="empty">Loading…</div></main>
<script>
const $ = s => document.querySelector(s);
const fmt = n => n==null ? "–" : n.toLocaleString();
const g = (n,d=4) => n==null ? "–" :
  (Math.abs(n)>=1e4||(n!==0&&Math.abs(n)<1e-3) ? n.toExponential(2) : +n.toFixed(d));
const enc = encodeURIComponent;

// the address bar follows the library path, so a level can be linked and the
// browser's back button does the obvious thing
addEventListener("popstate", () => render(location.hash.slice(1)));
function go(p){ history.pushState(null,"","#"+p); render(p); }

async function render(path){
  const d = await (await fetch(`/api/browse?path=${enc(path||"")}`)).json();
  $("#hFiles").textContent   = fmt(d.libfiles);
  $("#hSamples").textContent = fmt(d.libsamples);
  $("#hRoot").textContent    = d.root;
  crumbs(d);
  if(d.error){ $("#main").innerHTML = `<div class="empty">${d.error}</div>`; return; }
  if(d.childLevel === "run")   return renderRuns(d);
  if(d.childLevel === "model") return renderHome(d);
  return renderLevel(d);
}

function crumbs(d){
  const el = $("#crumbs");
  let h = `<a onclick="go('')">library</a>`;
  for(const c of d.crumbs) h += `<span class="sep">/</span><a onclick="go('${c.path}')">${c.name}</a>`;
  h += `<span class="lvl">${fmt(d.nfiles)} runs · ${fmt(d.nsamples)} samples`;
  if(d.childLevel !== "run") h += ` · choose a ${d.childLevel}`;
  h += `</span>`;
  el.innerHTML = h;
}

// Landing page: models as sections, their projects as cards. Everything the
// archive gains later lands under an existing model heading.
function renderHome(d){
  let h = "";
  for(const m of d.children){
    h += `<h2 class="sec">${m.name} <span class="dim" style="text-transform:none;
           font-weight:400">— ${m.nprojects} project${m.nprojects===1?"":"s"},
           ${fmt(m.nfiles)} runs, ${fmt(m.nsamples)} samples</span></h2>
          <div class="cards" id="m-${m.name}"></div>`;
  }
  $("#main").innerHTML = h;
  for(const m of d.children) fillProjects(m);
}

async function fillProjects(m){
  const d = await (await fetch(`/api/browse?path=${enc(m.path)}`)).json();
  const el = document.getElementById(`m-${m.name}`);
  if(!el) return;
  el.innerHTML = d.children.map(p => `
    <button class="card" onclick="go('${p.path}')">
      <div class="t">${p.name}</div>
      <div class="s">${fmt(p.nfiles)} runs · ${fmt(p.nsamples)} samples</div>
      <div class="s">${p.nlattices} lattice${p.nlattices===1?"":"s"} ·
           ${p.ntemperatures} temperature${p.ntemperatures===1?"":"s"}</div>
    </button>`).join("");
}

// Intermediate levels: lattices, parameter sets, sectors, temperatures.
function renderLevel(d){
  const cards = d.children.map(c => {
    let sub = "";
    if(d.childLevel === "lattice"){
      sub = `<div class="s">${c.nsites ?? "?"} sites · ${c.ntemperatures} T ·
             ${(c.bases||[]).join("/")}</div>
             <div class="s dim">${(c.couplings||[]).join(", ")}</div>`;
    } else if(d.childLevel === "parameters"){
      sub = `<div class="s">${Object.entries(c.parameters||{}).sort()
             .map(([k,v])=>`${k}=${g(v)}`).join("  ")}</div>
             <div class="s dim">${c.nsectors} sector${c.nsectors===1?"":"s"}</div>`;
    } else if(d.childLevel === "sector"){
      sub = `<div class="s">${Object.entries(c.sector||{}).sort()
             .map(([k,v])=>`${k} = ${v}`).join("  ")}</div>
             <div class="s dim">${c.ntemperatures} temperature${c.ntemperatures===1?"":"s"}</div>`;
    }
    return `<button class="card" onclick="go('${c.path}')">
        <div class="t">${c.name}</div>
        <div class="s">${fmt(c.nfiles)} run${c.nfiles===1?"":"s"} ·
             ${fmt(c.nsamples)} samples</div>${sub}</button>`;
  }).join("");

  // A sweep is dozens of temperatures. Cards would be a wall; a plain list
  // reads top to bottom and sorts naturally by T.
  const body = (d.childLevel === "temperature")
    ? `<div class="panel"><h2 class="sec" style="margin-top:0">temperatures
         <span class="dim" style="text-transform:none;font-weight:400">
         (${d.children.length})</span></h2>
       <div class="scroll" style="max-height:560px"><table><thead><tr>
         <th class="num">T</th><th class="num">β</th><th class="num">runs</th>
         <th class="num">samples</th><th>directory</th></tr></thead><tbody>
       ${d.children.map(c=>`<tr onclick="go('${c.path}')">
           <td class="num">${g(c.temperature)}</td>
           <td class="num">${g(c.beta)}</td>
           <td class="num">${fmt(c.nfiles)}</td>
           <td class="num">${fmt(c.nsamples)}</td>
           <td class="mono dim">${c.name}</td></tr>`).join("")}
       </tbody></table></div></div>`
    : `<div class="cards">${cards}</div>`;

  const latBox = d.lattice ? `<div class="panel" id="latpanel"><h2 class="sec"
        style="margin-top:0">lattice</h2><div id="latplot"></div>
        <div class="legend" id="latlegend"></div>
        <div class="dim mono" id="latmeta" style="margin-top:8px"></div></div>` : "";
  // The project page is the one place to put the project's own README: it
  // says who computed the runs and what is known about their quality, which is
  // what a reader needs before using any of it.
  const readmeBox = (d.childLevel === "lattice")
    ? `<div class="panel" id="readme"><h2 class="sec" style="margin-top:0">
         project readme</h2><div class="empty">Loading…</div></div>` : "";
  $("#main").innerHTML = `<div class="row">
      <div class="grow">${body}</div>
      ${d.lattice ? `<div class="grow" style="flex:0 1 770px">${latBox}</div>` : ""}
    </div>${readmeBox}`;
  if(d.lattice) loadLattice(d.lattice);
  if(d.childLevel === "lattice") loadReadme(d.path);
}

async function loadReadme(path){
  const el = $("#readme"); if(!el) return;
  const d = await (await fetch(`/api/readme?path=${enc(path)}`)).json();
  const head = `<h2 class="sec" style="margin-top:0">project readme
      <span class="dim mono" style="text-transform:none;font-weight:400">${path}/README.md</span></h2>`;
  el.innerHTML = head + (d.markdown
    ? `<div class="md">${md(d.markdown)}</div>`
    : `<div class="empty">${d.error || "this project has no README.md"}</div>`);
}

// A small markdown renderer: headings, tables, lists, rules, code, emphasis.
// The project READMEs use exactly that much, and pulling in a markdown library
// would break the "works on a node with no network" property.
function esc(s){ return s.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;"); }
function inline(s){
  return esc(s)
    .replace(/`([^`]+)`/g, '<code>$1</code>')
    .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
    .replace(/(^|[^*])\*([^*]+)\*/g, '$1<em>$2</em>')
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, '<a href="$2" rel="noreferrer">$1</a>');
}
function md(src){
  const out = [];
  const lines = src.split("\n");
  let i = 0, list = null;
  const closeList = () => { if(list){ out.push(`</${list}>`); list = null; } };
  while(i < lines.length){
    const L = lines[i];
    if(/^\s*$/.test(L)){ closeList(); i++; continue; }
    if(/^---+\s*$/.test(L)){ closeList(); out.push("<hr>"); i++; continue; }
    let m;
    if((m = L.match(/^(#{1,6})\s+(.*)$/))){
      closeList(); const n = m[1].length;
      out.push(`<h${Math.min(n+2,6)}>${inline(m[2])}</h${Math.min(n+2,6)}>`); i++; continue;
    }
    // a table: header row, separator, then body until a blank line
    if(/^\s*\|/.test(L) && i+1 < lines.length && /^\s*\|[\s:|-]+\|?\s*$/.test(lines[i+1])){
      closeList();
      const cells = r => r.trim().replace(/^\||\|$/g,"").split("|").map(c=>c.trim());
      let h = `<table><thead><tr>${cells(L).map(c=>`<th>${inline(c)}</th>`).join("")}</tr></thead><tbody>`;
      i += 2;
      while(i < lines.length && /^\s*\|/.test(lines[i])){
        h += `<tr>${cells(lines[i]).map(c=>`<td>${inline(c)}</td>`).join("")}</tr>`; i++;
      }
      out.push(h + "</tbody></table>"); continue;
    }
    if((m = L.match(/^\s*[-*]\s+(.*)$/))){
      if(list !== "ul"){ closeList(); out.push("<ul>"); list = "ul"; }
      out.push(`<li>${inline(m[1])}</li>`); i++; continue;
    }
    if((m = L.match(/^\s*\d+\.\s+(.*)$/))){
      if(list !== "ol"){ closeList(); out.push("<ol>"); list = "ol"; }
      out.push(`<li>${inline(m[1])}</li>`); i++; continue;
    }
    // paragraph: absorb following non-blank, non-structural lines
    closeList();
    let p = [L]; i++;
    while(i < lines.length && !/^\s*$/.test(lines[i]) &&
          !/^(#{1,6}\s|---+\s*$|\s*\|)/.test(lines[i]) &&
          !/^\s*([-*]|\d+\.)\s/.test(lines[i])){ p.push(lines[i]); i++; }
    out.push(`<p>${inline(p.join(" "))}</p>`);
  }
  closeList();
  return out.join("\n");
}

// Deepest level: the runs in one temperature directory. These are the same
// ensemble differing only in algorithm settings — usually just the seed — so
// the useful view is all of them on one axis, not one at a time. Every run is
// checked by default; unchecking is how you isolate one.
const PALETTE = ["#3a6ea5","#a8562a","#2f7d52","#8a5fbf","#b8902a","#3f9199",
                 "#c2566f","#6b8f3a","#8a6a4a","#5f7fbf","#a03b3b","#4a9a6a"];
let RUNS = [], CACHE = new Map();

function renderRuns(d){
  RUNS = d.children; CACHE = new Map();
  const c0 = RUNS[0] || {};
  const rows = RUNS.map((r,i) => {
    const bas = (r.collapse_bases||[]).map(b=>`<span class="tag ${b.toLowerCase()}">${b}</span>`).join("");
    return `<tr>
      <td><input type="checkbox" class="pick" data-i="${i}" checked></td>
      <td><span class="swatch" style="background:${PALETTE[i%PALETTE.length]}"></span></td>
      <td class="mono">${r.name}</td>
      <td class="num">${r.maxdim ?? "–"}</td><td class="num">${r.seed ?? "–"}</td>
      <td class="num">${g(r.tau)}</td><td class="num">${g(r.cutoff)}</td>
      <td>${bas}</td><td class="num">${fmt(r.nsamples)}</td>
      <td class="num">${(r.bytes/1024).toFixed(0)} kB</td></tr>`;
  }).join("");
  const names = [...new Set(RUNS.flatMap(r => r.observables||[]))].sort();
  const pref = names.includes("energy") ? "energy" : names[0];
  $("#main").innerHTML = `<div class="row">
    <div class="panel grow"><h2 class="sec" style="margin-top:0">runs
        <span class="dim" style="text-transform:none;font-weight:400">(${RUNS.length})</span></h2>
      <div class="dim mono" style="margin-bottom:9px">
        T = ${g(c0.temperature)} ·
        ${Object.entries(c0.parameters||{}).sort().map(([k,v])=>`${k}=${g(v)}`).join(" ")} ·
        ${Object.entries(c0.sector||{}).sort().map(([k,v])=>`${k}=${v}`).join(" ")}</div>
      <div class="scroll"><table><thead><tr>
        <th><input type="checkbox" id="pickall" checked title="all / none"></th>
        <th></th><th>tag</th><th class="num">maxdim</th><th class="num">seed</th>
        <th class="num">tau</th><th class="num">cutoff</th><th>basis</th>
        <th class="num">samples</th><th class="num">size</th></tr></thead>
        <tbody>${rows}</tbody></table></div>
    </div>
    <div class="grow" style="flex:0 1 770px"><div class="panel" id="latpanel">
      <h2 class="sec" style="margin-top:0">lattice</h2><div id="latplot"></div>
      <div class="legend" id="latlegend"></div>
      <div class="dim mono" id="latmeta" style="margin-top:8px"></div></div></div>
  </div>
  <div class="panel" id="seriespanel">
    <h2 class="sec" style="margin-top:0">time series</h2>
    <div class="row" style="align-items:center;gap:12px;margin-bottom:10px">
      <select id="obs">${names.map(n=>`<option ${n===pref?"selected":""}>${n}</option>`).join("")}</select>
      <label class="dim"><input type="checkbox" id="cum"> running mean</label>
      <span class="dim" id="plotnote"></span>
    </div>
    <div id="plot"><div class="empty">Loading…</div></div>
    <div id="runlegend" class="legend"></div></div>`;
  $("#pickall").onchange = e => {
    document.querySelectorAll(".pick").forEach(c => c.checked = e.target.checked);
    refreshSeries();
  };
  document.querySelectorAll(".pick").forEach(c => c.onchange = refreshSeries);
  $("#obs").onchange = refreshSeries;
  $("#cum").onchange = () => draw();
  if(d.lattice) loadLattice(d.lattice);
  refreshSeries();
}

// Fetch only the observable being plotted, and only for runs not already held.
// A temperature directory can hold 120 runs; pulling all twelve observables
// for each would be tens of megabytes to draw one curve.
async function refreshSeries(){
  const obs = $("#obs").value;
  const picked = [...document.querySelectorAll(".pick")]
    .filter(c => c.checked).map(c => RUNS[+c.dataset.i]);
  $("#plotnote").textContent = `${picked.length} of ${RUNS.length} runs`;
  if(!picked.length){ $("#plot").innerHTML = `<div class="empty">No run selected.</div>`;
                      $("#runlegend").innerHTML = ""; return; }
  const todo = picked.filter(r => !CACHE.has(r.runpath + "|" + obs));
  if(todo.length){
    $("#plot").innerHTML = `<div class="empty">Loading ${todo.length} run${todo.length===1?"":"s"}…</div>`;
    // a handful at a time: 120 simultaneous HDF5 reads helps nobody
    const queue = todo.slice();
    const worker = async () => {
      while(queue.length){
        const r = queue.shift();
        const d = await (await fetch(
          `/api/series?path=${enc(r.runpath)}&obs=${enc(obs)}`)).json();
        CACHE.set(r.runpath + "|" + obs, d.observables ? d.observables[obs] : null);
      }
    };
    await Promise.all(Array.from({length: Math.min(6, queue.length)}, worker));
  }
  draw();
}

function draw(){
  const obs = $("#obs").value, cum = $("#cum").checked;
  const picked = [...document.querySelectorAll(".pick")]
    .filter(c => c.checked).map(c => ({ r: RUNS[+c.dataset.i], i: +c.dataset.i }));
  const series = picked.map(({r,i}) => ({
      name: r.name, seed: r.seed, maxdim: r.maxdim, i,
      v: CACHE.get(r.runpath + "|" + obs) }))
    .filter(s => s.v && s.v.length);
  if(!series.length){ $("#plot").innerHTML =
      `<div class="empty">no run in this selection stores '${obs}'</div>`;
      $("#runlegend").innerHTML = ""; return; }
  drawSeries(series, obs, cum);
}

function drawSeries(series, name, cum){
  const W=1100, H=300, L=70, R=16, T=12, B=30;
  const ys = series.map(s => {
    if(!cum) return s.v;
    let acc=0; return s.v.map((x,i)=>{ acc+=x; return acc/(i+1); });
  });
  const nmax = Math.max(...series.map(s => s.v.length));
  const lo = Math.min(...ys.map(a=>Math.min(...a)));
  const hi = Math.max(...ys.map(a=>Math.max(...a)));
  const span = (hi-lo)||1;
  const X = i => L + (nmax<2?0:(i/(nmax-1))*(W-L-R));
  const Y = q => T + (1-(q-lo)/span)*(H-T-B);
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" role="img"
                  aria-label="${name} for ${series.length} runs">`;
  for(let k=0;k<=4;k++){
    const q = lo + span*k/4, yy = Y(q);
    svg += `<line x1="${L}" y1="${yy}" x2="${W-R}" y2="${yy}" stroke="var(--line)"/>
            <text x="${L-8}" y="${yy+4}" text-anchor="end" font-size="11"
              fill="var(--dim)" font-family="var(--mono)">${g(q,4)}</text>`;
  }
  series.forEach((s,k) => {
    const y = ys[k], n = y.length;
    // at most one point per pixel column: a 9,000-sample chain does not need
    // 9,000 path segments and the shape is what matters
    const step = Math.max(1, Math.floor(n/(W-L-R)));
    let dp = "";
    for(let i=0;i<n;i+=step) dp += (dp?"L":"M") + X(i).toFixed(1) + " " + Y(y[i]).toFixed(1);
    svg += `<path d="${dp}" fill="none" stroke="${PALETTE[s.i%PALETTE.length]}"
             stroke-width="1.2" opacity="${series.length>20?0.55:0.9}"/>`;
  });
  svg += `<text x="${L}" y="${H-8}" font-size="11" fill="var(--dim)">sample 0</text>
          <text x="${W-R}" y="${H-8}" text-anchor="end" font-size="11" fill="var(--dim)">${nmax-1}</text>
          <text x="${(L+W-R)/2}" y="${H-8}" text-anchor="middle" font-size="11"
            fill="var(--dim)">${name}${cum?" (running mean)":""} vs sample index —
            these runs used nwarm=0, so early samples are burn-in</text></svg>`;
  $("#plot").innerHTML = svg;
  // a legend of 120 entries is noise; past a dozen just give the spread
  if(series.length <= 14){
    $("#runlegend").innerHTML = series.map((s,k) => {
      const m = s.v.reduce((a,b)=>a+b,0)/s.v.length;
      return `<span><i style="background:${PALETTE[s.i%PALETTE.length]}"></i>
        seed ${s.seed ?? "?"}${s.maxdim?` · D=${s.maxdim}`:""}
        <span class="mono">${g(m,5)}</span></span>`;
    }).join("");
  } else {
    const means = series.map(s => s.v.reduce((a,b)=>a+b,0)/s.v.length);
    const mm = means.reduce((a,b)=>a+b,0)/means.length;
    const sd = Math.sqrt(means.reduce((a,b)=>a+(b-mm)**2,0)/Math.max(1,means.length-1));
    $("#runlegend").innerHTML =
      `<span>${series.length} runs · mean of per-run means
       <span class="mono">${g(mm,5)}</span> · spread between chains
       <span class="mono">${g(sd,4)}</span></span>`;
  }
}

async function loadLattice(p){
  const lat = await (await fetch(`/api/lattice?path=${enc(p)}`)).json();
  const el = $("#latplot"); if(!el) return;
  if(lat.error || !lat.coordinates){ el.innerHTML = `<div class="empty">${lat.error||"no lattice"}</div>`; return; }
  const pts = lat.coordinates.map(c => [c[0], c.length>1 ? c[1] : 0]);
  const xs = pts.map(p=>p[0]), ys = pts.map(p=>p[1]);
  const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  // padding leaves room for the site labels, which sit outside the markers
  const W=730, pad=26;
  const sx = (x1>x0) ? (W-2*pad)/(x1-x0) : 1;
  const H = Math.max(130, Math.min(580, (y1-y0)*sx + 2*pad));
  const sy = (y1>y0) ? (H-2*pad)/(y1-y0) : 1;
  const s = Math.min(sx, sy);
  const X = x => pad + (x-x0)*s, Y = y => H - pad - (y-y0)*s;
  const cols = ["var(--accent)","var(--accent2)","var(--good)","var(--warn)","#8a7fd0"];
  const names = [...new Set((lat.bonds||[]).map(b=>b.coupling))];
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="lattice">`;
  for(const b of (lat.bonds||[])){
    const c = cols[names.indexOf(b.coupling)%cols.length];
    const [i,j] = b.sites;
    if(i>=pts.length||j>=pts.length) continue;
    // a bond that wraps the periodic direction would cross the whole picture;
    // draw it dashed rather than lying about the geometry
    const dx=Math.abs(pts[i][0]-pts[j][0]), dy=Math.abs(pts[i][1]-pts[j][1]);
    const long = (dx>(x1-x0)*0.5)||(dy>(y1-y0)*0.5);
    svg += `<line x1="${X(pts[i][0])}" y1="${Y(pts[i][1])}" x2="${X(pts[j][0])}" y2="${Y(pts[j][1])}"
             stroke="${c}" stroke-width="${long?1:1.6}" ${long?'stroke-dasharray="3 3" opacity="0.45"':''}/>`;
  }
  // Site indices annotated beside each marker, not inside it: the stored states
  // are ordered by site index, so reading a configuration means knowing which
  // site is which — but putting the number in the circle forces the marker far
  // larger than the lattice wants.
  const near = pts.length > 1 ? Math.min(...pts.slice(1).map((p,k)=>
      Math.hypot((p[0]-pts[k][0])*s, (p[1]-pts[k][1])*s)).filter(v=>v>0)) : 20;
  const fs = Math.max(7, Math.min(11, near*0.34));
  pts.forEach((p,i)=>{
    svg += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="3.2"
        fill="var(--panel)" stroke="var(--ink)" stroke-width="1.1"/>
      <text x="${(X(p[0])+4.5).toFixed(1)}" y="${(Y(p[1])-4.5).toFixed(1)}"
        font-size="${fs.toFixed(1)}" font-family="var(--mono)" fill="var(--dim)"
        style="pointer-events:none">${i}</text>`;
  });
  svg += `</svg>`;
  el.innerHTML = svg;
  $("#latlegend").innerHTML = names.map((n,k)=>
    `<span><i style="background:${cols[k%cols.length]}"></i>${n}</span>`).join("")
    + `<span><i style="background:var(--dim);opacity:.5"></i>wraps the boundary</span>`;
  $("#latmeta").textContent =
    `${lat.name} — ${lat.nsites} sites, ${lat.bonds.length} bonds, couplings ${(lat.couplings||[]).join(", ")}`;
}

render(location.hash.slice(1));
</script>
</body>
</html>
"""
