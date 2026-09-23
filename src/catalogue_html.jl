# The catalogue page, embedded in the package so `serve()` has no asset paths
# to resolve and works from a read-only depot. Deliberately self-contained: no
# CDN, no fonts, no charting library, because this is often opened on a cluster
# node with no outbound network. Drawing is done with SVG, which is enough for
# a lattice and a time series and costs nothing.
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
  --bg:#f7f7f5; --panel:#ffffff; --ink:#1b1c1a; --dim:#6a6c66; --line:#dcdcd6;
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
header{padding:14px 20px;border-bottom:1px solid var(--line);display:flex;
       align-items:baseline;gap:16px;flex-wrap:wrap;background:var(--panel)}
header h1{margin:0;font-size:16px;font-weight:650;letter-spacing:-.01em}
header .stat{color:var(--dim);font-size:13px}
header .stat b{color:var(--ink);font-variant-numeric:tabular-nums}
#root{display:grid;grid-template-columns:290px 1fr;gap:0;height:calc(100vh - 53px)}
#side{border-right:1px solid var(--line);overflow:auto;background:var(--panel)}
#main{overflow:auto;padding:18px 20px 40px}
.proj{padding:10px 14px;border-bottom:1px solid var(--line)}
.proj h2{margin:0 0 2px;font-size:13px;font-weight:650}
.proj .meta{color:var(--dim);font-size:12px;font-variant-numeric:tabular-nums}
.lat{display:block;width:100%;text-align:left;border:0;background:none;color:inherit;
     font:inherit;padding:7px 14px 7px 24px;cursor:pointer;border-bottom:1px solid var(--line)}
.lat:hover{background:color-mix(in srgb,var(--accent) 10%,transparent)}
.lat.sel{background:color-mix(in srgb,var(--accent) 18%,transparent);
         box-shadow:inset 3px 0 0 var(--accent)}
.lat .n{font-family:var(--mono);font-size:12px}
.lat .m{color:var(--dim);font-size:11.5px;font-variant-numeric:tabular-nums}
h3{margin:0 0 10px;font-size:14px;font-weight:650}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:6px;
       padding:14px;margin-bottom:16px}
.row{display:flex;gap:16px;flex-wrap:wrap;align-items:flex-start}
.grow{flex:1 1 380px;min-width:0}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{text-align:left;padding:4px 8px;border-bottom:1px solid var(--line);white-space:nowrap}
th{color:var(--dim);font-weight:600;position:sticky;top:0;background:var(--panel)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;font-family:var(--mono)}
tbody tr{cursor:pointer}
tbody tr:hover{background:color-mix(in srgb,var(--accent) 9%,transparent)}
tbody tr.sel{background:color-mix(in srgb,var(--accent) 18%,transparent)}
.scroll{max-height:340px;overflow:auto;border:1px solid var(--line);border-radius:5px}
.tag{display:inline-block;padding:1px 6px;border-radius:99px;font-size:11px;
     border:1px solid var(--line);color:var(--dim);margin-right:4px}
.tag.x{border-color:var(--accent);color:var(--accent)}
.tag.z{border-color:var(--accent2);color:var(--accent2)}
select,button.ctl{font:inherit;font-size:12.5px;padding:3px 8px;border:1px solid var(--line);
                  border-radius:5px;background:var(--panel);color:var(--ink)}
.dim{color:var(--dim)}
.mono{font-family:var(--mono);font-size:12px}
.empty{color:var(--dim);padding:28px 4px}
svg{display:block;max-width:100%}
.legend{display:flex;gap:12px;flex-wrap:wrap;font-size:12px;color:var(--dim);margin-top:6px}
.legend i{display:inline-block;width:16px;height:3px;vertical-align:middle;margin-right:5px}
</style>
</head>
<body>
<header>
  <h1>METTS product-state library</h1>
  <span class="stat"><b id="hFiles">–</b> runs</span>
  <span class="stat"><b id="hSamples">–</b> samples</span>
  <span class="stat mono" id="hRoot"></span>
</header>
<div id="root">
  <nav id="side"></nav>
  <main id="main"><div class="empty">Choose a lattice on the left.</div></main>
</div>
<script>
const $ = s => document.querySelector(s);
const fmt = n => n==null ? "–" : n.toLocaleString();
const g = (n,d=4) => n==null ? "–" :
  (Math.abs(n)>=1e4||(n!==0&&Math.abs(n)<1e-3) ? n.toExponential(2) : +n.toFixed(d));
let SUM=null, CUR=null, RUNS=[], SELROW=null;

async function boot(){
  SUM = await (await fetch("/api/summary")).json();
  $("#hFiles").textContent = fmt(SUM.nfiles);
  $("#hSamples").textContent = fmt(SUM.nsamples);
  $("#hRoot").textContent = SUM.root;
  const side = $("#side");
  for(const p of SUM.projects){
    const d = document.createElement("div");
    d.className = "proj";
    d.innerHTML = `<h2>${p.model} / ${p.project}</h2>
      <div class="meta">${fmt(p.nfiles)} runs · ${fmt(p.nsamples)} samples</div>`;
    side.appendChild(d);
    const lats = Object.values(p.lattices).sort((a,b)=>a.name.localeCompare(b.name));
    for(const l of lats){
      const b = document.createElement("button");
      b.className = "lat";
      b.innerHTML = `<div class="n">${l.name}</div>
        <div class="m">${l.nsites??"?"} sites · ${fmt(l.nfiles)} runs ·
        ${l.ntemperatures} T · ${l.bases.join("/")}</div>`;
      b.onclick = () => { document.querySelectorAll(".lat").forEach(e=>e.classList.remove("sel"));
                          b.classList.add("sel"); openLattice(p, l); };
      side.appendChild(b);
    }
  }
}

async function openLattice(p, l){
  CUR = {p, l};
  $("#main").innerHTML = `<div class="empty">Loading ${l.name}…</div>`;
  const [lat, ens] = await Promise.all([
    fetch(`/api/lattice?project=${encodeURIComponent(p.project)}&lattice=${encodeURIComponent(l.name)}`).then(r=>r.json()),
    fetch(`/api/ensembles?project=${encodeURIComponent(p.project)}&lattice=${encodeURIComponent(l.name)}`).then(r=>r.json())
  ]);
  RUNS = ens.runs;
  $("#main").innerHTML = `
    <div class="row">
      <div class="panel grow"><h3>${l.name}</h3>
        <div id="latplot"></div><div class="legend" id="latlegend"></div>
        <div class="dim" style="margin-top:8px">${lat.nsites} sites ·
          ${lat.bonds?lat.bonds.length:0} bonds · couplings ${(lat.couplings||[]).join(", ")}</div>
      </div>
      <div class="panel grow"><h3>Runs <span class="dim">(${RUNS.length})</span></h3>
        <div class="scroll"><table id="runs"><thead><tr>
          <th class="num">T</th><th>parameters</th><th>sector</th>
          <th class="num">maxdim</th><th class="num">seed</th><th>basis</th>
          <th class="num">samples</th></tr></thead><tbody></tbody></table></div>
        <div class="dim" style="margin-top:8px">Click a run to plot its per-sample observables.</div>
      </div>
    </div>
    <div class="panel" id="seriespanel"><h3>Time series</h3>
      <div class="empty">No run selected.</div></div>`;
  drawLattice(lat);
  fillRuns();
}

function drawLattice(lat){
  const el = $("#latplot");
  if(lat.error || !lat.coordinates){ el.innerHTML = `<div class="empty">${lat.error||"no lattice"}</div>`; return; }
  const pts = lat.coordinates.map(c => [c[0], c.length>1 ? c[1] : 0]);
  const xs = pts.map(p=>p[0]), ys = pts.map(p=>p[1]);
  const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  const W=560, pad=18;
  const sx = (x1>x0) ? (W-2*pad)/(x1-x0) : 1;
  const H = Math.max(120, Math.min(420, (y1-y0)*sx + 2*pad));
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
    // A bond that wraps the periodic direction would cross the whole picture;
    // draw it dashed and short rather than lying about the geometry.
    const dx=Math.abs(pts[i][0]-pts[j][0]), dy=Math.abs(pts[i][1]-pts[j][1]);
    const long = (dx>(x1-x0)*0.5)||(dy>(y1-y0)*0.5);
    svg += `<line x1="${X(pts[i][0])}" y1="${Y(pts[i][1])}" x2="${X(pts[j][0])}" y2="${Y(pts[j][1])}"
             stroke="${c}" stroke-width="${long?1:1.6}" ${long?'stroke-dasharray="3 3" opacity="0.45"':''}/>`;
  }
  pts.forEach((p,i)=>{ svg += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="3.1"
      fill="var(--panel)" stroke="var(--ink)" stroke-width="1.1"><title>site ${i}</title></circle>`; });
  svg += `</svg>`;
  el.innerHTML = svg;
  $("#latlegend").innerHTML = names.map((n,k)=>
    `<span><i style="background:${cols[k%cols.length]}"></i>${n}</span>`).join("")
    + `<span><i style="background:var(--dim);opacity:.5"></i>wraps the boundary</span>`;
}

function fillRuns(){
  const tb = document.querySelector("#runs tbody");
  tb.innerHTML = "";
  for(const r of RUNS){
    const tr = document.createElement("tr");
    const pars = Object.entries(r.parameters).sort().map(([k,v])=>`${k}=${g(v)}`).join(" ");
    const sec  = Object.entries(r.sector).sort().map(([k,v])=>`${k}=${v}`).join(" ");
    const bas  = (r.collapse_bases||[]).map(b=>`<span class="tag ${b.toLowerCase()}">${b}</span>`).join("");
    tr.innerHTML = `<td class="num">${g(r.temperature)}</td><td class="mono">${pars}</td>
      <td class="mono">${sec}</td><td class="num">${r.maxdim??"–"}</td>
      <td class="num">${r.seed??"–"}</td><td>${bas}</td>
      <td class="num">${fmt(r.nsamples)}</td>`;
    tr.onclick = () => { if(SELROW) SELROW.classList.remove("sel");
                         tr.classList.add("sel"); SELROW = tr; openSeries(r); };
    tb.appendChild(tr);
  }
}

async function openSeries(r){
  const p = $("#seriespanel");
  p.innerHTML = `<h3>Time series</h3><div class="empty">Loading…</div>`;
  const d = await (await fetch(`/api/series?path=${encodeURIComponent(r.path)}`)).json();
  if(d.error){ p.innerHTML = `<h3>Time series</h3><div class="empty">${d.error}</div>`; return; }
  const names = Object.keys(d.observables).sort();
  if(!names.length){ p.innerHTML = `<h3>Time series</h3><div class="empty">this run stores no observables</div>`; return; }
  const pref = names.includes("energy") ? "energy" : names[0];
  p.innerHTML = `<h3>Time series <span class="dim mono">${r.path.split("/").pop()}</span></h3>
    <div class="row" style="align-items:center;gap:10px;margin-bottom:10px">
      <select id="obs">${names.map(n=>`<option ${n===pref?"selected":""}>${n}</option>`).join("")}</select>
      <span class="dim">${fmt(d.nsamples)} samples · T=${g(d.temperature)} · ${d.nsites} sites</span>
      <label class="dim"><input type="checkbox" id="cum"> running mean</label>
    </div>
    <div id="plot"></div><div class="dim" id="stats" style="margin-top:6px"></div>`;
  const redraw = () => drawSeries(d.observables[$("#obs").value], $("#obs").value, $("#cum").checked);
  $("#obs").onchange = redraw; $("#cum").onchange = redraw;
  redraw();
}

function drawSeries(v, name, cum){
  const W=1000, H=260, L=64, R=14, T=12, B=30;
  const n = v.length;
  let y = v;
  if(cum){ let s=0; y = v.map((x,i)=>{ s+=x; return s/(i+1); }); }
  const lo=Math.min(...y), hi=Math.max(...y), span=(hi-lo)||1;
  const X = i => L + (n<2?0:(i/(n-1))*(W-L-R));
  const Y = q => T + (1-(q-lo)/span)*(H-T-B);
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" role="img" aria-label="${name}">`;
  for(let k=0;k<=4;k++){
    const q = lo + span*k/4, yy = Y(q);
    svg += `<line x1="${L}" y1="${yy}" x2="${W-R}" y2="${yy}" stroke="var(--line)"/>
            <text x="${L-8}" y="${yy+4}" text-anchor="end" font-size="11"
              fill="var(--dim)" font-family="var(--mono)">${g(q,4)}</text>`;
  }
  // one point per pixel column at most: a 9,000-sample chain does not need
  // 9,000 path segments, and the shape is what matters
  const step = Math.max(1, Math.floor(n/(W-L-R)));
  let dpath = "";
  for(let i=0;i<n;i+=step) dpath += (dpath?"L":"M") + X(i).toFixed(1) + " " + Y(y[i]).toFixed(1);
  svg += `<path d="${dpath}" fill="none" stroke="var(--accent)" stroke-width="1.3"/>`;
  svg += `<text x="${L}" y="${H-8}" font-size="11" fill="var(--dim)">sample 0</text>
          <text x="${W-R}" y="${H-8}" text-anchor="end" font-size="11" fill="var(--dim)">${n-1}</text>
          <text x="${(L+W-R)/2}" y="${H-8}" text-anchor="middle" font-size="11"
            fill="var(--dim)">${name}${cum?" (running mean)":""} vs sample index</text></svg>`;
  $("#plot").innerHTML = svg;
  const mean = v.reduce((a,b)=>a+b,0)/n;
  const sd = Math.sqrt(v.reduce((a,b)=>a+(b-mean)**2,0)/Math.max(1,n-1));
  $("#stats").innerHTML =
    `mean <b class="mono">${g(mean,5)}</b> · sd <b class="mono">${g(sd,4)}</b> ·
     min <b class="mono">${g(Math.min(...v),5)}</b> · max <b class="mono">${g(Math.max(...v),5)}</b>
     <span class="dim">— these runs used nwarm=0, so early samples are burn-in</span>`;
}

boot();
</script>
</body>
</html>
"""
