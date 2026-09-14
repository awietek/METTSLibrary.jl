#!/bin/bash
# Splice data.json into the page template -> catalogue.html
set -eu
W=${METTS_SCAN_DIR:-$PWD/scan}
HERE=$(cd "$(dirname "$0")" && pwd)

[ -s "$W/data.json" ] || { echo "data.json missing - run catalogue_data.jl"; exit 1; }

python3 - "$HERE/catalogue_page.html" "$W/data.json" "$W/catalogue.html" <<'PY'
import sys, pathlib, json
tpl, data, out = (pathlib.Path(p) for p in sys.argv[1:4])
t = tpl.read_text(); d = data.read_text()
assert '__DATA__' in t, 'template has no __DATA__ placeholder'
out.write_text(t.replace('__DATA__', d))
j = json.loads(d)
print(f"{len(j['fams'])} families, {sum(len(f[19]) for f in j['fams'])} ensembles")
print(f"wrote {out} ({out.stat().st_size} bytes)")
PY
