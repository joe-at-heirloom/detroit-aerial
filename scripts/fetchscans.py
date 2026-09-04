#!/usr/bin/env python
"""Download DTE negatives from Wayne State's IIIF endpoint into a stable store.

The scans were once kept in a session scratchpad under /tmp and were lost when the
machine cleared it, taking every COLMAP solve with them. They live in scans/ now,
beside the repo and gitignored, so a reboot costs nothing.

Selection is by catalogue (data/dte_catalogue.json), so you ask for what you want
in the collection's own terms:

  ./.venv/bin/python scripts/fetchscans.py --year 1949 --roll ha-17
  ./.venv/bin/python scripts/fetchscans.py --year 1961 --county "Wayne County"
  ./.venv/bin/python scripts/fetchscans.py --pointers 452,607,1408
"""
import sys, os, json, time, threading, queue, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                    'AppleWebKit/537.36 Chrome/126.0 Safari/537.36'}
# Ask for native resolution. Requesting a fixed width fails with HTTP 501 on any
# negative scanned narrower than it -- the endpoint is IIIF level1 and will not
# upscale -- and the widths here run 5352 to 5400. A fixed 5366 is how the first
# scan store came to be missing frames. Heights are the stable dimension (~5100),
# and colmap_block normalises on height, so native costs nothing.
WIDTH = 'full'
MINBYTES = 50_000


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def select(items):
    y, roll, county = arg('--year'), arg('--roll'), arg('--county')
    if '--pointers' in sys.argv:
        want = {int(x) for x in arg('--pointers').replace(' ', '').split(',') if x}
        return [i for i in items if i['pointer'] in want]
    out = items
    if y: out = [i for i in out if i['year'] == y]
    if county: out = [i for i in out if i['county'] == county]
    if roll: out = [i for i in out if i['section'] and i['section'].rsplit('-', 1)[0] == roll]
    return out


ERRORS = {}


def fetch(rec, out, width, tries=5):
    p = os.path.join(out, f"{rec}.jpg")
    if os.path.exists(p) and os.path.getsize(p) > MINBYTES:
        return 'have'
    size = 'full' if str(width) == 'full' else f"{width},"
    url = f"https://wayne.contentdm.oclc.org/digital/iiif/dte-aerial/{rec}/full/{size}/0/default.jpg"
    for i in range(tries):
        try:
            d = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=240).read()
            if len(d) < MINBYTES:
                return 'small'
            tmp = p + '.part'
            open(tmp, 'wb').write(d); os.replace(tmp, p)
            return 'got'
        except Exception as exc:
            ERRORS[f"{type(exc).__name__}: {exc}"[:90]] = ERRORS.get(f"{type(exc).__name__}: {exc}"[:90], 0) + 1
            if i == tries - 1:
                return 'fail'
            time.sleep(5 * (i + 1))      # a library server, not a CDN: back off properly


def main():
    out = os.path.join(ROOT, arg('--out', 'scans/fullres'))
    width = arg('--width', WIDTH)
    workers = int(arg('--workers', 2))
    items = json.load(open(os.path.join(ROOT, 'data', 'dte_catalogue.json')))['items']
    sel = select(items)
    if not sel:
        raise SystemExit('nothing selected; check --year / --roll / --county / --pointers')
    os.makedirs(out, exist_ok=True)
    print(f"{len(sel)} frames -> {out} at width {width}", flush=True)
    q = queue.Queue()
    for i in sel: q.put(i['pointer'])
    tally, lock = {'got': 0, 'have': 0, 'small': 0, 'fail': 0}, threading.Lock()

    def work():
        while True:
            try: rec = q.get_nowait()
            except queue.Empty: return
            r = fetch(rec, out, width)
            with lock:
                tally[r] += 1
                n = sum(tally.values())
                if n % 10 == 0 or n == len(sel):
                    print(f"  {n}/{len(sel)}  {tally}", flush=True)
            time.sleep(1.0)

    ts = [threading.Thread(target=work, daemon=True) for _ in range(workers)]
    for t in ts: t.start()
    for t in ts: t.join()
    print(f"done: {tally}", flush=True)
    for e, n in sorted(ERRORS.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d} x {e}", flush=True)
    have = [f for f in os.listdir(out) if f.endswith('.jpg')]
    print(f"store now holds {len(have)} scans, "
          f"{sum(os.path.getsize(os.path.join(out,f)) for f in have)/2**30:.2f} GiB", flush=True)


if __name__ == '__main__':
    main()
