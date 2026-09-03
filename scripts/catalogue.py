#!/usr/bin/env python
"""Harvest the Wayne State DTE Aerial collection catalogue.

The collection is 13,500-odd scanned negatives across several counties and
several missions. Each record's title carries the county, the year and a section
code (e.g. "Wayne County, 1956, Section ga-15-22"), and the pointer is what the
IIIF endpoint wants. This pulls the lot once so we can ask what exists before
downloading anything.

Usage:  ./.venv/bin/python scripts/catalogue.py [--out data/dte_catalogue.json]
"""
import sys, os, json, time, re, urllib.request

BASE = "https://wayne.contentdm.oclc.org/digital/bl/dmwebservices/index.php"
UA = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                    'AppleWebKit/537.36 Chrome/126.0 Safari/537.36'}
PAGE = 1000
TITLE = re.compile(r'DTE Aerial Photograph,\s*([^,]+?),\s*(\d{4}),\s*Section\s*(\S+)', re.I)


def get(url, tries=4):
    for i in range(tries):
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=120).read())
        except Exception as exc:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def main():
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else 'data/dte_catalogue.json'
    q = f"{BASE}?q=dmQuery/dte-aerial/0/title!date/nosort/{PAGE}/%d/1/0/0/0/0/json"
    first = get(q % 0)
    total = int(first['pager']['total'])
    print(f"collection holds {total} records", flush=True)
    recs, start = [], 0
    while start < total:
        d = first if start == 0 else get(q % start)
        rows = d.get('records') or []
        if not rows:
            break
        recs.extend(rows)
        start += len(rows)
        print(f"  {start}/{total}", flush=True)
        time.sleep(0.4)
    items = []
    for r in recs:
        t = r.get('title') or ''
        m = TITLE.search(t)
        items.append(dict(pointer=r.get('pointer'), title=t,
                          county=(m.group(1).strip() if m else None),
                          year=(m.group(2) if m else (r.get('date') or None)),
                          section=(m.group(3) if m else None)))
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    json.dump(dict(total=total, items=items), open(out, 'w'))
    print(f"wrote {out}: {len(items)} records", flush=True)


if __name__ == '__main__':
    main()
