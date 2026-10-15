#!/usr/bin/env python3
"""Resolve search leads through DOI-indexed metadata; verify downloaded XML identity."""
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

base = Path('resources/literature/darkmatter-2026-10-13')
dest = base / 'verified'
dest.mkdir(exist_ok=True)

def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'P4-literature-audit/1.0'})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.read()

def one(doi):
    key = hashlib.sha256(doi.lower().encode()).hexdigest()[:16]
    info = {'doi': doi, 'status': 'unresolved', 'path': '', 'title': '', 'year': '', 'error': ''}
    cache = dest / (key + '.json')
    if cache.exists():
        old = json.loads(cache.read_text())
        if old.get('status') == 'fulltext_doi_verified':
            return old
    try:
        url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': 'DOI:"' + doi + '"', 'format': 'json', 'resultType': 'core'})
        data = json.loads(fetch(url))
        hits = [x for x in data['resultList']['result'] if x.get('doi', '').lower() == doi.lower()]
        if hits:
            hit = hits[0]
            info.update(title=hit.get('title', ''), year=hit.get('pubYear', ''), metadata=hit, status='doi_metadata_verified')
            if hit.get('pmcid'):
                src = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + hit['pmcid'] + '/fullTextXML'
                try:
                    blob = fetch(src)
                    tree = ET.fromstring(blob)
                    ids = [x.text.lower() for x in tree.findall('.//article-meta/article-id') if x.get('pub-id-type') == 'doi' and x.text]
                    if doi.lower() not in ids:
                        raise ValueError('Downloaded full text has different DOI: ' + repr(ids))
                    p = dest / (key + '.xml'); p.write_bytes(blob)
                    info.update(status='fulltext_doi_verified', path=str(p), source_url=src, sha256=hashlib.sha256(blob).hexdigest())
                except Exception as e:
                    info['error'] = 'Fulltext: ' + str(e)
        else:
            url = 'https://api.crossref.org/works/' + urllib.parse.quote(doi, safe='')
            hit = json.loads(fetch(url))['message']
            if hit['DOI'].lower() != doi.lower():
                raise ValueError('Crossref DOI mismatch')
            info.update(title='; '.join(hit.get('title', [])), year=str(hit.get('published', {}).get('date-parts', [['']])[0][0]), metadata=hit, status='doi_metadata_verified')
    except Exception as e:
        info['error'] = str(e)
    cache.write_text(json.dumps(info, indent=2))
    print(doi, info['status'], info['title'][:90], flush=True)
    return info

dois = set()
for manifest in (base / 'papers').rglob('*.tsv'):
    if 'manifest' not in manifest.name.lower():
        continue
    for row in csv.DictReader(manifest.open(), delimiter='\t'):
        doi = (row.get('doi') or row.get('DOI') or '').strip()
        if doi.startswith('10.'):
            dois.add(doi)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(one, sorted(dois)))
with (base / 'verified_sources.tsv').open('w') as f:
    cols = ['doi', 'title', 'year', 'status', 'path', 'source_url', 'sha256', 'error']
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(results)
