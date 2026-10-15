#!/usr/bin/env python3
"""Retain primary-source full texts and provenance for hypothesis planning only."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'resources/literature/darkmatter-2026-10-15-family-hypotheses'
SOURCES = [
 ('wu2025_pet_prostt5', 'PMC12228687'),
 ('singh2009_light_families', 'PMC2612456'),
 ('pushkarev2018_heliorhodopsins', 'PMC11128463'),
 ('bulzu2022_rhodopsin_kaleidoscope', 'PMC9601134'),
 ('czech2019_ectc', 'PMC6344544'),
 ('pairoh2025_carbonic_anhydrase', 'PMC12635333'),
 ('heinzinger2024_prostt5', 'PMC11616678'),
 ('trehalose2013_metagenome', 'PMC3797794'),
 ('rao2006_t6pp', 'PMC2242562'),
 ('heliorhodopsin2021_context', 'PMC8612252'),
 ('imminger2024_desert', 'PMC11519504'),
 ('xu2021_short_photolyases', 'PMC8476157'),
]

def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'ResearchLiteratureArchive/1.0'})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.read(), response.geturl()

def retain(item):
    key, pmcid = item
    row = {'key': key, 'pmcid': pmcid, 'retrieved_date': '2026-10-15', 'errors': []}
    query = urllib.parse.urlencode({'query': 'PMCID:' + pmcid, 'format': 'json', 'resultType': 'core'})
    try:
        raw, url = fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + query)
        (OUT / (key + '.metadata.json')).write_bytes(raw)
        result = json.loads(raw)['resultList']['result'][0]
        row.update({k: result.get(k, '') for k in ('title', 'doi', 'authorString', 'pubYear', 'journalTitle')})
    except Exception as exc:
        row['errors'].append('metadata: ' + str(exc))
    try:
        raw, url = fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/' + pmcid + '/fullTextXML')
        tree = ET.fromstring(raw)
        if tree.find('.//body') is None:
            raise ValueError('XML has no article body')
        path = OUT / (key + '.xml')
        path.write_bytes(raw)
        blocks = [' '.join(x.itertext()) for x in tree.iter() if x.tag in ('article-title', 'title', 'p')]
        (OUT / (key + '.txt')).write_text('\n\n'.join(blocks))
        row.update({'fulltext_path': str(path.relative_to(ROOT)), 'fulltext_url': url,
                    'sha256': hashlib.sha256(raw).hexdigest(), 'fulltext_status': 'retained_full_xml'})
    except Exception as exc:
        row['errors'].append('fulltext: ' + str(exc))
        row['fulltext_status'] = 'unavailable'
    return row

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(retain, SOURCES))
    # Preserve independently acquired full texts and their retrieval records.
    manifest = OUT / 'manifest.json'
    prior = {r['key']: r for r in json.loads(manifest.read_text())} if manifest.exists() else {}
    for row in rows:
        old = prior.get(row['key'], {})
        if row['fulltext_status'] == 'unavailable' and old.get('fulltext_status', '').startswith('retained'):
            old['retry_errors'] = row['errors']
            row.clear()
            row.update(old)
    rows.extend(r for k, r in prior.items() if k not in {x['key'] for x in rows})
    (OUT / 'manifest.json').write_text(json.dumps(rows, indent=2) + '\n')
    for row in rows:
        print(json.dumps(row))
