#!/usr/bin/env python3
"""Extract source identity independently of model-generated download filenames."""
import csv
import hashlib
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import fitz

root = Path('resources/literature/darkmatter-2026-10-13')
with (root / 'source_identity.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['path', 'bytes', 'sha256', 'source_title', 'source_doi', 'source_year', 'first_page_excerpt', 'status'])
    for p in sorted((root / 'papers').rglob('*')):
        if p.suffix not in ('.pdf', '.xml', '.html'):
            continue
        title = doi = year = excerpt = ''
        state = 'identity_requires_review'
        try:
            if p.suffix == '.xml':
                a = ET.parse(p).getroot().find('.//article-meta')
                title = ' '.join(''.join(x.itertext()) for x in a.findall('title-group/article-title'))
                doi = ';'.join(x.text or '' for x in a.findall('article-id') if x.get('pub-id-type') == 'doi')
                year = ';'.join(sorted({x.text for x in a.findall('pub-date/year') if x.text}))
                state = 'source_metadata_extracted'
            elif p.suffix == '.pdf':
                with fitz.open(p) as doc:
                    excerpt = doc[0].get_text()[:2500].replace('\n', ' ')
                    text = '\n'.join(page.get_text() for page in doc)
                    doi = ';'.join(sorted(set(re.findall(r'10\.\d{4,9}/[^\s<>]+', text)))[:12])
                    title = str(doc.metadata.get('title', ''))
                    state = 'pdf_parsed_identity_requires_review'
            w.writerow([str(p), p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest(), title, doi, year, excerpt, state])
        except Exception as e:
            w.writerow([str(p), p.stat().st_size, '', '', '', '', '', 'ERROR:' + str(e)])
