from pathlib import Path
import json,csv,re,urllib.request,concurrent.futures,hashlib
R=Path(__file__).resolve().parents[1];I=R/'inputs';tree=json.loads((I/'microtrait-hmm-tree.json').read_text());sha=tree['sha'];rules={r['microtrait_rule-name']:r['microtrait_rule-boolean'] for r in csv.DictReader((I/'microtrait_rules.txt').open(),delimiter='\t')}
selected=['hsp10','hsp24','hsp40','hsp60','hsp70','hsp100_key','hsp100_or','lon_or','heatinducibleTF','sodA','ahpCF','TF_redox','trehalose_biosynthesis','betaine_biosynthesis','mannosylglycerate_biosynthesis','proline_biosynthesis','glycine betaine/proline transport system','glycine betaine transporter','choline/glycine/proline betaine transport','proline transport','opu_and','mechanosensitive channel','dpp-1_full','dpp-1_key','dpp-2_full','dpp-2_key','opp-1_full','opp-1_key','aap_full','aap_key','gln_full','glt_full','amt','urease','ammonia assimilation to glutamate','ammonia assimilation to glutamine-GS','nitrogen fixation_full','nitrogen fixation_key','ammonia oxidation to nitrite_full','lycopene biosynthesis','spoIIGA']
paths={Path(t['path']).stem:t['path'] for t in tree['tree'] if t['path'].startswith('data/out/hmm/') and t['path'].endswith('.hmm')}
known={r['microtrait_hmm-name'] for r in csv.DictReader((I/'microtrait_hmm.txt').open(),delimiter='\t')};used={};leaves=set()
def expand(n,stack=()):
 if n in known:leaves.add(n);return repr(n)
 assert n in rules,(n,'missing rule/profile');assert n not in stack,('cycle',n)
 expr=re.sub(r"'([^']+)'",lambda m:'('+expand(m[1],stack+(n,))+')',rules[n]);used[n]=expr;return expr
panels={n:expand(n) for n in selected}
# Keep published rules verbatim; document permissive markers rather than calling them complete pathways.
(I/'panels.json').write_text(json.dumps({'source_repo':'ukaraoz/microtrait','rules_commit':json.loads((I/'microtrait-tree.json').read_text())['sha'],'hmm_commit':sha,'rules':panels,'notes':'Published rule-level endpoints; urease is any-subunit evidence, spoIIGA is a single marker, transporter substrates can overlap. No broad stress tolerance or completed sporulation inference.'},indent=2))
(I/'missing_hmms.json').write_text(json.dumps(sorted(leaves-set(paths))))
leaves=leaves & set(paths)
D=I/'hmms';D.mkdir(exist_ok=True)
def fetch(n):
 p=D/(n+'.hmm');u=f'https://raw.githubusercontent.com/ukaraoz/microtrait-hmm/{sha}/{paths[n]}'
 if not p.exists():p.write_bytes(urllib.request.urlopen(u,timeout=90).read())
 b=p.read_bytes();assert b.startswith(b'HMMER3'),n
 fields={line.split()[0]:line.split()[1:] for line in b.decode().splitlines() if line.startswith(('NAME ','TC ','GA '))}
 assert 'TC' in fields,(n,'missing trusted cutoff')
 return dict(leaf=n,path=paths[n],url=u,sha256=hashlib.sha256(b).hexdigest(),name=fields['NAME'][0],tc=fields['TC'])
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:manifest=list(ex.map(fetch,sorted(leaves)))
(I/'hmm-manifest.json').write_text(json.dumps(manifest,indent=2));(I/'panels.hmm').write_bytes(b''.join((D/(n+'.hmm')).read_bytes() for n in sorted(leaves)))
print('Frozen',len(selected),'rules',len(leaves),'profiles',flush=True)
