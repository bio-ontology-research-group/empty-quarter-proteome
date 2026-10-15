from common import *
import urllib.request,json,subprocess
p=N/'literature';p.mkdir(exist_ok=True)
attempts=[]
for url in ['https://pmc.ncbi.nlm.nih.gov/articles/PMC7049168/','https://europepmc.org/articles/PMC7049168?fulltext=1','https://juser.fz-juelich.de/record/875138/files/4131.full.pdf']:
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
        data=urllib.request.urlopen(req,timeout=90).read()
        if data.startswith(b'%PDF'):
            target=p/'kovalev2020_heliorhodopsin.pdf';target.write_bytes(data)
            subprocess.run(['pdftotext','-layout',str(target),str(target.with_suffix('.txt'))],check=True)
            assert '1915888117'in target.with_suffix('.txt').read_text()
        else:
            content=data.decode();assert '1915888117'in content and 'High-resolution structural insights'in content and len(content)>50000
            target=p/'kovalev2020_heliorhodopsin.html';target.write_bytes(data)
        (p/'kovalev2020_source.json').write_text(json.dumps({'url':url,'doi':'10.1073/pnas.1915888117','PDB':'6SU3','sha256':sha(target),'job':os.environ['SLURM_JOB_ID'],'prior_attempts':attempts},indent=2));break
    except Exception as e:attempts.append({'url':url,'error':repr(e)})
else:
    (p/'kovalev2020_unavailable.json').write_text(json.dumps(attempts,indent=2));raise RuntimeError(attempts)
