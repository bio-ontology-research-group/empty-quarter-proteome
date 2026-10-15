import csv,os,sys
p='/ibex/scratch/projects/c2014/rob/metagenomics/metadata/sample_metadata.tsv'
w=csv.writer(sys.stdout,delimiter='\t');w.writerow(['sample','role','path','resolved_path','bytes','status'])
for r in csv.DictReader(open(p),delimiter='\t'):
 for role in ('contigs','R1','R2','SE'):
  x=r.get(role,'')
  try:
   s=os.stat(x);w.writerow([r['sample'],role,x,os.path.realpath(x),s.st_size,'present_metadata_only'])
  except OSError as e:w.writerow([r['sample'],role,x,'','','missing:'+str(e)])
