#parsing the 

filename='result_mc3.csv'

segments={}
enumerates={}

def get_dic_from_line(s):
	tagss, valss= s.split(';')
	tags=tagss.split('-')
	if tags[0]=='':
		return {}
	vals=map(float, valss.split('-'))
	return {t:v for t,v in zip(tags, vals)}

with open(filename, 'r') as f:
	first_line=f.readline().split(',')
	fd={l:d for d,l in enumerate(first_line)}
	index=fd['admixtures']
	print index
	for i,r in enumerate(f.readlines()):
		ad=r.split(',')[index]
		print ad
		news=get_dic_from_line(ad)
		for n, val in news.items():
			if n in segments:
				segments[n].append(val)
				enumerates[n].append(i)
			else:
				segments[n]=[val]
				enumerates[n]=[i]			

print segments
print enumerates

import os
for seg, vals in segments.items():
	with open(os.path.join('segments',seg+'.txt'), 'w') as f:
		f.write(' '.join(map(str,enumerates[seg]))+'\n')
		f.write(' '.join(map(str,vals)))
