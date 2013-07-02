import os, csv

filename = 'analysis_results.csv'

csv_file = csv.DictReader(open(filename, 'rb'), delimiter=',', quotechar='"')

counts_types = {}
counts_clon = {}

for line in csv_file:
	print line['plasmid_clonability_g']
	try:
		counts_clon[line['plasmid_clonability_g']] += 1
	except:
		counts_clon[line['plasmid_clonability_g']] = 1
	for k in line['subtypes'].split('|'):
		try:
			counts_types[k] += 1
		except:
			counts_types[k] = 1
	

out=open('count_types.txt', 'w')
for i in sorted(counts_types, key=counts_types.get, reverse=True):
	out.write(str(i)+'\t'+str(counts_types[i])+'\n')
	print i, counts_types[i]
out.close()

out=open('count_clon.txt', 'w')
for i in sorted(counts_clon, key=counts_clon.get, reverse=True):	
	out.write(str(i)+'\t'+str(counts_clon[i])+'\n')
	print i, counts_clon[i]
out.close()


csv_file2 = csv.DictReader(open(filename, 'rb'), delimiter=',', quotechar='"')

normal_breakdown = {}
decreased_breakdown = {}
hitchhiker_breakdown = {}
unclonable_breakdown = {}
na_breakdown = {}

for line in csv_file2:
	for k in counts_types:
		if k in line['subtypes'].split('|'):
			if 'normal' in line['plasmid_clonability_g']:
				try:
					normal_breakdown[k] += 1
				except:
					normal_breakdown[k] = 1
			elif 'decreased' in line['plasmid_clonability_g']:
				try:
					decreased_breakdown[k] += 1
				except:
					decreased_breakdown[k] = 1
			elif 'hitchhiker' in line['plasmid_clonability_g']:
				try:
					hitchhiker_breakdown[k] += 1
				except:
					hitchhiker_breakdown[k] = 1
			elif 'unclonable' in line['plasmid_clonability_g']:
				try:
					unclonable_breakdown[k] += 1
				except:
					unclonable_breakdown[k] = 1
			elif 'n/a' in line['plasmid_clonability_g']:
				try:
					na_breakdown[k] += 1
				except:
					na_breakdown[k] = 1
			else:
				print 'FAILED', line['plasmid_clonability_g']
	
print normal_breakdown
print decreased_breakdown
print hitchhiker_breakdown
print unclonable_breakdown

out1=open('normal_breakdown.txt', 'w')
for i in sorted(normal_breakdown, key=normal_breakdown.get, reverse=True):
	out1.write(str(i)+'\t'+str(normal_breakdown[i])+'\n')
	print i, normal_breakdown[i]
out1.close()

out2=open('decreased_breakdown.txt', 'w')
for i in sorted(decreased_breakdown, key=decreased_breakdown.get, reverse=True):
	out2.write(str(i)+'\t'+str(decreased_breakdown[i])+'\n')
	print i, decreased_breakdown[i]
out2.close()

out3=open('hitchhiker_breakdown.txt', 'w')
for i in sorted(hitchhiker_breakdown, key=hitchhiker_breakdown.get, reverse=True):
	out3.write(str(i)+'\t'+str(hitchhiker_breakdown[i])+'\n')
	print i, hitchhiker_breakdown[i]
out3.close()

out4=open('unclonable_breakdown.txt', 'w')
for i in sorted(unclonable_breakdown, key=unclonable_breakdown.get, reverse=True):
	out4.write(str(i)+'\t'+str(unclonable_breakdown[i])+'\n')
	print i, unclonable_breakdown[i]
out4.close()

out5=open('na_breakdown.txt', 'w')
for i in sorted(na_breakdown, key=na_breakdown.get, reverse=True):
	out5.write(str(i)+'\t'+str(na_breakdown[i])+'\n')
	print i, na_breakdown[i]
out5.close()