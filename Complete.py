import os, sys, subprocess, json, pickle
import mechanize
from bs4 import BeautifulSoup
import csv

"""

Complete Program for AMIA Poster
Written by Nikhil Gopal
ngopal@u.washington.edu

"""

class Ontology:
	def __init__(self):
		self.mapping = {
		'coding':'SO_0000316',
		'binding_site':'SO_0000409',
		'stem_loop':'SO_0000313',
		'poly_A':'SO_0000553',
		'primer_binding':'SO_0005850',
		'coding_end':'SO_0000327',
		'coding_start':'SO_0000323',
		'operator':'SO_0000057',
		'promoter':'SO_0000167',
		'conserved':'SO_0000856',
		'polypeptide':'SO_0000104',
		'ribosome_entry_site':'SO_0000139'
		}
		
class Dawg:
	def __init__(self, install_path, endpoint):
		self.install_path = install_path
		self.endpoint = endpoint
	def properties(self):
		return {'install':self.install_path, 'endpoint':self.endpoint} 
	def run_test_query(self):
		test = [self.install_path, 'query', '-c', self.endpoint, '-q', '\"select * {?part a ?type; :name ?name; :description ?desc; :dnaSequence ?dna . ?dna :nucleotides ?nucs } limit 20\"']
	def custom_query(self,query):
		query = '\"'+str(query)+'\"'
		command = [self.install_path, 'query', '-c', self.endpoint, '-f', 'JSON', '-q', query]
		return subprocess.check_output(command)
	def custom_query_dict(self,query,delim="##DELIM##"):
		query = '\"'+str(query)+'\"'
		command = [self.install_path, 'query', '-c', self.endpoint, '-f', 'JSON', '-q', query]
		if delim not in query:
			sys.exit("Please include delimiter ##DELIM## at end of query")
		print ' '.join(command)
		js = subprocess.check_output(command)
		return json.loads(js.split(delim)[1])
	def save_query(self,dict_query,out='query_results.p'):
		pickle.dump( dict_query, open( out, "wb" ) )
		return 0
	def load_query(self, infile='query_results.p'):
		""" Returns a dict """
		return pickle.load( open( infile, "rb" ) )
	def open_results_from_csv(self,path):
		""" HACKED TOGETHER TO READ FROM CSV FROM EVREN """
		dat = []
		first_line = 0
		for i in open(path, 'r').readlines():
			if not first_line:
				first_line = 1
				continue
			else:
				line = i.strip('\r\n').split(',')
				pdict = {'part':line[0]},
				ndict = {'nucs':line[1]},
				tdict = {'type':line[2].split(' ')}
				dat.append([pdict,ndict,tdict])
		return dat
		
""" QUERY FUNCTIONS """
def hammertime(seq_query, file_name, out_dir, type='blastn'):
	try:
		str(seq_query)
		str(file_name)
		str(out_dir)
		#if not os.path.exists(d):
		#	os.mkdir(out_dir)
	except:
		return 'ERROR: Query sequence or file name not a valid string'
	url2 = "http://exploration.weizmann.ac.il/perl/pandatox1_0.mperl?submit=blast"
	px = mechanize.Browser()
	px.set_handle_robots(False)
	px.open(url2)
	px.select_form(name="browse")
	px["sequence"] = str(seq_query)
	px["program"] = [str(type)] # blastn, blastx, and another option
	px["collection"] = ["allGenes"]
	pxres = px.submit()
	pxcont = pxres.read()
	with open(out_dir+os.sep+str(file_name)+".html", "w") as f:
		f.write(pxcont)
	return 0

def hammedFileToDict(file_name):
	soup = BeautifulSoup(open(file_name).read())
	head = [i.get('id') for i in soup.find_all('th')]
	big_data = []
	for k in soup.find_all('tr')[1:]:
		datar = {}
		for n in range(0,17):
			if 'significance' in head[n]:
				datar[head[n]] = float(k.text.encode("ascii","ignore").split('\n')[n])
			else:
				datar[head[n]] = k.text.encode("ascii","ignore").split('\n')[n]
		big_data.append(datar)
	return big_data
	
def hammedFileToDictWithSubtypes(file_name,typ):
	soup = BeautifulSoup(open(file_name).read())
	head = [i.get('id') for i in soup.find_all('th')]
	big_data = []
	for k in soup.find_all('tr')[1:]:
		datar = {}
		for n in range(0,17):
			if 'significance' in head[n]:
				datar[head[n]] = float(k.text.encode("ascii","ignore").split('\n')[n])
			else:
				datar[head[n]] = k.text.encode("ascii","ignore").split('\n')[n]
		datar['subtypes'] = '|'.join(typ)
		big_data.append(datar)
	return big_data
	
def lowestEValue(list_to_be_sorted):
	return sorted(list_to_be_sorted, key=lambda k: k['significance'])

def stringToTriples(line):
	return line.strip('\n').split(' ')
""" END """


""" Read Results """
def lowestEValue(list_to_be_sorted):
	return sorted(list_to_be_sorted, key=lambda k: k['significance'])
""" END """


if __name__ == "__main__":

	""" OBTAIN SEQUENCES FROM SPBPKB  """
	#D = Dawg("/Users/nikhilgopal/Downloads/stardog-1.1.3/stardog", "http://ec2-174-129-47-60.compute-1.amazonaws.com:5822/SBPkb")
	D = Dawg("/Users/nikhilgopal/Downloads/stardog-1.1.3/stardog", '\"snarl://ec2-174-129-47-60.compute-1.amazonaws.com:5820/SBPkb;reasoning=RDFS\"')
	print D.properties()
	print ""
	
	##Run Query ##
	#data = D.custom_query_dict("select ?part ?nucs { {?part a prtype:cds; :name ?name; :description ?desc; :dnaSequence ?dna; pr:status 'Deleted' . ?dna :nucleotides ?nucs} UNION {?part a prtype:coding; :name ?name; :description ?desc; :dnaSequence ?dna; pr:status 'Deleted' . ?dna :nucleotides ?nucs } UNION {?part a so:SO_0000316; :name ?name; :description ?desc; :dnaSequence ?dna; pr:status 'Deleted' . ?dna :nucleotides ?nucs} } ##DELIM##")
	#data = D.custom_query_dict("select distinct ?part ?nucs {?part a so:SO_0000316; :dnaSequence ?dna . ?dna :nucleotides ?nucs} ##DELIM##")
	#D.save_query(data)
	data = D.open_results_from_csv('cds.csv')
		
	##Load Results from Previous Query ##
	#data = D.load_query()
	#for i in data['results']['bindings']:
	#	print i['part']['value'], i['nucs']['value']
		
		
	""" RUN SBPKB DATA AGAINST PANDATOX """
	out = 'html_files'
	pout = 'pickles'
	for i in data:
		name = i[0][0]['part'].split('/')[-1]
		seq = i[1][0]['nucs']
		types = i[2]['type']
		file = out+os.sep+name+'.html'
		pfile = pout+os.sep+name+'.p'
		print file, pfile
		print name
		print seq
		print types
		if os.path.isfile(pfile):
			print "SKIPPED EXISTS", name
			continue
		else:
			print name, hammertime(seq, name, out)
			results = lowestEValue(hammedFileToDictWithSubtypes(file,types)) #lowestEValue(hammedFileToDict(file))
			pickle.dump( results, open( pfile, "wb" ) )
		
	#""" RUN SBPKB DATA AGAINST PANDATOX """
	#out = 'html_files'
	#pout = 'pickles'
	#for i in data['results']['bindings']:
	#	print i['part']['value'], i['nucs']['value']
	#	name = i['part']['value'].split('/')[-1]
	#	seq = i['nucs']['value']
	#	file = out+os.sep+name+'.html'
	#	pfile = pout+os.sep+name+'.p'
	#	print file, pfile
	#	if os.path.isfile(pfile):
	#		print "SKIPPED EXISTS", name
	#		continue
	#	else:
	#		print name, hammertime(seq, name, out)
	#		results = lowestEValue(hammedFileToDict(file))
	#		pickle.dump( results, open( pfile, "wb" ) )
				
	""" READ RESULTS AND SAVE TO CSV """
	csvfile = open('analysis_results.csv', 'w')
	header_written = 0
	for i in os.listdir('pickles'):
		try: #setup to get header output
			data = pickle.load( open('pickles'+os.sep+i, 'rb'))
			data_to_print = lowestEValue(data)[0]
			data_to_print['part'] = str(i)
			if not header_written:
				w = csv.DictWriter(csvfile, data_to_print.keys())
				w.writeheader()
				header_written = 1
			print i, data_to_print
			w.writerow(data_to_print)
		except:
			print 'SKIPPED'
	csvfile.close()