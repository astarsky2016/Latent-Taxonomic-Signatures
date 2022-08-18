"""a copy of NCBI "nr" in FASTA format and mappings of protein accessions to NCBI taxIds have to be downloaded and extracted locally: https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
this piece of code parses NCBI "nr" database protein records in order to map protein sequences to NCBI Taxonomy taxId taxonomic identifiers
it creates a fast key-value storage library, using LevelDB that provides an ordered mapping from string keys to string values. For this it needs Plyvel,a fast and feature-rich Python interface.
its good idea to manually make several copies of acc2tax_levDB LevelDB database (in this code time of execution, 4 were made), once it's been created in order to speed up the NR database creation
"""
import glob, logging, pickle, random, shelve, sys, os, itertools, plyvel, random, csv
from itertools import islice, tee
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from multiprocessing import Pool
ncbi = NCBITaxa()
def splitter(lst, num):
	size = int(len(lst)/num)
	cnt = 0
	for i in range(num):
		cnt = i + 1
		if cnt == 4:
			yield lst[i*size:]
		else:
			yield lst[i*size:cnt*size]
def makeTaxonomyDB_plyvel():
	print ("clean start")
	db = plyvel.DB('./acc2tax_levDB', create_if_missing=True)
	wb = db.write_batch()
	counter = 0
	total = 0
	acc_counter = 0
	with open("prot.accession2taxid.FULL", "r") as mapping:
		print (next(mapping))
		for line in mapping:
			data = [el.strip() for el in line.split()]
			wb.put(bytes(data[0].split(":")[0], 'utf8'), bytes(data[1], 'utf8'))
			counter += 1
			if counter == 100000000:
				wb.write()
				total += counter
				print (total)
				wb.clear()
				counter = 0
	print ("exiting")
	if counter:
		print ("final entries")
		wb.write()
		print ("commited")
		wb.clear()
	db.close()
	print ("done!")
def makeAcc2seq_plyvel():
	print ("clean start")
	db = plyvel.DB('./acc2seq_levDB', create_if_missing=True)
	wb = db.write_batch()
	counter = 0
	total = 0
	acc_counter = 0
	with open("nr", "r") as nr:
		for record in SeqIO.parse(nr, "fasta"):
			seq_fasta = str(record.seq).lower()
			multiple_descriptions = record.description.split('\x01')
			accs = [el.split()[0].strip() for el in multiple_descriptions]
			for acc in accs:
				wb.put(bytes(acc, 'utf8'), bytes(seq_fasta, 'utf8'))
			counter += 1
			if counter == 10000000:
				wb.write()
				total += counter
				print (total)
				wb.clear()
				counter = 0
	print ("exiting")
	if counter:
		print ("final entries")
		wb.write()
		print ("commited")
		wb.clear()
	db.close()
	print ("done!")		
def worker(rec):
	idx = rec[0]
	acc2tax = plyvel.DB('./acc2tax_levDB_{}'.format(idx))
	new_rec = []
	cnt = 0
	for seq in rec[1]:
		for d in seq[1]:
			t = acc2tax.get(bytes(d, 'utf8'))
			if t:
				new_rec.append((bytes("{}_{}".format(int(t), d), 'utf8'), bytes(seq[0], 'utf8')))
				cnt += 1
	acc2tax.close()
	return (new_rec, cnt)
def parseNR():
	print ("clean start")
	print ("parsing NR")
	tax2prot = plyvel.DB('./tax2prot_levDB', create_if_missing=True)
	wb = tax2prot.write_batch()
	with open("nr", "r") as nr:
		tmp = []
		batch = 0
		true_counter = 0
		total = 0
		for record in SeqIO.parse(nr, "fasta"):
			true_counter += 1
			batch += 1
			seq_fasta = str(record.seq).lower()
			multiple_descriptions = record.description.split('\x01')
			accs = [el.split()[0].strip() for el in multiple_descriptions]
			true_acc = [el for el in accs if el]
			if true_acc:
				tmp.append((seq_fasta, accs))
			if batch == 10000000:
				print ("Splitting in parallel")
				with Pool(4) as pool:  #instead of 4, you can put a different number here, depending on how many copies of acc2tax database you make
					dataset = []
					for pos, el in enumerate(splitter(tmp, 4)):
						dataset.append((pos + 1, el))
					results = pool.map(worker, dataset)
					for r in results:
						total += r[1]
						for tup in r[0]:
							wb.put(tup[0], tup[1])
					wb.write()
					wb.clear()
					print(total, true_counter)
					tmp = []
					batch = 0
		if tmp:
			print ("last batch")	
			with Pool(4) as pool:
				dataset = []
				for pos, el in enumerate(splitter(tmp, 4)):
					dataset.append((pos + 1, el))
				results = pool.map(worker, dataset)
				for r in results:
					total += r[1]
					for tup in r[0]:
						wb.put(tup[0], tup[1])
				wb.write()
				wb.clear()
				wb.write()
				wb.clear()
				print(total, true_counter)
	tax2prot.close()
	print ("done!")
def extractTrigrams(data):
	species = data[0]
	sequences = data[1]
	texts = []
	for p in sequences:
		trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
		texts += ["".join(el) for el in trigrams]
	return (bytes(str(species), 'utf8'), bytes("|".join(texts), 'utf8'))
def prepareSpeciesTexts():
	db_text = plyvel.DB('./tax2text_levDB', create_if_missing=True)
	wb = db_text.write_batch()
	upper_limit = {2157: 15000, 10239: 3000, 2: 15000, 2759: 100000}
	lower_limit = {2157: 300, 10239: 10, 2: 300, 2759: 1000}
	db = plyvel.DB('./tax2prot_levDB')
	candidates = pickle.load(open("taxa_candidates_kingdoms_level.p", "rb"))
	docs = []
	for superkingdom in candidates.keys():
		print (superkingdom)
		batch = []
		for tx in candidates[superkingdom]:
			proteome = set()
			upper = upper_limit[superkingdom]
			lower = lower_limit[superkingdom]
			proteome = set([protein.decode("utf8") for protein in db.iterator(prefix=bytes("{}_".format(tx), 'utf8'), include_key=False)])
			size = len(proteome)
			if size >= lower:
				docs.append(tx)
				selection = list(proteome)
				if size > upper:
					selection = random.sample(selection, upper)
				batch.append((tx, selection))
			if len(batch) == 100:
				with Pool(15) as pool:
					results = pool.map(extractTrigrams, batch)
					for r in results:
						wb.put(r[0], r[1])
				wb.write()
				wb.clear()
				print (len(docs))
				batch = []
		if batch:
			with Pool(15) as pool:
				results = pool.map(extractTrigrams, batch)
				for r in results:
					wb.put(r[0], r[1])
			wb.write()
			wb.clear()
	pickle.dump(docs, open("species.documents", "wb"))
	db.close()
	db_text.close()
	print ("texts ready", len(docs), "species")
def recordSpecies():
	docs = pickle.load(open("orphan_species.p", "rb"))
	ofile = open("orphan_organisms.csv", "w")
	writer = csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
	classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
	row = ["organism taxId"] + classes
	writer.writerow(row)
	cnt = 0
	for d in docs:
		lineage = ncbi.get_lineage(d)
		taxid2name = ncbi.get_taxid_translator(lineage)
		ranks = ncbi.get_rank(lineage)
		inverted_ranks = dict([(r, i) for i, r in ranks.items()])
		row = [d]
		for rnk in classes:
			otu = inverted_ranks.get(rnk, "Not assigned")
			row.append(taxid2name.get(otu, "no name assigned"))
		writer.writerow(row)
		cnt += 1
		print (cnt)
	ofile.close()
	print (cnt, "done!")
def select_species():
	cnt = 0
	minimal_proteome = {2157: 300, 10239: 10, 2: 300, 2759: 1000}
	db = plyvel.DB('./acc2tax_levDB_1')
	selected_taxa = defaultdict(int)
	for key, val in db:
		selected_taxa[val] += 1
	print ("DB processed", len(selected_taxa.keys()))
	candidates = defaultdict(list)
	for tx in selected_taxa.keys():
		try:
			lineage = ncbi.get_lineage(int(tx))
			rank = ncbi.get_rank(lineage)
		except:
			continue
		print ("inverting")
		inverted_rank = dict([(v, k) for k,v in rank.items()])
		superkingdom = inverted_rank.get("superkingdom", None)
		species = inverted_rank.get("species", None)
		if superkingdom and species:
			total = selected_taxa[tx]
			if total >= minimal_proteome[superkingdom]:
				cnt += 1
				candidates[superkingdom].append(int(tx))
	pickle.dump(dict(candidates), open("taxa_candidates_kingdoms_level.p", "wb"))
	print ("total", cnt)
	for kingdom in candidates.keys():
		print (kingdom, len(candidates[kingdom]))
	db.close()
if __name__ == '__main__':
	makeTaxonomyDB_plyvel()
	select_species()
	parseNR()
	prepareSpeciesTexts()
	#recordSpecies()

