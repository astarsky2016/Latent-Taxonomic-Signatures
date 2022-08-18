"""This is the code that was used to generate vecotor represenations out of venomous and non-venomous animal proteomes. 
"""
from collections import defaultdict
import pickle, xlsxwriter, sys, random, plyvel, itertools, statistics, glob, os, re
from gensim import corpora, models, similarities
from scipy.spatial import distance
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
ncbi = NCBITaxa()
def loadLSI():
    loaded_lsi = models.LsiModel.load("species.lsi")
    loaded_tfidf = models.TfidfModel.load('species.tfidf')
    dictionary = corpora.Dictionary.load('species.dictionary')
    index = similarities.MatrixSimilarity.load("species.index")
    return (loaded_lsi, loaded_tfidf, dictionary, index)
def addOutgroup(ToxProt, num, taxa_included):
	lsi, tfidf, dictionary, index = loadLSI()
	db = plyvel.DB('./tax2text_levDB')
	ToxProt["bacteria"] = {}
	docs = set(pickle.load(open("species.documents", "rb")))
	bacteria = []
	for d in docs:
		lineage = ncbi.get_lineage(d)
		if 2 in lineage:
			bacteria.append(d)
	selected = random.sample(bacteria, num)
	for cnt, d in enumerate(selected):
		ToxProt["bacteria"][cnt] = d
		taxa_included.add(d)
	print ("collected")
	db.close()
	return (ToxProt, taxa_included)
def distances():
	dist_group = defaultdict(list)
	ToxProt = pickle.load(open("toxprot_clustering.p", "rb"))
	for d1 in ToxProt["venomous animals"]:
		v1 = [el[1] for el in d1]
		for d2 in ToxProt["non_venomous animals"]:
			v2 = [el[1] for el in d2]
			cos_sim = 1 - distance.cosine(v1, v2)
			dist_group["mixed"].append(cos_sim)
	for k, v in ToxProt.items():
		for combo in list(itertools.combinations(list(v), 2)):
			v1 = [el[1] for el in combo[0]]
			v2 = [el[1] for el in combo[1]]
			cos_sim = 1 - distance.cosine(v1, v2)
			dist_group[k].append(cos_sim)
	for k in dist_group.keys():
		print (k, min(dist_group[k]), max(dist_group[k]), sum(dist_group[k])/len(dist_group[k]), statistics.stdev(dist_group[k]))
def indexLineages(docs, existing_proteome):
	filter_out = set()
	for t in existing_proteome:
		lineage = ncbi.get_lineage(t)
		lineage.reverse()
		for f in lineage[0:3]:
			filter_out.add(f)
	all_venom = set([1329911, 6073, 7400, 6854])
	indexed_by_lineage = defaultdict(list)
	for el in docs:
		lineage = ncbi.get_lineage(el)
		if 2759 in lineage:
			if set(lineage).intersection(all_venom):
				continue
			else:
				lineage.reverse()
				for pos, l in enumerate(lineage):
					if l in filter_out:
						break
					else:
						indexed_by_lineage[l].append(el)
						if pos == 11:
							break
	pickle.dump(indexed_by_lineage, open("lineage_indexed.p", "wb"))
	return (indexed_by_lineage, filter_out)
def getFastas():
	db = plyvel.DB('./tax2prot_levDB')
	dataset = pickle.load(open("toxprot_clustering.p", "rb"))
	venomous = []
	non_venomous = defaultdict(list)
	bacteria = []
	cnt = 0
	for cat in dataset.keys():
		if cat == "non_venomous animals":
			for d in range(3):
				for tx in dataset[cat][d].keys():
					for acc, prot in db.iterator(prefix=bytes("{}_".format(tx), 'utf8')):
						non_venomous[d].append(SeqRecord(Seq(prot.decode("utf8")), id="{}_{}".format(tx, acc.decode("utf8")), description=""))
					cnt += 1
					print (cnt)
			print ("collected non venomous")
		if cat == "venomous animals":
			for tx in dataset[cat].keys():
				for acc, prot in db.iterator(prefix=bytes("{}_".format(tx), 'utf8')):
					venomous.append(SeqRecord(Seq(prot.decode("utf8")), id="{}_{}".format(tx, acc.decode("utf8")), description=""))
				cnt += 1
				print (cnt)
			print ("collected venomous")
		if cat == "bacteria":
			for tx in dataset[cat].keys():
				for acc, prot in db.iterator(prefix=bytes("{}_".format(tx), 'utf8')):
					bacteria.append(SeqRecord(Seq(prot.decode("utf8")), id="{}_{}".format(tx, acc.decode("utf8")), description=""))
				cnt += 1
				print(cnt)
	db.close()
	SeqIO.write(bacteria, open("ANIMALS/bacteria_toxprot.fasta", "w"), "fasta")
	for d in non_venomous.keys():
		SeqIO.write(non_venomous[d], open("ANIMALS/non_venomous_toxprot_{}.fasta".format(d+1), "w"), "fasta")
	SeqIO.write(venomous, open("ANIMALS/venomous_toxprot.fasta", "w"), "fasta")
	fastas = glob.glob("ANIMALS/*_toxprot*.fasta")
	for f in fastas:
		if f != "ANIMALS/venomous_toxprot.fasta":
			d = f.replace(".fasta", ".dmnd")
			o = f.replace(".fasta", ".tab")
			cmd = "/home/astar/diamond blastp -d {} -q {} --outfmt 6 qseqid sseqid evalue -o {}".format(d, f, o)
			os.system(cmd)
			print (f)
	print ("singles done")
	for combo in list(itertools.combinations(fastas, 2)):
		d = combo[0].replace(".fasta", ".dmnd")
		p = combo[1]
		o = combo[0].replace(".fasta", "_{}".format(os.path.basename(combo[1]).replace(".fasta", ".tab")))
		cmd = "/home/astar/diamond blastp -d {} -q {} --outfmt 6 qseqid sseqid evalue -o {}".format(d, p, o)
		os.system(cmd)
		print (o)
	for f in glob.glob("ANIMALS/*toxprot*.fasta"):
		if "bacteria" not in f:
			d = "ANIMALS/bacteria_toxprot.dmnd"
			o = f.replace(".fasta", "_bacteria.tab")
			cmd = "/home/astar/diamond blastp -d {} -q {} --outfmt 6 qseqid sseqid evalue -o {}".format(d, f, o)
			os.system(cmd)
			print (o)
def calculateCoverage():
	proteomes = defaultdict(set)
	fastas = glob.glob("ANIMALS/*toxprot*.fasta")
	for f in fastas:
		cnt = 0
		for rec in SeqIO.parse(f, "fasta"):
			tx = rec.id.split("_")[0]
			proteomes[tx].add(rec.id)
	pickle.dump(proteomes, open("proteomes_toxprot.p", "wb"))
	print ("proteomes counted")
	proteomes = pickle.load(open("proteomes_toxprot.p", "rb"))
	blast_outputs = glob.glob("ANIMALS/*toxprot*.tab")
	grouped_results = defaultdict(dict)
	for out in blast_outputs:
		if "bacteria" not in out:
			splitter = out.count('venomous')
			if splitter == 2:
				homologs = defaultdict(set)
				base = os.path.basename(out).replace(".tab", "")
				annotation_type = base.count('non_venomous')
				if annotation_type == 2:
					annotation = "non_venomous"
				else:
					annotation = "mixed"
				with open(out, "r") as f_in:
					for line in f_in:
						data = line.split()
						q = data[0]
						s = data[1]
						e_val = float(data[2])
						if e_val <= 1e-3:
							query_animal = q.split("_")[0]
							subject_animal = s.split("_")[0]
							homologs[query_animal].add(q)
							homologs[subject_animal].add(s)
				grouped_results[annotation][base] = dict(homologs)
				print (base, annotation)
			else:
				homologs = defaultdict(set)
				base = os.path.basename(out).replace(".tab", "")
				if "non_venomous" in base:
					annotation= "non_venomous"
				else:
					annotation = "venomous"
				with open(out, "r") as f_in:
					for line in f_in:
						data = line.split()
						q = data[0]
						s = data[1]
						e_val = float(data[2])
						if e_val <= 1e-3:
							query_animal = q.split("_")[0]
							subject_animal = s.split("_")[0]
							if query_animal != subject_animal:
								homologs[query_animal].add(q)
								homologs[subject_animal].add(s)
				grouped_results[annotation][base] = dict(homologs)
				print (base, annotation)
		else:
			print (out)
			homologs = defaultdict(set)
			base = os.path.basename(out).replace(".tab", "")
			if 'non_venomous' in base:
				annotation = "non_venomous|bacteria"
			else:
				annotation = "venomous|bacteria"
			with open(out, "r") as f_in:
				for line in f_in:
					data = line.split()
					q = data[0]
					e_val = float(data[2])
					if e_val <= 1e-3:
						query_animal = q.split("_")[0]
						homologs[query_animal].add(q)
			grouped_results[annotation][base] = dict(homologs)
	pickle.dump(dict(grouped_results), open("grouped_toxprot.p", "wb"))
	coverages = {}
	for annotation in grouped_results.keys():
		annotation_data = {}
		for dataset in grouped_results[annotation].keys():
			dataset_taxa = {}
			for tx in grouped_results[annotation][dataset].keys():
				total = len(proteomes[tx])
				dataset_taxa[tx] = (len(grouped_results[annotation][dataset][tx])/total)*100
			annotation_data[dataset] = dataset_taxa
		coverages[annotation] = annotation_data
	pickle.dump(coverages, open("coverages_toxprot.p", "wb"))
if __name__ == '__main__':
	getFastas()
	calculateCoverage()
	workbook = xlsxwriter.Workbook('venomous_vs_non_venomous.xlsx')
	worksheet = workbook.add_worksheet()
	tox_data = defaultdict(list)
	with open("toxprot_latest.tab", "r") as tox_in:
		print (next(tox_in))
		for line in tox_in:
			data = [el.strip() for el in line.split("\t") if el.strip()]
			tox_data[int(data[5])].append(data)
	print (len(tox_data))
	existing_proteome = set()
	for el in tox_data.keys():
		for ln in tox_data[el]:
			if len(ln) > 6:
				existing_proteome.add(el)
	print (len(existing_proteome))
	docs = set(pickle.load(open("species.documents", "rb")))
	existing_proteome = set(tox_data.keys()).intersection(docs)
	indexed_by_lineage, filter_out = indexLineages(docs, existing_proteome)
	proteomes = {"venomous animals":{}, "non_venomous animals":{}}
	cnt = 0
	non_venomous_check = set()
	taxa_included = set(tox_data.keys()).intersection(docs)
	for el in existing_proteome:
		lineage = ncbi.get_lineage(el)
		taxa_included.update(set(lineage))
		lineage.reverse()
		ranks = ncbi.get_rank(lineage)
		for pos, l in enumerate(lineage[1:]):
			exists = indexed_by_lineage.get(l, None)
			if exists:
				if len(exists) >= 2 and pos < 12:
					taxa_included.update(set(exists))
					#taxid2name.update(ncbi.get_taxid_translator([el] + exists))
					iterations = 0
					non_venom = []
					while True:
						selected = random.choice(exists)
						if selected not in non_venomous_check:
							non_venom.append(selected)
							non_venomous_check.add(selected)
						else:
							iterations += 1
						if iterations == 1000 or len(non_venom) == 3:
							break
					if len(non_venom) == 3:
						animal_info = []
						for non_venomous_animal in non_venom:
							shared = set(lineage).intersection(set(ncbi.get_lineage(non_venomous_animal)))
							rank_info = None
							for low in lineage:
								if low in shared:
									rank_info = (low, ranks[low])
									break
							animal_info.append((non_venomous_animal, rank_info))
						proteomes["venomous animals"][cnt] = el
						proteomes["non_venomous animals"][cnt] = animal_info
						cnt += 1
						break
	print (cnt)
	proteomes, taxa_included = addOutgroup(proteomes, cnt, taxa_included)
	taxid2name = ncbi.get_taxid_translator(taxa_included)
	columns = [["ToxProt animals", "Non-venomous relatives 1", "Shared Taxonomy level", "Non-venomous relatives 2", "Shared Taxonomy level", "Non-venomous relatives 3", "Shared Taxonomy level", "Outgroup bacteria"]]
	groups = ["venomous animals", "non_venomous animals", "bacteria"]
	for p in range(cnt):
		row = []
		for group in groups:
			if group != "non_venomous animals":
				taxa = proteomes[group][p]
				row.append("{} - NCBI TaxId:{}".format(taxid2name[taxa], taxa))
			else:
				for el in proteomes[group][p]:
					taxa = el[0]
					row.append("{} - NCBI TaxId:{}".format(taxid2name[taxa], taxa))
					shared = el[1]
					row.append("{} - Level:{}".format(taxid2name[shared[0]], shared[1]))
		columns.append(row)
	row = 0
	for c in columns:
		for col, el in enumerate(c):
			worksheet.write(row, col, el)
		row += 1
	workbook.close()
	db = plyvel.DB('./tax2text_levDB')
	lsi, tfidf, dictionary, index = loadLSI()
	vecs = {"venomous animals":{}, "non_venomous animals":defaultdict(dict), "bacteria":{}}
	cnt = 0
	for cat in proteomes.keys():
		for pos in proteomes[cat]:
			if cat != "non_venomous animals":
				txId = proteomes[cat][pos]
				trigrams = db.get(bytes(str(txId), "utf8"))
				texts = trigrams.decode("utf8").split("|")
				folded_vec = lsi[tfidf[dictionary.doc2bow(texts)]]
				vecs[cat][txId] = folded_vec
				cnt += 1
				print (cnt)
			else:
				for p, el in enumerate(proteomes[cat][pos]):
					txId = el[0]
					trigrams = db.get(bytes(str(txId), "utf8"))
					texts = trigrams.decode("utf8").split("|")
					folded_vec = lsi[tfidf[dictionary.doc2bow(texts)]]
					vecs[cat][p][txId] = folded_vec
					cnt += 1
					print(cnt)
pickle.dump(vecs, open("toxprot_clustering.p", "wb"))
db.close()





