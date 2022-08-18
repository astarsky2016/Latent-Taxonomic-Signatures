"""This is the code that was used for LSA vs BLAST taxonomy benchmarking tests. Here, the smaller 500 protein excluded LSA model has been created and used in comparison
to BLAST alignment-based search with the same species and query proteins. A Diamond-based blast implementation has to be installed and present in the $PATH for this code to execute properly.
Also, a Python module for NCBI Taxonomy handling, ETE3 is a requirement that needs to be present.
"""
import glob, logging, pickle, random, shelve, sys, os, itertools, plyvel, random
from itertools import islice, tee
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from multiprocessing import Pool
from gensim import corpora, models, similarities
from multiprocessing import Pool
ncbi = NCBITaxa()
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

def blastSearch():
	f_q = []
	taxa = defaultdict(list)
	fastas = defaultdict(list)
	for rec in SeqIO.parse("out_500_query.fasta", "fasta"):
		tax = int(rec.id.split("_")[0])
		fastas[tax].append(rec)
	print ("collected", len(fastas))
	for tx in fastas.keys():
		lineage = ncbi.get_lineage(tx)
		ranks = ncbi.get_rank(lineage)
		inverted_rank = dict([(el[1], el[0]) for el in ranks.items()])
		superkingdom = inverted_rank.get("superkingdom", None)
		if superkingdom:
			taxa[superkingdom].append(tx)
	for superkingdom in taxa.keys():
		print (superkingdom, len(taxa[superkingdom]))
		if len(taxa[superkingdom]) >= 1000:
			selected = random.sample(taxa[superkingdom], 1000)
			query = "{}_500_query.fasta".format(superkingdom)
			f_q.append(query)
			with open(query, "w") as f_out:
				sample_1000 = []
				for t in selected:
					sample_1000 += fastas[t]
				SeqIO.write(sample_1000, f_out, "fasta")
			print (query)
		else:
			query = "{}_500_query.fasta".format(superkingdom)
			print (len(taxa[superkingdom]))
			f_q.append(query)
			with open(query, "w") as f_out:
				sample = []
				for t in taxa[superkingdom]:
					sample += fastas[t]
				SeqIO.write(sample, f_out, "fasta")
	for q in f_q:
		cmd = "/home/astar/diamond blastp -c1 -g300 -d main_500.dmnd -q {} --outfmt 6 qseqid sseqid evalue -o {}".format(q, q.replace(".fasta", ".tab"))
		res = os.system(cmd)
	print ("finished!")
def parseBlast():
	mapped_taxonomy = pickle.load(open("mapped_taxonomy.p", "rb"))
	rank_numbers = defaultdict(int)
	for tx in mapped_taxonomy.keys():
		ranks = mapped_taxonomy[tx]
		for v in ranks.values():
			rank_numbers[v] += 1
	kingdoms = [2157, 2, 2759, 10239]
	subsets = {500:defaultdict(list), 100:defaultdict(list), 50:defaultdict(list)}
	queries = defaultdict(dict)
	docs = []
	blast_out = {500:{}, 100:{}, 50:{}}
	for k in kingdoms:
		subsets = defaultdict(dict)
		for rec in SeqIO.parse("{}_500_query.fasta".format(k), "fasta"):
			subsets[int(rec.id.split("_")[0])][rec.id.split()[0]] = 1
			docs.append(int(rec.id.split("_")[0]))
		queries[k][500] = subsets
	print(len(docs))
	for k in queries.keys():
		subsample100 = {}
		subsample50 = {}
		for tx in queries[k][500].keys():
			sample_500 = list(queries[k][500][tx].keys())
			random.shuffle(sample_500)
			subsample100[tx] = dict([(el, 1) for el in sample_500[0:100]])
			subsample50[tx] = dict([(el, 1) for el in sample_500[100:150]])
		queries[k][100] = subsample100
		queries[k][50] = subsample50
	for k in kingdoms:
		with open("ANIMALS/{}_500_query.tab".format(k), "r") as tab_in:
			results_500 = defaultdict(list)
			results_100 = defaultdict(list)
			results_50 = defaultdict(list)
			for line in tab_in:
				data = [el.strip() for el in line.split()]
				if queries[k][100][int(data[0].split("_")[0])].get(data[0], None):
					results_100[int(data[0].split("_")[0])].append((int(data[1].split("_")[0]), float(data[2])))
				if queries[k][50][int(data[0].split("_")[0])].get(data[0], None):
					results_50[int(data[0].split("_")[0])].append((int(data[1].split("_")[0]), float(data[2])))
				results_500[int(data[0].split("_")[0])].append((int(data[1].split("_")[0]), float(data[2])))
			blast_out[500][k] = dict(results_500)
			blast_out[100][k] = dict(results_100)
			blast_out[50][k] = dict(results_50)
		print (k, "done!")
	pickle.dump(blast_out, open("blast_500_results.p", "wb"))
	blast_out = pickle.load(open("blast_500_results.p", "rb"))
	docs = set()
	benchmark_results = {}
	for subgroup in blast_out.keys():
		for k in blast_out[subgroup].keys():
			for tx in blast_out[subgroup][k].keys():
				docs.add(tx)
		single = bestHitBlast(mapped_taxonomy, blast_out[subgroup])
		voting = votingHitsBlast(mapped_taxonomy, rank_numbers, blast_out[subgroup])
		benchmark_results[subgroup] = (single, voting)
	pickle.dump(docs, open("limited_500_species.p", "wb"))
	print (len(docs))
	pickle.dump(benchmark_results, open("500_results_kingdom_blast.p", "wb"))
	
def bestHitBlast(taxonomy, data):
	results_ranks = {2157:defaultdict(lambda: defaultdict(int)), 2:defaultdict(lambda: defaultdict(int)), 2759:defaultdict(lambda: defaultdict(int)), 10239:defaultdict(lambda: defaultdict(int))}
	results_total = defaultdict(int)
	total = 0
	total_tested = 0
	for kingdom in data.keys():
		for target in data[kingdom].keys():
			sorted_hits = sorted(data[kingdom][target], key=lambda tup: tup[1])
			total += 1
			best_hit = sorted_hits[0][0]
			ranks_hit = taxonomy.get(best_hit, None)
			ranks_target = taxonomy.get(target, None)
			if target == best_hit:
				results_total["T"] += 1
			else:
				results_total["F"] += 1
			total_tested += 1
			if ranks_target:
				if ranks_hit:
					for rank in set(ranks_target.keys()).intersection(set(ranks_hit.keys())):
						if ranks_target[rank] == ranks_hit[rank]:
							results_ranks[kingdom][rank]["T"] += 1
						else:
							results_ranks[kingdom][rank]["F"] += 1
	depickled = {}
	for superkingdom in results_ranks.keys():
		depickled[superkingdom] = dict(results_ranks[superkingdom])
	print (results_total, total)
	return ((results_total, total), (depickled, total_tested))
def votingHitsBlast(taxonomy, acceptable, data):
	results_ranks = {2157:defaultdict(lambda: defaultdict(int)), 2:defaultdict(lambda: defaultdict(int)), 2759:defaultdict(lambda: defaultdict(int)), 10239:defaultdict(lambda: defaultdict(int))}
	results_total = defaultdict(int)
	total = 0
	total_voted = 0
	for kingdom in data.keys():
		for target in data[kingdom].keys():
			total += 1
			sorted_hits = sorted(data[kingdom][target], key=lambda tup: tup[1])
			top5 = sorted_hits[0:5]
			ranks_target = taxonomy.get(target, None)
			match = False
			for hit in top5:
				if target == hit[0]:
					match = True
					break
			if match:
				results_total["T"] += 1
			else:
				results_total["F"] += 1
			votes = defaultdict(lambda: defaultdict(int))
			matches = {}
			for hit in top5:
				ranks_hit = taxonomy.get(hit[0], None)
				if ranks_hit:
					for rank in set(ranks_target.keys()).intersection(set(ranks_hit.keys())):
						votes[rank][ranks_hit[rank]] += 1
						matches[rank] = ranks_target[rank]
			if ranks_target:
				if votes:
					total_voted += 1
					for rank in votes.keys():
						sorted_hits = sorted(votes[rank].items(), key=lambda x: x[1], reverse=True)
						if len(sorted_hits) > 1:
							if sorted_hits[0][1] > sorted_hits[1][1]:
								if acceptable[matches[rank]] >= 3:
									if sorted_hits[0][0] == matches[rank]:
										results_ranks[kingdom][rank]["T"] += 1
									else:
										results_ranks[kingdom][rank]["F"] += 1
						else:
							if sorted_hits[0][0] == matches[rank]:
								results_ranks[kingdom][rank]["T"] += 1
							else:
								results_ranks[kingdom][rank]["F"] += 1
	depickled = {}
	for superkingdom in results_ranks.keys():
		depickled[superkingdom] = dict(results_ranks[superkingdom])
	print (results_total, total)
	return ((results_total, total),  (depickled, total_voted))
def bestHit(taxonomy, docs, data):
	results_ranks = {2157:defaultdict(lambda: defaultdict(int)), 2:defaultdict(lambda: defaultdict(int)), 2759:defaultdict(lambda: defaultdict(int)), 10239:defaultdict(lambda: defaultdict(int))}
	results_total = defaultdict(int)
	total = 0
	total_tested = 0
	for tup in data:
		total += 1
		target = tup[0]
		best_hit = docs[tup[1][0][0]]
		ranks_target = taxonomy.get(target, None)
		ranks_hit = taxonomy.get(best_hit, None)
		if target == best_hit:
			results_total["T"] += 1
		else:
			results_total["F"] += 1
		if ranks_target:
			superkingdom = ranks_target["superkingdom"]
			if ranks_hit:
				total_tested += 1
				for rank in set(ranks_target.keys()).intersection(set(ranks_hit.keys())):
					if ranks_target[rank] == ranks_hit[rank]:
						results_ranks[superkingdom][rank]["T"] += 1
					else:
						results_ranks[superkingdom][rank]["F"] += 1
	depickled = {}
	for superkingdom in results_ranks.keys():
		depickled[superkingdom] = dict(results_ranks[superkingdom])
	return ((results_total, total), (depickled, total_tested))
def votingHits(taxonomy, docs, acceptable, data):
	results_ranks = {2157:defaultdict(lambda: defaultdict(int)), 2:defaultdict(lambda: defaultdict(int)), 2759:defaultdict(lambda: defaultdict(int)), 10239:defaultdict(lambda: defaultdict(int))}
	results_total = defaultdict(int)
	total = 0
	total_voted = 0
	for tup in data:
		total += 1
		target = tup[0]
		top5 = [docs[el[0]] for el in tup[1]]
		ranks_target = taxonomy.get(target, None)
		match = False
		for hit in top5:
			if target == hit:
				match = True
				break
		if match:
			results_total["T"] += 1
		else:
			results_total["F"] += 1
		if ranks_target:
			superkingdom = ranks_target["superkingdom"]
			votes = defaultdict(lambda: defaultdict(int))
			matches = {}
			for hit in top5:
				ranks_hit = taxonomy.get(hit, None)
				if ranks_hit:
					for rank in set(ranks_target.keys()).intersection(set(ranks_hit.keys())):
						votes[rank][ranks_hit[rank]] += 1
						matches[rank] = ranks_target[rank]
			if votes:
				total_voted += 1
				for rank in votes.keys():
					sorted_hits = sorted(votes[rank].items(), key=lambda x: x[1], reverse=True)
					if len(sorted_hits) > 1:
						if sorted_hits[0][1] > sorted_hits[1][1]:
							if acceptable[matches[rank]] >= 3:
								if sorted_hits[0][0] == matches[rank]:
									results_ranks[superkingdom][rank]["T"] += 1
								else:
									results_ranks[superkingdom][rank]["F"] += 1
					else:
						if sorted_hits[0][0] == matches[rank]:
							results_ranks[superkingdom][rank]["T"] += 1
						else:
							results_ranks[superkingdom][rank]["F"] += 1
	depickled = {}
	for superkingdom in results_ranks.keys():
		depickled[superkingdom] = dict(results_ranks[superkingdom])
	return ((results_total, total),  (depickled, total_voted))
def count_kingdoms():
	#docs, docs_500 = pickle.load(open("species_removed500.documents", "rb"))
	docs_500 = pickle.load(open("limited_500_species.p", "rb"))
	kingdom_count = defaultdict(int)
	kingdoms = set([2157, 2, 2759, 10239])
	for tax in docs_500:
		lineage = set(ncbi.get_lineage(tax))
		superkingdom = list(lineage.intersection(kingdoms))[0]
		kingdom_count[superkingdom] += 1
	print (kingdom_count)
	print (sum(kingdom_count.values()))
def score():
	docs, docs_500 = pickle.load(open("species_removed500.documents", "rb"))
	mapped_taxonomy = pickle.load(open("mapped_taxonomy.p", "rb"))
	data = pickle.load(open("sims_control_500_limited.p", "rb"))
	rank_numbers = defaultdict(int)
	for tx in mapped_taxonomy.keys():
		ranks = mapped_taxonomy[tx]
		for v in ranks.values():
			rank_numbers[v] += 1
	results = defaultdict(dict)
	for subgroup in data.keys():
		single = bestHit(mapped_taxonomy, docs, data[subgroup])
		voting = votingHits(mapped_taxonomy, docs, rank_numbers, data[subgroup])
		results[subgroup] = (single, voting)
	pickle.dump(results, open("500_results_kingdom_control.p", "wb"))
def foldIt(data):
	index = data[0]
	tx = data[1]
	vec = data[2]
	sims = index[vec]
	sorted_sims = sorted(enumerate(sims), key=lambda item: -item[1])
	return ((tx, sorted_sims[0:5]))
def worker(data):
	lsi = data[2]
	tfidf = data[3]
	dictionary = data[4]
	folded_vec = lsi[tfidf[dictionary.doc2bow(data[1].decode("utf8").split("|"))]]
	return ((data[0], folded_vec))
def loadLSI():
	print ("Loading matrix!")
	loaded_lsi = models.LsiModel.load("species_500_control.lsi")
	loaded_tfidf = models.TfidfModel.load('species_500_control.tfidf')
	dictionary = corpora.Dictionary.load('species_500_control.dictionary')
	index = similarities.MatrixSimilarity.load("species_500_control.index")
	return (loaded_lsi, loaded_tfidf, dictionary, index)
def fillBatch(db, docs_500, lsi, tfidf, dictionary):
	vecs = []
	batch = []
	cnt = 0
	with Pool(32) as pool:
		for tx in docs_500:
			trigrams = db.get(bytes(str(tx), "utf8"))
			batch.append((tx, db.get(bytes(str(tx), "utf8")), lsi, tfidf, dictionary))
			if len(batch) == 1500:
				results = pool.map(worker, batch)
				vecs += [r for r in results]
				cnt += 1500
				print (cnt)
				batch = []
		if batch:
			results = pool.map(worker, batch)
			vecs += [r for r in results]
	return vecs
def benchmark():
	lsi, tfidf, dictionary, index = loadLSI()
	docs, docs_500 = pickle.load(open("species_removed500_control.documents", "rb"))
	test = pickle.load(open("limited_500_species.p", "rb"))
	dbs = ["sub500-", "sub100-", "sub50-", "in500-", "in100-", "in50-"]
	vec_set = {}
	main_db = plyvel.DB('./tax2text_500_levDB')
	total = len(docs_500)
	for subset in dbs:
		db = main_db.prefixed_db(bytes(subset, "utf8"))
		vec_set[subset] = fillBatch(db, docs_500, lsi, tfidf, dictionary)
		print (subset, "done!")
	pickle.dump(vec_set, open("folded_vecs_500_control.p", "wb"))
	vec_set = pickle.load(open("folded_vecs_500_control.p", "rb"))
	similarities = defaultdict(list)
	for subset in vec_set.keys():
		with Pool(32) as pool:
			batch = []
			for tup in vec_set[subset]:
				if tup[0] in test:
					batch.append((index, tup[0], tup[1]))
					if len(batch) == 1000:
						results = pool.map(foldIt, batch)
						similarities[subset] += [r for r in results]
						batch = []
			if batch:
				results = pool.map(foldIt, batch)
				similarities[subset] += [r for r in results]
			print (subset, "done!")
	pickle.dump(dict(similarities), open("sims_control_500_limited.p", "wb"))
	main_db.close()
def extractLineage(docs):
	tax_ranks = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
	mapped_taxonomy = {}
	for d in docs:
		try:
			lineage = ncbi.get_lineage(d)
			rank = ncbi.get_rank(lineage)
			inverted_rank = dict([(v, k) for k,v in rank.items() if v in tax_ranks])
			mapped_taxonomy[d] = inverted_rank
		except:
			print ("error", d)
	print (len(mapped_taxonomy), "mapped")
	pickle.dump(mapped_taxonomy, open("mapped_taxonomy.p", "wb"))
class MyCorpus(object):
	def __init__(self, dictionary, docs, db):
		self.dictionary = dictionary
		self.documents = docs
		self.cursor = db
	def __iter__(self):
		for species in self.documents:
			trigrams = self.cursor.get(bytes(str(species), "utf8"))
			texts = trigrams.decode("utf8").split("|")
			yield self.dictionary.doc2bow(texts)
class MyDictionary(object):
	def __init__(self, docs, db):
		self.documents = docs
		self.cursor = db
	def __iter__(self):
		for species in self.documents:
			trigrams = self.cursor.get(bytes(str(species), "utf8"))
			texts = trigrams.decode("utf8").split("|")
			yield texts
def matrix():
	main_db = plyvel.DB('./tax2text_500_levDB')
	db = main_db.prefixed_db(b'main-')
	docs, docs_500 = pickle.load(open("species_removed500_control.documents", "rb"))
	print ("starting")
	print (len(docs), "species collected in total")
	dictionary_mem = MyDictionary(docs, db)
	dictionary = corpora.Dictionary(s for s in dictionary_mem)
	dictionary.save("species_500_control.dictionary")
	dictionary = corpora.Dictionary.load("species_500_control.dictionary")
	print ("dictionary created")
	corpus = MyCorpus(dictionary, docs, db)
	corpora.MmCorpus.serialize('species_500_control.corpus', corpus)
	corpus = corpora.MmCorpus("species_500_control.corpus")
	print ("corpus created")
	tfidf = models.TfidfModel(corpus)
	corpus_tfidf = tfidf[corpus]
	tfidf.save("species_500_control.tfidf")
	print ("tfidf model created")
	lsi_model = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=400)
	lsi_model.save("species_500_control.lsi")
	corpus_lsi = lsi_model[corpus_tfidf]
	index = similarities.MatrixSimilarity(corpus_lsi)
	index.save("species_500_control.index")
	print ("finished matrix creation")
	main_db.close()
	extractLineage(docs)
if __name__ == '__main__':
	matrix()
	benchmark()
	blastSearch()
	score()
	count_kingdoms()
	parseBlast()

