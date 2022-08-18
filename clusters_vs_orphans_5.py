#presuming you have wget installed and in your PATH
"""This is the code that uses previously created LSA species model in order to establish connection between taxonomically restricted proteins (orphans),
and the rest of the "shared" proteome. For this purpose, NCBI Clusters database has been used to provide Clusters dataset, which was used to also create candidate 
orphan protein dataset, so that Clusters can be compared to orphans in a taxonomy benchmarking scenario.
"""
import glob, pickle, sqlite3, random, sys, os, wget, plyvel
from collections import defaultdict
from Bio import SeqIO, Entrez
from ete3 import NCBITaxa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import islice, tee
from gensim import corpora, models, similarities
Entrez.email = "astar@pbf.hr"
ncbi = NCBITaxa()
def extractLineage():
	docs = pickle.load(open("species.documents", "rb"))
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
def cluster_information(outdir):
	cluster_tags = ["CHL", "CLSC", "CLSM", "CLSN", "CLSP", "CLSZ", "MTH", "PCLA", "PHA", "PLN", "PTZ"]
	cluster_info = defaultdict(set)
	protein_clusters = []
	#Download CLuster info from NCBI ftp
	for tag in cluster_tags:
		protein_cluster = "{}_proteins.txt".format(tag)
		url = "https://ftp.ncbi.nih.gov/genomes/CLUSTERS/" + protein_cluster
		wget.download(url, out=outdir)
		protein_clusters.append(protein_cluster)
	for f in protein_clusters:
		with open(f, "r") as handle:
			next(handle)
			for line in handle:
				data = line.split()
				cluster = data[0]
				acc = data[1].split(".")[0]
				taxid = int(data[-2])
				cluster_info[taxid].add(acc)
	print (len(cluster_info))
	sorted_clusters = sorted(cluster_info.items(), key=lambda tup: len(tup[1]), reverse = True)
	pickle.dump(cluster_info, open("cluster_info.p", "wb"))
	return cluster_info
def extract_clusters():
	lower_limit = {2157: 300, 10239: 10, 2: 300, 2759: 1000}
	docs = pickle.load(open("species.documents", "rb"))
	db = plyvel.DB('./tax2prot_levDB')
	path = os.getcwd()
	cnt = 0
	new_path = os.path.join(path + "NCBI_CLUSTERS")
	try:
		os.mkdir(new_path)
	except OSError:
		print ("Creation of the directory {} failed".format(new_path))
	else:
		print ("Successfully created the directory {}".format(new_path))
	cnt = 0
	cluster_info = cluster_information(new_path)
	for taxon in cluster_info.keys():
		print (type(taxon))
		cluster_accs = cluster_info[taxon]
		if taxon in docs:
			lineage = ncbi.get_lineage(taxon)
			rank = ncbi.get_rank(lineage)
			inverted_rank = dict([(v, k) for k,v in rank.items()])
			superkingdom = inverted_rank["superkingdom"]
			acc2prot = {}
			for key, protein in db.iterator(prefix=bytes("{}_".format(tx), 'utf8')):
				acc = "".join(key.decode("utf8").split("_")[1:])
				protein = protein.decode("utf8")
				acc2prot[acc] = protein
			proteome_acc = set(acc2prot.keys())
			obtained_cluster = cluster_accs.intersection(proteome_acc)
			no_cluster = proteome_acc.difference(cluster_accs)
			if len(obtained_cluster) >= lower_limit[superkingdom]:
				if len(no_cluster) >= lower_limit[superkingdom]:
					cnt += 1
					print(cnt)
					cluster_fastas = [SeqRecord(Seq(acc2prot[ac]), id=ac, description="taxa {}".format(taxon)) for ac in obtained_cluster]
					no_cluster_fastas = [SeqRecord(Seq(records[ac]), id=ac, description="taxa {}".format(taxon)) for ac in no_cluster]
					with open("CLUSTERS_FASTA/{}_cluster.fasta".format(taxon), "w") as out:
						SeqIO.write(cluster_fastas, out, "fasta")	
					with open("CLUSTERS_FASTA/{}_out.fasta".format(taxon), "w") as out:
						SeqIO.write(no_cluster_fastas, out, "fasta")	
	db.close()
def diamond_screen():
	cnt = 0
	taxons, queries = pickle.load(open("core_set.p", "rb"))
	for class_ in queries.keys():
		cnt_join = 0
		tmp_file = "/disk2/NCBI/CLUSTERS_FASTA/{}_tmp.fasta".format(class_)
		with open(tmp_file, "a") as tmp_join:
			for tx in queries[class_]:
				tx_fasta = "/disk2/NCBI/CLUSTERS_FASTA/{}_out.fasta".format(tx)
				sequences = []
				for rec in SeqIO.parse(tx_fasta, "fasta"):
					rec.id = "{}_{}".format(tx, rec.id)
					sequences.append(rec)
				SeqIO.write(sequences, tmp_join, "fasta")
		cnt_join += 1
		print ("joined: ", cnt_join)
		bl_out = "/disk2/NCBI/CLUSTERS_FASTA/SUBSAMPLED/{}.tab".format(class_)
		ref = "/disk2/NCBI/CLUSTERS_FASTA/SUBSAMPLED/{}.dmnd".format(class_)
		cmd = '/home/astar/diamond blastp -d nr.dmnd -q {} -e 0.1 --very-sensitive --outfmt 6 qseqid sseqid evalue -o {} >/dev/null 2>&1'.format(tmp_file, bl_out)
		os.system(cmd)
		cnt += 1
		os.remove(tmp_file)
		print (cnt)
def bestHit(taxonomy, docs, hits):
	results_ranks = defaultdict(lambda: defaultdict(int))
	results_total = defaultdict(int)
	total_tested = 0
	for target in hits.keys():
		best_hit = docs[hits[target][0][0]]
		ranks_target = taxonomy.get(target, None)
		ranks_hit = taxonomy.get(best_hit, None)
		if target == best_hit:
			results_total["T"] += 1
		else:
			results_total["F"] += 1
		if ranks_target:
			if ranks_hit:
				total_tested += 1
				for rank in set(ranks_target.keys()).intersection(set(ranks_hit.keys())):
					if ranks_target[rank] == ranks_hit[rank]:
						results_ranks[rank]["T"] += 1
					else:
						results_ranks[rank]["F"] += 1
	depickled = {}
	for rnk in results_ranks.keys():
		depickled[rnk] = dict(results_ranks[rnk])
	return (results_total, depickled, total_tested)
def votingHits(taxonomy, docs, acceptable, hits):
	results_total = defaultdict(int)
	results_ranks = defaultdict(lambda: defaultdict(int))
	total_voted = 0
	for target in hits.keys():
		top5 = [docs[el[0]] for el in hits[target]]
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
									results_ranks[rank]["T"] += 1
								else:
									results_ranks[rank]["F"] += 1
					else:
						if sorted_hits[0][0] == matches[rank]:
							results_ranks[rank]["T"] += 1
						else:
							results_ranks[rank]["F"] += 1
	depickled = {}
	for rnk in results_ranks.keys():
		depickled[rnk] = dict(results_ranks[rnk])
	return (results_total, depickled, total_voted)
def loadLSI():
	print ("Loading matrix!")
	loaded_lsi = models.LsiModel.load("species.lsi")
	loaded_tfidf = models.TfidfModel.load('species.tfidf')
	dictionary = corpora.Dictionary.load('species.dictionary')
	index = similarities.MatrixSimilarity.load("species.index")
	return (loaded_lsi, loaded_tfidf, dictionary, index)
def queryMatrix():
	docs = pickle.load(open("species.documents", "rb"))
	mapped_taxonomy = pickle.load(open("mapped_taxonomy.p", "rb"))
	rank_numbers = defaultdict(int)
	for tx in mapped_taxonomy.keys():
		ranks = mapped_taxonomy[tx]
		for v in ranks.values():
			rank_numbers[v] += 1
	lsi, tfidf, dictionary, index = loadLSI()
	clusters_fastas = set([int(os.path.basename(f).split("_")[0]) for f in glob.glob("CLUSTERS_FASTA/*.fasta")])
	orphans_fastas = set([int(os.path.basename(f).split("_")[0]) for f in glob.glob("ORPHANS_FASTA/*.fasta")])
	shared = clusters_fastas.intersection(orphans_fastas)
	print (len(shared), len(clusters_fastas), len(orphans_fastas))
	query_cluster = defaultdict(set)
	query_orphan = defaultdict(set)
	for tx in shared:
		orphans = "ORPHANS_FASTA/{}_true_orphan.fasta".format(tx)
		clusters = "CLUSTERS_FASTA/{}_cluster.fasta".format(tx)
		for rec in SeqIO.parse(orphans, "fasta"):
			query_orphan[tx].add(str(rec.seq))
		for rec in SeqIO.parse(clusters, "fasta"):
			query_cluster[tx].add(str(rec.seq))
	cnt = 0
	hits = {"orphan":{}, "cluster":{}}
	range_len = defaultdict(list)
	f_orphan = open("orphans.fasta", "a")
	f_clusters = open("clusters.fasta", "a")
	for tx in query_orphan.keys():
		limit_orphan = len(query_orphan[tx])
		limit_cluster = len(query_cluster[tx])
		if limit_orphan >= 100 and limit_cluster >= 100:
			normalized_orphan = query_orphan[tx]
			normalized_cluster = query_cluster[tx]
			range_len["orphan"].append(limit_orphan)
			cnt += 1
			texts = []
			if limit_orphan > limit_cluster:
				normalized_orphan = random.sample(list(query_orphan[tx]), limit_cluster)
			else:
				normalized_cluster = random.sample(list(query_cluster[tx]), limit_orphan)
			sequences_orphan = []
			num = 1
			for seq in normalized_orphan:
				sequences_orphan.append(SeqRecord(Seq(seq), id="{}_{}".format(tx, num), description="orphan protein"))
				trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(str(seq), 3)))))
				texts += ["".join(el) for el in trigrams]
				num += 1
			SeqIO.write(sequences_orphan, f_orphan, "fasta")
			folded_orphan = lsi[tfidf[dictionary.doc2bow(texts)]]
			hits_orphan = sorted(enumerate(index[folded_orphan]), key=lambda item: -item[1])[0:5] 
			texts = []
			range_len["cluster"].append(limit_cluster)
			sequences_cluster = []
			num = 1
			for seq in normalized_cluster:
				sequences_cluster.append(SeqRecord(Seq(seq), id="{}_{}".format(tx, num), description="clusters protein"))
				trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(str(seq), 3)))))
				texts += ["".join(el) for el in trigrams]
				num += 1
			SeqIO.write(sequences_cluster, f_clusters, "fasta")
			folded_cluster = lsi[tfidf[dictionary.doc2bow(texts)]]
			hits_cluster = sorted(enumerate(index[folded_cluster]), key=lambda item: -item[1])[0:5] 
			print (cnt)
			hits["cluster"][tx] = hits_cluster
			hits["orphan"][tx] = hits_orphan
	print (cnt, "done!")
	pickle.dump((hits,dict(range_len)), open("compare_hits.p", "wb"))
	f_clusters.close()
	f_orphan.close()
	hits, range_len = pickle.load(open("compare_hits.p", "rb"))
	benchmarking = defaultdict(dict)
	for k, v in hits.items():
		single = bestHit(mapped_taxonomy, docs, v)
		voting = votingHits(mapped_taxonomy, docs, rank_numbers, v)
		benchmarking[k]["single"] = single
		benchmarking[k]["votes"] = voting
	pickle.dump(benchmarking, open("cluster_orphan_benchmarking.p", "wb"))
if __name__ == '__main__':
	extractLineage()
	queryMatrix()




