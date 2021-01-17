import sqlite3, random, pickle, os, glob, shlex, subprocess, sys, math
from gensim import models, similarities, corpora
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import islice, tee
from ete3 import NCBITaxa
from functools import reduce
from multiprocessing import Pool, Value, cpu_count
from decimal import Decimal
#This program requires NCBI BLAST suite of programs and HMMER3 to be installed and in the executable path
ncbi = NCBITaxa()
class MyCorpus_stringent(object):
	def __init__(self, docs, fastas, dictionary):
		self.seq = fastas
		self.documents = docs
		self.dictionary = dictionary
	def __iter__(self):
		for t in self.documents:
			species_doc = []
			for p in self.seq[t]:
				trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
				trigram_str = ["".join(el) for el in trigrams]
				species_doc += trigram_str
			yield self.dictionary.doc2bow(species_doc)
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    return None
def parseHmmer(fls):
	fastas = []
	cnt = 0
	for f in fls:
		fasta = f.replace(".tbl", "")
		accs2seq = dict([(rec.id, rec) for rec in SeqIO.parse(fasta, "fasta")])
		with open(f, "r") as handle:
			homologs = set()
			for line in handle:
				if line.startswith("#"):
					continue
				data = line.split()
				homologs.add(data[0])
		diff = set(accs2seq.keys()).difference(homologs)
		if len(diff) >= 100:
			print (len(diff), cnt)
			cnt += 1
			out_file = fasta.split(".")[0] + "_orphan_candidate.fasta"
			with open(out_file, "w") as out:
				SeqIO.write([accs2seq[k] for k in diff], out, "fasta")
			fastas.append(out_file)
	return fastas
def makeReferences():
	ref_taxa = {}
	for rec in SeqIO.parse("500_proteome.fasta", "fasta"):
		data = rec.description.split()
		ref_taxa[data[0]] = data[-1]
	pickle.dump(ref_taxa, open("acc2seq.p", "wb"))
	return (ref_taxa)
def parseBlast(fls, ref_taxa):	
	stringent = []
	for f in fls:
		target = os.path.basename(f).split("_")[0]
		fasta = f.replace(".tab", "")
		accs2seq = dict([(rec.id, rec) for rec in SeqIO.parse(fasta, "fasta")])
		hits = defaultdict(list)
		with open(f, "r") as handle:
			for line in handle:
				data = line.split()
				e_val = float(data[2])
				hits[data[0]].append((ref_taxa[data[1]], e_val))
		orphan_candidates = set()
		homologs = set()
		for acc in hits.keys():
			spread = list(set([el[0] for el in hits[acc]]))
			if (len(spread) == 1 and spread[0] == target):
				orphan_candidates.add(acc) #if the sequence is only being matched within it's original taxa, it becomes an orphan candidate
			else:
				if min([el[1] for el in hits[acc]]) >= 1.0:
					orphan_candidates.add(acc)
			homologs.add(acc)
		diff = set(accs2seq.keys()).difference(homologs)
		candidates = orphan_candidates.union(diff)
		if len(candidates) >= 50:
			out_file = fasta.split(".")[0] + "_orphan_stringent.fasta"
			with open(out_file, "w") as out:
				SeqIO.write([accs2seq[k] for k in candidates], out, "fasta")
			stringent.append(out_file)
	return (stringent, orphan_candidates)
def runBlast(fls):
	blast_results = []
	for cnt, f in enumerate(fls):
		cmd = 'blastp -db 500_proteome.fasta -query {} -num_threads 32 -outfmt "6 qseqid sseqid evalue" -word_size 3 -out {}.tab'.format(f, f)
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.DEVNULL)
		print (p.wait(), cnt, "done!")
		blast_results.append("{}.tab".format(f))
	return blast_results
def runHmmer(fls):
	outFiles = []
	for cnt, f in enumerate(fls):
		cmd = "hmmsearch --cpu 8 --noali --tblout {}.tbl ../PFAM_PROTEOMES/HMM_MODELS/Pfam-A.hmm {}".format(f, f)
		args = shlex.split(cmd)
		p = subprocess.Popen(args, stdout=subprocess.DEVNULL)
		print (p.wait(),cnt, "done!")
		outFiles.append("{}.tbl".format(f))
	return outFiles
def queryDB(cur, spec):
	cur.execute("SELECT acc, fasta_seq FROM proteome WHERE taxon_id=?", (spec,)) 
	proteome = [el for el in cur.fetchall()]
	proteome_fastas = [SeqRecord(Seq(rec[1]), id="{}_{}".format(rec[0], spec)) for rec in proteome]
	f_name = "STRINGENT_ORPHAN/{}_proteome.fasta".format(spec)
	with open(f_name, "w") as proteome_out:
		SeqIO.write(proteome_fastas, proteome_out, "fasta")	
	return f_name
def testSingle(data):
    catch_tx = {}
    d = data[0]
    best_hit = data[1]
    taxonomy = data[2]
    results_single = defaultdict(lambda:{"T":0, "F":0})
    lineage_hit = ncbi.get_lineage(int(best_hit))
    if d == best_hit:
        results_single["strain"]["T"] += 1
        catch_tx["strain"] = d
    else:
        results_single["strain"]["F"] += 1
    for level in taxonomy.keys():
        otu = taxonomy[level]["rnk"]
        if otu in lineage_hit:
            results_single[level]["T"] += 1
            catch_tx[level] = otu
        else:
            results_single[level]["F"] += 1
    pickable_results = {}
    for r in results_single.keys():
        pickable_results[r] = dict(results_single[r])
    return (pickable_results, catch_tx)
def testVote(data):
    d = data[0]
    print ("tu je", d)
    top5 = data[1]
    taxonomy = data[2]
    catch_taxa = defaultdict(dict)
    results_vote = defaultdict(lambda:{"T":0, "F":0})
    OTU_matches = defaultdict(lambda: defaultdict(int))
    lineage_target = ncbi.get_lineage(d)
    for t in top5:
        lineage_hit = ncbi.get_lineage(int(t))
        ranks = ncbi.get_rank(lineage_hit)
        for l in lineage_hit:
            OTU_matches[ranks[l]][l] += 1
    for level in taxonomy.keys():
        otu = taxonomy[level]["rnk"]
        scoring = taxonomy[level]["scoring"]
        if "vote" in scoring:
            if OTU_matches.get(level, None):
                votes = sorted(OTU_matches[level].items(), key=lambda tup: tup[1], reverse=True)
                if len(votes) > 1:
                    if votes[0][1] > votes[1][1]:
                                if otu == votes[0][0]:
                                        results_vote[level]["T"] += 1
                                        catch_taxa[level][d] = otu
                                else:
                                        results_vote[level]["F"] += 1
                    else:
                        results_vote[level]["F"] += 1
                else:
                    if votes[0][0] == otu:
                        results_vote[level]["T"] += 1
                        catch_taxa[level][d] = otu
                    else:
                        results_vote[level]["F"] += 1
            else:
                results_vote[level]["F"] += 1
    pickable_results = {}
    for r in results_vote.keys():
        pickable_results[r] = dict(results_vote[r])
    return (pickable_results, dict(catch_taxa))
def method_selection(docs):
    classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    class_counter = defaultdict(dict)
    method = defaultdict(dict)
    for doc in docs:
        try:
            lineage = ncbi.get_lineage(doc)
            rank = ncbi.get_rank(lineage)
            for l in lineage:
            	if rank[l] in classes:
            		class_counter[l][doc] = rank[l]
        except:
            print ("no lineage data!")
            continue   
    for otu in class_counter.keys():
        total_docs = len(class_counter[otu].keys())
        if total_docs >= 3:
            for doc in class_counter[otu].keys():
                method[doc][class_counter[otu][doc]] = {"scoring":["vote", "single"], "rnk":otu}
        else:
            for doc in class_counter[otu].keys():
                method[doc][class_counter[otu][doc]] = {"scoring":["single"], "rnk":otu}
    pickle_method = {}
    for d in method.keys():
        pickle_method[d] = dict(method[d])
    return pickle_method
def train_LSA_stringent(matrix_data):
	docs = list(matrix_data.keys())
	pickle.dump(docs, open("docs_stringent.p", "wb"))
	dictionary = corpora.Dictionary.load('species_LSA.dict')
	corpus = MyCorpus_stringent(docs, matrix_data, dictionary)
	corpora.MmCorpus.serialize("species_stringent.corp", corpus)
	corpus = corpora.MmCorpus("species_stringent.corp")
	tfidf = models.TfidfModel.load('species_LSA.tfidf')
	corpus_tfidf = tfidf[corpus]
	lsi = models.LsiModel.load("species_LSA.lsi")
	index = similarities.MatrixSimilarity(lsi[corpus_tfidf])
	index.save("species_stringent.index")
	index = similarities.MatrixSimilarity.load("species_stringent.index")
	return (lsi, tfidf, dictionary, index)	
def orphanStringentLSA(orphan_accs):
	matrix_proteins = defaultdict(list)
	for rec in SeqIO.parse("500_proteome.fasta", "fasta"):
		taxa = rec.description.split()[-1]
		if rec.id not in orphan_accs:
			matrix_proteins[taxa].append(str(rec.seq))
	print ("starting LSI construction work!")
	matrices = train_LSA_stringent(matrix_proteins)
	return matrices
def make_bates(f, mixed):
    proteins_sets = {"actual":[str(rec.seq) for rec in SeqIO.parse(f, "fasta")], "mixed":mixed}
    text = {"actual":[], "mixed":[]}
    for sample, proteins in proteins_sets.items():
    	for p in proteins:
        	trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
        	trigram_protein = ["".join(el) for el in trigrams]
        	text[sample] += trigram_protein
    shuffled_aa = [e for e in "".join(text["actual"])]
    random.shuffle(shuffled_aa)
    text_shuffle = ["".join(shuffled_aa[x:x+3]) for x in range(0,len(shuffled_aa),3)]
    alphabet = list(set(shuffled_aa))
    random_str = "".join([random.choice(alphabet) for _ in range(len(shuffled_aa))])
    text_random = [random_str[x:x+3] for x in range(0,len(random_str),3)]
    return (text["actual"], text["mixed"], text_shuffle, text_random)

def mixSamples(fls):
	mixer = {}
	heap = []
	for f in fls:
		target = int(os.path.basename(f).split("_")[0])
		if target in targets:
			proteins = [str(rec.seq) for rec in SeqIO.parse(f, "fasta")]
			mixer[target] = len(proteins)
			heap += proteins
	random.shuffle(heap)
	for k, v in mixer.items():
		mixer[k] = random.sample(heap, v)
	return mixer
def benchmark(fls, matrices):
	classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
	taxonomy_ranking = defaultdict(set)
	catch_vote = defaultdict(lambda: defaultdict(int))
	catch_single = defaultdict(lambda: defaultdict(int))
	tags = ["actual", "mixed", "shuffled", "random"]
	docs = [int(i) for i in pickle.load(open("docs_stringent.p", "rb"))]
	results =  {"actual":(defaultdict(lambda: {"T":0, "F":0}), defaultdict(lambda: {"T":0, "F":0})), 
				"mixed":(defaultdict(lambda: {"T":0, "F":0}), defaultdict(lambda: {"T":0, "F":0})),
	            "shuffled":(defaultdict(lambda: {"T":0, "F":0}), defaultdict(lambda: {"T":0, "F":0})),
	            "random":(defaultdict(lambda: {"T":0, "F":0}), defaultdict(lambda: {"T":0, "F":0}))}
	targets = []
	for f in fls:
		target = int(os.path.basename(f).split("_")[0])
		targets.append(target)
	print (len(targets), "queries")
	mixed_samples = mixSamples(fls)
	method = method_selection(targets)
	lsi, tfidf, dictionary, index = matrices
	cnt = 0
	hits_per_rank = defaultdict(lambda: defaultdict(lambda: {"T":0, "F":0}))
	all_cmp = set()
	confirmed = set()
	confirmed_per_rank = defaultdict(int)
	for target in targets:
		ranks = method[target]
		lineage = ncbi.get_lineage(target)
		rank = ncbi.get_rank(lineage)
		check_ranks = set(rank.values())
		cnt += 1
		sent, mixed, decoy, rand = make_bates(f, mixed_samples[target])
		inverted_ranks = dict([(el[1], el[0]) for el in rank.items()])
		taxonomy_ranking["strain"].add(target)
		for c in classes:
			level_up = inverted_ranks.get(c, None)
			if level_up:
				taxonomy_ranking[c].add(level_up)
		for data_set, species_doc in enumerate([sent, mixed, decoy, rand]):
			data = tags[data_set]
			folded_vec = lsi[tfidf[dictionary.doc2bow(species_doc)]]
			sims = index[folded_vec]
			sorted_sims = sorted(enumerate(sims), key=lambda item: -item[1])
			top5 = [docs[el[0]] for el in [el for el in sorted_sims if docs[el[0]] in existing][0:5]]
			single, catch_single_partial = testSingle((target, top5[0], ranks))
			vote, catch_vote_partial = testVote((target, top5, ranks))
			for pos, r in enumerate((single, vote)):
				for rank in r.keys():
					for match in r[rank].keys():
						results[data][pos][rank][match] += r[rank][match]
			if data == "actual":
				for rnk in single.keys():
					for outcome in single[rnk].keys():
						all_cmp.add(target)
						hits_per_rank[rnk][target][outcome] += single[rnk][outcome]
						if single[rnk]["T"] == 1:
							confirmed_per_rank[rnk] += 1
							confirmed.add(target)
				for rank, tx in catch_single_partial.items():
					catch_single[rank][tx] += 1
				for k, v in catch_vote_partial.items():
					for tx in v.keys():
						catch_vote[k][tx] = catch_vote_partial[k][tx]
		if cnt == 100:
			print ("Threshold of 100 organisms processed achieved!") #we are processing only 100 organisms, however this number can be altered
			break
	results_for_pickle = {}
	for rnk in hits_per_rank.keys():
		tr = 0
		fl = 0
		for tx in hits_per_rank[rnk]:
			tr += hits_per_rank[rnk][tx]["T"]
			fl += hits_per_rank[rnk][tx]["F"]
		print ("True:", tr, "False:", fl)
		print ("#############################")
	print ("SINGLE BEST HIT results")
	for k in results.keys():
		results_for_pickle[k] = (dict(results[k][0]), dict(results[k][1]))
		print (results[k][0])
	print ("#############################")
	print ("VOTING SCHEME results")
	for k in results.keys():
		print (results[k][1])
	pickle.dump(results_for_pickle, open("result_stringent_100_taxa.p", "wb"))
if __name__ == '__main__':
	pfamA_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
	wget.download(pfamA_ftp)
	f = gzip.open('Pfam-A.hmm.gz', 'rb')
	bindata = f.read()
	with open('Pfam-A.hmm', 'w') as f_out:
		f_out.write(bindata)
	acc2seq = makeReferences()
	classes = set(["species", "genus", "family", "order", "class", "phylum", "superkingdom"])
	refseqDB = create_connection("nr_species.db")
	c = refseqDB.cursor()
	species = pickle.load(open("docs_above1000.p", "rb"))
	existing_lineage = defaultdict(list)
	for s in species:
		try:
			#here we test to make sure that we include all kingdoms of life in our sample and that all included organisms have most of the taxonomic ranks
			lineage = ncbi.get_lineage(s)
			ranks = ncbi.get_rank(lineage)
			reverse_ranks = dict([(el[1], el[0]) for el in ranks.items()])
			superkingdom = reverse_ranks.get("superkingdom", None)
			if reverse_ranks:
				if len(set(reverse_ranks.keys()).intersection(classes)) == 7:
			existing_lineage[superkingdom].append(s)
		except:
			continue
	path = os.getcwd()
	new_path = path + "/STRINGENT_ORPHAN"
	try:
		os.mkdir(path)
	except OSError:
		print ("Creation of the directory {} failed".format(new_path))
	selected_species = []
	total_stringent = []
	while True:
		while True:
			for kingdom in existing_lineage.keys():
				random_selection = random.choice(existing_lineage)
				if random_selection not in selected_species:
					selected_species.append(random_selection)
			if len(selected_species) > 100:
				print ("collecting species", iteration)
				break
		f_in = []
		for pos, spec in enumerate(selected_species):
			fl = queryDB(c, spec)
			f_in.append(fl)
		print ("running HMMER")
		f_hmmer = runHmmer(f_in)
		stage1_candidates = parseHmmer(["{}.tbl".format(f) for f in f_in])
		print ("running BLASTp")
		f_blast = runBlast(f_in)
		stage2_candidates, orphan_accs = parseBlast(f_blast, acc2seq)
		total_stringent += stage2_candidates
		if total_stringent >= 100:
			break
	print ("Collected {} organisms represented with {} stringent orphan sequences".format(len(total_stringent)), len(orphan_accs))
	matrices = orphanStringentLSA(orphan_accs)
	benchmark(total_stringent, matrices)

