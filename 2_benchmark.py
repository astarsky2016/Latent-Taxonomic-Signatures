import sqlite3, glob, logging, pickle, random, shelve, sys, os, re, csv
from Bio import SeqIO
from ete3 import NCBITaxa
from collections import defaultdict
from gensim import corpora, models, similarities
import multiprocessing as mp
from itertools import islice, tee
from gensim.test.utils import datapath
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
ncbi = NCBITaxa()
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return None

def create_table(c, create_table_sql):
    try:
        c.execute(create_table_sql)
    except:
        print("error")
def makeDocs(sizes):
    #Were we create separate LSA train and taxonomy benchmark test sequence datasets for specified sizes for all taxa with more than 1,000 sequences
    docs = pickle.load(open("docs.p", "rb"))
	proteome_limit = 5000
    cmds = []
    for size in sizes:
    	sql_create_documents_table_test = """CREATE TABLE IF NOT EXISTS texts_test{} (
    									taxon_id integer PRIMARY KEY,
    	                                trigrams text NOT NULL
    	                            );""".format(size)
    	cmds.append(sql_create_documents_table_test)
    sql_create_documents_table_train = """CREATE TABLE IF NOT EXISTS texts_train (
                                        taxon_id integer PRIMARY KEY,
                                        trigrams text NOT NULL
                                    );""".format(size)
	conn = create_connection("nr_text.db")
    cmds.append(sql_create_documents_table_train)
	store = conn.cursor()
    for cmd in cmds:
        create_table(store, cmd)
	conn.commit()
    sequencesDB = create_connection("nr_species.db")
    cur = sequencesDB.cursor()
	tables = ["texts_train", "texts_test500", "texts_test100", "texts_test50"]
	cnt = 0
	map_benchmark = {}
    docs_above1000 = []
    size = max(sizes) #we are subsequently undersizing the maximum sized sample
	for id in docs:
		cur.execute("SELECT acc, fasta_seq FROM proteome WHERE taxon_id=?", (id,)) 
		records = dict([(el[0], el[1]) for el in cur.fetchall()])
        if len(records) >= 1000:
    		accs = list(records.keys())
    		if len(accs) >= proteome_limit:
    			selected_proteome = random.sample(accs, proteome_limit)
    		else:
    			selected_proteome = accs
    		query_accs = random.sample(selected_proteome, size)
    		fastas500 = list([records[el] for el in query_accs])
    		subject_accs = set(selected_proteome).difference(set(query_accs))
    		train_fastas = list([records[el] for el in subject_accs])
    		subsample100_accs = random.sample(query_accs, 100)
    		subsample50_accs = random.sample(subsample100_accs, 50)
    		fastas100 = list([records[el] for el in subsample100_accs])
    		fastas50 = list([records[el] for el in subsample50_accs])
    		fastas = [train_fastas, fastas500, fastas100, fastas50]
    		print ("working", cnt)
    		for pos, selection in enumerate(fastas):
    			text = ""
    			for s in selection:
    				trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(s, 3)))))
    				trigram_protein = "|".join(["".join(el) for el in trigrams]) + "|"
    				text += trigram_protein
    			store.execute("INSERT INTO {} VALUES (?, ?)".format(tables[pos]), (id, text[0:-1]))
    			map_benchmark[id] = {"test500":query_accs, "test100":subsample100_accs, "test50":subsample50_accs}
    		conn.commit()
    		cnt += 1
    		print(cnt)
            docs_above1000.append(id)
	print ("texts extracted")
	store.execute("CREATE INDEX taxon_idx_test500 ON texts_test500(taxon_id)")
	store.execute("CREATE INDEX taxon_idx_test100 ON texts_test100(taxon_id)")
	store.execute("CREATE INDEX taxon_idx_test50 ON texts_test50(taxon_id)")
	store.execute("CREATE INDEX taxon_idx_model ON texts_train(taxon_id)")
	conn.commit()
	pickle.dump(map_benchmark, open("benchmark_accs2id.p", "wb"))
    store.close()
    cur.close()
    pickle.dump(docs_above1000, open("docs_above1000.p", "wb"))
	return (docs)
class dictionary_stream(object):
	def __init__(self, db, docs, tag):
		self.store = db
		self.docs = docs
		self.tag = tag
	def __iter__(self):
		for id in self.docs:
			self.store.execute("SELECT trigrams FROM texts_{} WHERE taxon_id=?".format(self.tag), (id,))
			yield self.store.fetchall()[0][0].split("|")
class MyCorpus(object):
	def __init__(self, dctn, docs, db, tag):
		self.dictionary = dctn
		self.store = db
		self.docs = docs
		self.tag = tag
	def __iter__(self):
		for t in self.docs:
			self.store.execute("SELECT trigrams FROM texts_{} WHERE taxon_id=?".format(self.tag), (t,))
			text = self.store.fetchall()[0][0].split("|")
			yield self.dictionary.doc2bow(text)
def loadLSI():
    print ("Loading matrix!")
    loaded_lsi = models.LsiModel.load("species_benchmark.lsi")
    loaded_tfidf = models.TfidfModel.load('species_benchmark.tfidf')
    dictionary = corpora.Dictionary.load('species_benchmark.dict')
    index = similarities.MatrixSimilarity.load("species_benchmark.index")
    return (loaded_lsi, loaded_tfidf, dictionary, index)
def make_bates(text):
    decoy_conc = ""
    for el in text:
        decoy_conc += el
    decoy_chars = [el for el in decoy_conc]
    random.shuffle(decoy_chars)
    decoy_string = "".join(decoy_chars)
    alphabet = list(set(decoy_chars))
    decoy = [decoy_string[x:x+3] for x in range(0,len(decoy_string),3)]
    true_rand_str = "".join([random.choice(alphabet) for _ in range(len(decoy_string))])
    true_rand = [true_rand_str[x:x+3] for x in range(0,len(true_rand_str),3)]
    return (decoy, true_rand)
def testSingle(data): #SINGLE BEST HIT method
    d = data[0]
    best_hit = data[1]
    taxonomy = data[2]
    results_single = defaultdict(lambda:{"T":0, "F":0})
    lineage_hit = ncbi.get_lineage(int(best_hit))
    if d == best_hit:
        results_single["strain"]["T"] += 1
    else:
        results_single["strain"]["F"] += 1
    for level in taxonomy.keys():
        otu = taxonomy[level]["rnk"]
        if otu in lineage_hit:
            results_single[level]["T"] += 1
        else:
            results_single[level]["F"] += 1
    pickable_results = {}
    for r in results_single.keys():
        pickable_results[r] = dict(results_single[r])
    return pickable_results
def testVote(data): #VOTING SCHEME method
    d = data[0]
    top5 = data[1]
    taxonomy = data[2]
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
                                else:
                                        results_vote[level]["F"] += 1
                    else:
                        results_vote[level]["F"] += 1
                else:
                    if votes[0][0] == otu:
                        results_vote[level]["T"] += 1
                    else:
                        results_vote[level]["F"] += 1
            else:
                results_vote[level]["F"] += 1
    pickable_results = {}
    for r in results_vote.keys():
        pickable_results[r] = dict(results_vote[r])
    return pickable_results
def analyze_n_grams(data):
    tags = ["text", "decoy", "rand"]
    text = data["text"]
    decoy, rand = make_bates(text)
    lsi, tfidf, dictionary, index = data["lsi"]
    docs = data["docs"]
    doc = data["doc"]
    ranks = data["taxonomy"]
    allowed = data["allowed"]
    results = {"single":{}, "vote":{}}
    n_grams = [text, decoy, rand]
    for pos, t in enumerate(n_grams):
        run = tags[pos]
        folded_vec = lsi[tfidf[dictionary.doc2bow(t)]]
        sims = index[folded_vec]
        sorted_sims = sorted(enumerate(sims), key=lambda item: -item[1])
        cnt = 0
        top5 = []
        best_hit = None
        for pos, el in enumerate(sorted_sims):
        	if docs[el[0]] in allowed:
        		cnt += 1
        		if cnt == 1:
        			best_hit = docs[el[0]]
        			top5.append(docs[el[0]])
        		else:
        			top5.append(docs[el[0]])
        		if cnt == 5:
        			break
        	if pos == 100:
        		break
        if best_hit:
        	results["single"][run] = testSingle((doc, best_hit, ranks))
        if len(top5) == 5:
        	results["vote"][run]  = testVote((doc, top5, ranks))
    return(results)
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
def benchmark(docs):
    pool = mp.Pool(processes=36) #if you have more or less cpu, adjust this number
    all_results = {}
    texts = create_connection("nr_text.db")
    cur = texts.cursor()
    method = method_selection(docs)
    allowed = list(method.keys())
    random.shuffle(allowed)
    lsi, tfidf, dictionary, index = loadLSI()
    print ("performing query!")
    fractions = [500, 100, 50]          
    for num in fractions:
        fraction_results_text = {"single":defaultdict(lambda: {"T":0, "F":0}), "vote":defaultdict(lambda: {"T":0, "F":0})}
        fraction_results_decoy = {"single":defaultdict(lambda: {"T":0, "F":0}), "vote":defaultdict(lambda: {"T":0, "F":0})}
        fraction_results_random = {"single":defaultdict(lambda: {"T":0, "F":0}), "vote":defaultdict(lambda: {"T":0, "F":0})}
        fraction_results = {"text":fraction_results_text, "decoy":fraction_results_decoy, "rand":fraction_results_random}
        batch = []
        cnt = 0
        print ("working on:", num)
        for doc in allowed:       
            data = {}
            data["lsi"] = (lsi, tfidf, dictionary, index)
            data["doc"] = doc
            data["docs"] = docs
            data["taxonomy"] = method[doc]
            data["allowed"] = allowed
            rec = cur.execute("SELECT trigrams FROM documents_test{} WHERE taxon_id=?".format(num), (doc,))
            real = rec.fetchall()[0][0].split("|")
            data["text"] = real
            batch.append(data)        
            if len(batch) == 100:
                cnt += 100
                print ("pooling")
                results = pool.map(analyze_n_grams, batch)
                for r in results:
                    for method_type in r.keys():
                        for run in r[method_type].keys():
                            for otu in r[method_type][run].keys():
                                for hit in r[method_type][run][otu].keys():
                                    fraction_results[run][method_type][otu][hit] += r[method_type][run][otu][hit]
                print (num, cnt, "finished!")
                batch = []
        if batch:
            results = pool.map(analyze_n_grams, batch)
            for r in results:
                for method_type in r.keys():
                    for run in r[method_type].keys():
                        for otu in r[method_type][run].keys():
                            for hit in r[method_type][run][otu].keys():
                                fraction_results[run][method_type][otu][hit] += r[method_type][run][otu][hit]
        all_results[num] = fraction_results
    pickable_results = {}
    for fraction in all_results.keys():
        pickable_results[fraction] = {}
        fraction_result = all_results[fraction]
        for run in fraction_result.keys():
            pickable_results[fraction][run] = {}
            for method_type in fraction_result[run].keys():
                pickable_results[fraction][run][method_type] = {}
                for otu in fraction_result[run][method_type].keys():
                    pickable_results[fraction][run][method_type][otu] = dict(fraction_result[run][method_type][otu])
    pool.close()
    pool.join()
    pickle.dump(pickable_results, open("results_benchmark.p", "wb"))
    texts.close()
if __name__ == '__main__':
    docs = makeDocs([500, 100, 50])
	conn = create_connection("nr_text.db")
	texts = conn.cursor()
	tag = "train"
	mem_friendly = dictionary_stream(texts, docs, tag)
	dictionary.save("species_benchmark.dict")
	dictionary = corpora.Dictionary.load("species_benchmark.dict")
	corpus = MyCorpus(dictionary, docs, texts, tag)
	corpora.MmCorpus.serialize("species_benchmark.corp", corpus)
	corpus = corpora.MmCorpus("species_benchmark.corp")
	tfidf = models.TfidfModel(corpus)
	tfidf.save("species_benchmark.tfidf")
	tfidf = models.TfidfModel.load("species_benchmark.tfidf")
	corpus_tfidf = tfidf[corpus]
	lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=400)
	lsi.save("species_benchmark.lsi")
	index = similarities.MatrixSimilarity(lsi[corpus_tfidf])
	index.save("species_benchmark.index")
	texts.close()
    benchmark()

