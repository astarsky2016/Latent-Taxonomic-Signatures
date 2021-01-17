import glob, pickle, sqlite3, random, sys, os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gensim import models, similarities, corpora
from gensim.matutils import cossim
import pickle, shelve, random, sys, glob, sqlite3, ntpath
from itertools import islice, tee
from collections import defaultdict
from Bio import SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import numpy as np
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except:
        print ("no DB connection")
def testSingle(data):
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

def testVote(data):
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
def make_bates(f):
    #Here we create random and shuffled control samples
    proteins = [str(rec.seq) for rec in SeqIO.parse(f, "fasta")]
    text = []
    for p in proteins:
        trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
        trigram_protein = ["".join(el) for el in trigrams]
        text += trigram_protein
    shuffled_aa = [e for e in "".join(text)]
    random.shuffle(shuffled_aa)
    text_shuffle = ["".join(shuffled_aa[x:x+3]) for x in range(0,len(shuffled_aa),3)]
    alphabet = list(set(shuffled_aa))
    random_str = "".join([random.choice(alphabet) for _ in range(len(shuffled_aa))])
    text_random = [random_str[x:x+3] for x in range(0,len(random_str),3)]
    return (text, text_shuffle, text_random)
def loadLSI(tag):
    docs = pickle.load(open("docs_{}_relaxed.p".format(tag), "rb"))
    loaded_lsi = models.LsiModel.load("species_LSA.lsi")
    loaded_tfidf = models.TfidfModel.load('species_LSA.tfidf')
    dictionary = corpora.Dictionary.load('species_LSA.dict')
    index = similarities.MatrixSimilarity.load("species_{}_relaxed.index".format(tag))
    return (docs, loaded_lsi, loaded_tfidf, dictionary, index)
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
def foldIn_matrixofLife(files, tag, docs):
    texts = []
    docs_specific = []
    for f in files:
        doc = re.findall("(\d+)", os.path.basename(f))
        if int(doc[0]) in docs:
            records = SeqIO.parse(f, "fasta")
            species_doc = []
            for rec in records:
                trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(str(rec.seq), 3)))))
                trigram_str = ["".join(el) for el in trigrams]
                species_doc += trigram_str
            docs_specific.append(int(doc[0]))
        texts.append(species_doc)
    pickle.dump(docs_specific, open("docs_{}_relaxed.p".format(tag), "wb"))
    dictionary = corpora.Dictionary.load('species_LSA.dict')
    corpus = [dictionary.doc2bow(text) for text in texts]
    tfidf = models.TfidfModel.load('species_LSA.tfidf')
    corpus_tfidf = tfidf[corpus]
    lsi = models.LsiModel.load("species_LSA.lsi")
    index = similarities.MatrixSimilarity(lsi[corpus_tfidf])
    index.save("species_{}_relaxed.index".format(tag))   
if __name__ == '__main__':
    results_single = defaultdict(lambda: defaultdict(lambda:{"T":0, "F":0}))
    results_vote = defaultdict(lambda: defaultdict(lambda:{"T":0, "F":0}))
    cluster_sets = glob.glob("RELAXED_ORPHAN/*_cluster.fasta")
    orphan_sets = glob.glob("RELAXED_ORPHAN/*_orphan.fasta")
    tags = ["text", "decoy", "rand"]
    cluster_dict = dict([(int(os.path.basename(f).split("_")[0]), f) for f in cluster_sets]) 
    orphan_dict = dict([(int(os.path.basename(f).split(".")[0]), f) for f in orphan_sets]) 
    joined_dicts = {"cluster":cluster_dict, "orphan":orphan_dict}
    docs = pickle.load(open("docs_relaxed.p", "rb"))
    existing_lineage = set()
    for d in docs:
        try:
            l = ncbi.get_lineage(d)
            existing_lineage.add(d)
        except:
            continue
    #create appropriate LSA indices by folding-in of both cluster and orphan sequence datasets
    foldIn_matrixofLife(cluster_sets, "cluster", docs)
    foldIn_matrixofLife(orphan_sets, "orphan", docs)
    method = method_selection(docs)
    testing = ["orphan", "cluster"]
    for tag in testing:
        docs_specific, lsi, tfidf, dictionary, index = loadLSI(tag)
        for cnt, doc in enumerate(docs_specific):
            sent, decoy, rand = make_bates(joined_dicts[tag][doc])
            ranks = method[doc]
            for pos, t in enumerate([sent, decoy, rand]):
                run = tags[pos]
                folded_vec = lsi[tfidf[dictionary.doc2bow(t)]]
                sims = index[folded_vec]
                sorted_sims = sorted(enumerate(sims), key=lambda item: -item[1])
                top5 = [docs_specific[el[0]] for el in sorted_sims[0:5] if docs_specific[el[0]] in existing_lineage]
                best_hit = docs_specific[sorted_sims[0][0]]
                if best_hit in existing_lineage:
                    single = testSingle((doc, best_hit, ranks))
                    if single:
                        for rank in single.keys():
                            for s in single[rank]:
                                results_single[run][rank][s] += single[rank][s]
                else:
                    print ("skipping")
                vote = testVote((doc, top5, ranks))
                if vote:
                    for rank in vote.keys():
                        for s in vote[rank]:
                            results_vote[run][rank][s] += vote[rank][s]
            print (cnt, "done") 
        print ("finished searching {} LSA vector space".format(tag))
        singles = {}
        for run in results_single.keys():
            otu_dict = {}
            for otu in results_single[run].keys():
                otu_dict[otu] = dict(results_single[run][otu])
            singles[run] = otu_dict
        votes = {}
        for run in results_vote.keys():
            otu_dict = {}
            for otu in results_vote[run].keys():
                otu_dict[otu] = dict(results_vote[run][otu])
            votes[run] = otu_dict
        pickle.dump((singles, votes), open("LSA_{}_vector_space.p".format(tag), "wb"))
        print ("total species processed", len(docs_specific))
