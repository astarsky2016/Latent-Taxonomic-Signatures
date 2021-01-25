from gensim import models, similarities, corpora
import multiprocessing as mp
from gensim.matutils import cossim
import pickle, shelve, random, sys, glob, sqlite3, ntpath, subprocess, os
from itertools import islice, tee, permutations
from collections import defaultdict
from Bio import SeqIO
from ete3 import NCBITaxa
from check_lineage_class import helper_class
ncbi = NCBITaxa()
import numpy as np
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield "".join(lst[i:i + n])
def chunks_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except:
        print ("no DB connection")
def generateFromHMM(d):
    rnk = d[0]
    pf = d[1]
    species = d[2]
    text_family = []
    for i in range(100):
        cmd = ["hmmemit", "-p", "PFAM/HMM_MODELS/{}.hmm".format(pf)]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        emmited = out.decode('utf-8')
        seq_emmited = "".join(emmited.split()[1:]).lower()
        trigrams_em = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(seq_emmited, 3)))))
        text_family += ["".join(el) for el in trigrams_em]
    return (rnk, pf, species, text_family)
def compare_sim():
    #These are the results of both selfish and selfless outcomes for the pfam protein families HMM, actual and mixed datasets
    hmm = pickle.load(open("results_selfish_5_HMM.p", "rb"))
    actual = pickle.load(open("results_selfish_5_actual.p", "rb"))
    mixed = pickle.load(open("results_selfish_5_mixed.p", "rb"))
    data = {"HMM":hmm, "actual":actual, "mixed":mixed}
    for d_Set in data.keys():
        for r in data[d_Set].keys(): 
            selfish = []
            selfless = []
            for pf in data[d_Set][r].keys():
                for res in data[d_Set][r][pf]:
                    selfish.append(res["selfish"])
                    selfless.append(res["selfless"])
            pickle.dump(selfish, open("{}_{}_selfish.p".format(d_Set, r), "wb"))
            pickle.dump(selfless, open("{}_{}_selfless.p".format(d_Set, r), "wb"))
def para_fam(data):
    #this is where we calculate cosine between different combinations of protein family and species vectors
    lsi, tfidf, dictionary, index = data["lsi"]
    docs = data["docs"]
    texts_data = data["texts"]
    rnk = data["rank"]
    pf_family = data["pfam"]
    results = {"rank":rnk, "pfam":pf_family, "hits":[]}
    species = [el[0] for el in texts_data]
    texts = [el[1] for el in texts_data]
    spec_dict = {}
    for pos, s in enumerate(species):
        spec_dict[s] = pos
    SBHs = []
    vecs = []
    for pos, sp in enumerate(species):
        folded_vec = lsi[tfidf[dictionary.doc2bow(texts[pos])]]
        vecs.append(folded_vec)
        sims = index[folded_vec]
        sorted_sims = sorted(enumerate(sims), key=lambda item: -item[1])
        SBHs.append(sorted_sims[0])
    combinations_of_2 = set([tuple(sorted(el)) for el in permutations(species, 2)])
    for combo in combinations_of_2:
        selfish = cossim(vecs[spec_dict[combo[0]]], vecs[spec_dict[combo[1]]])
        selfless1 = SBHs[spec_dict[combo[0]]][1]
        selfless2 = SBHs[spec_dict[combo[1]]][1]
        avg_selfless = (selfless1 + selfless2)/2
        hit1 = docs[SBHs[spec_dict[combo[0]]][0]]
        hit2 = docs[SBHs[spec_dict[combo[1]]][0]]
        if selfish < selfless1:
            if selfish < selfless2:
                results["hits"].append({"query":(combo[0], combo[1]), "subject":(hit1, hit2), "selfish":selfish, "selfless":avg_selfless})
            else:
                results["hits"].append({"query": (combo[0], combo[1]), "subject":(hit1, None), "selfish":selfish, "selfless":avg_selfless})
        else:
            if selfish < selfless2:
                results["hits"].append({"query":(combo[0], combo[1]), "subject":(None, hit2), "selfish":selfish, "selfless":avg_selfless})
            else:
                results["hits"].append({"query":(combo[0], combo[1]), "subject":None, "selfish":selfish, "selfless":avg_selfless})
    return results
def simulate_HMM():
    #This is where we use actual data to create HMM emitted sequences acting as selfish controls, beside mixed and actual sequences collected previously
    parent_dir = od.getcwd()
    path = os.path.join(parent_dir, "PFAM/HMM_MODELS") 
    try: 
        os.mkdir(path) 
    except OSError as error:
        print (error)
    mp.set_start_method("spawn")
    data = pickle.load(open("pfam_hypothesis_actual.p", "rb"))
    cnt = 0
    total = len(data["species"]) + len(data["genus"])
    print ("simulating HMM data")
    tasks = []
    for rnk in data.keys():
        for pf_family in data[rnk].keys():
            species = list(data[rnk][pf_family].keys())
            tasks += [(rnk, pf_family, s) for s in species]
    experiment = {"species":defaultdict(list), "genus":defaultdict(list)}
    cnt = 0
    with mp.get_context("spawn").Pool(32) as pool:
        for r in pool.imap_unordered(generateFromHMM, tasks):
                experiment[r[0]][r[1]].append((r[2], r[3]))
                cnt += 1
    pickable_experiment = {}
    print ("HMM data simulated")
    for rnk in experiment.keys():
        pickable_experiment[rnk] = dict(experiment[rnk])
    pickle.dump(pickable_experiment, open("pfam_hypothesis_HMM.p", "wb"))
def families_test(tag):
    mp.set_start_method("spawn")
    lsi, tfidf, dictionary, index = loadLSI()
    docs = pickle.load(open("../docs.p", "rb"))
    results = defaultdict(lambda: defaultdict(list))
    if tag == "HMM":
        experiment = pickle.load(open("pfam_hypothesis_HMM.p", "rb"))
    if tag == "mixed"
        experiment = pickle.load(open("pfam_hypothesis_mixed.p", "rb"))
    else:
        experiment = pickle.load(open("pfam_hypothesis_actual.p", "rb"))
    lsi_batch = []
    for rnk in experiment.keys():
        for pf in experiment[rnk].keys():
            if tag == "HMM":
                texts = experiment[rnk][pf]
            if tag == "mixed"
                texts = list(experiment[rnk][pf].items())
            else:
                texts == list(experiment[rnk][pf].items())
            info = {"rank":rnk, "pfam":pf, "texts":texts, "lsi":(lsi, tfidf, dictionary, index), "docs":docs}
            lsi_batch.append(info)
    with mp.get_context("spawn").Pool(32) as pool:
        for r in pool.imap_unordered(para_fam, lsi_batch):
            results[r["rank"]][r["pfam"]] += r["hits"]
    pickable_results = {}
    for rnk in results.keys():
        pickable_results[rnk] = dict(results[rnk])
    pickle.dump(pickable_results, open("results_selfish_5_{}.p".format(tag), "wb"))

def loadLSI(size):
    loaded_lsi = models.LsiModel.load("species_LSA.lsi")
    loaded_tfidf = models.TfidfModel.load('species_LSA.tfidf')
    dictionary = corpora.Dictionary.load('species_LSA.dict')
    index = similarities.MatrixSimilarity.load("species_LSA.index")
    return (loaded_lsi, loaded_tfidf, dictionary, index)


if __name__ == '__main__':
    simulate_HMM()
    tags = ["actual", "mixed", "HMM"]
    for tag in tags:
        families_test(tag)
    compare_sim()


