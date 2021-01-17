import glob, pickle, sys, os, re
from collections import defaultdict
from Bio import SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()
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
def assesOTU(orphanDB, clusterDB, docs):
    classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    analyze_data = [orphanDB, clusterDB]
    compare = ["orphan", "cluster"]
    results = {}
    for pos, t in enumerate(analyze_data):
        check = defaultdict(lambda:{"T":0, "F":0})
        results_single = defaultdict(lambda:{"T":0, "F":0})
        results_vote = defaultdict(lambda:{"T":0, "F":0})
        results_tuple = [results_single, results_vote]
        algorithm = compare[pos]
        method = method_selection(docs)
        for doc in docs:  
            ranks = method[doc]
            top5 = t[doc]
            best_hit = top5[0]
            hit_results = [testSingle((doc, best_hit, ranks)), testVote((doc, top5, ranks))]
            for pos2, r in enumerate(hit_results):
                if r:
                    for rank in r.keys():
                        for s in r[rank]:
                            results_tuple[pos2][rank][s] += r[rank][s]
        singles = {}
        for tax in results_single.keys():
            singles[tax] = dict(results_single[tax])
        votes = {}
        for tax in results_vote.keys():
            votes[tax] = dict(results_vote[tax])
        results[algorithm] = (singles, votes)
    return results
def makeBlastDB():
    homolog_sets = glob.glob("RELAXED_ORPHAN/*_cluster.fasta")
    orphan_sets = glob.glob("RELAXED_ORPHAN/*_orphan.fasta")
    shared_species = set([el.split(".")[0] for el in orphan_sets]).intersection(set([el.split("_")[0] for el in homolog_sets]))
    tags = ["cluster", "orphan"]
    sets = [homolog_sets, orphan_sets]
    #First we create and index blast databases for orphan and cluster sequences
    for pos, s in enumerate(sets):
    	if pos == 0:
    		db = "RELAXED_ORPHAN/100_random_{}.fasta".format(tags[1])
    	else:
    		db = "RELAXED_ORPHAN/100_random_{}.fasta".format(tags[0])
    	fastas = []
    	for spec in shared_species:
    		if pos == 0:
    			f = "{}_cluster.fasta".format(spec)
    		else:
    			f = "{}_orphan.fasta".format(spec)            
    		fastas += [rec for rec in SeqIO.parse(f, "fasta")]
        with open(db, "w") as output_handle:
            SeqIO.write(fastas, output_handle, "fasta")
        os.system("makeblastdb -in {} -input_type fasta -dbtype prot".format("ORPHANS/100_random_{}.fasta".format(tags[pos])))
    #Next we repeat the process but this time performing blastp search over the databases we created in a bi-directional fashion
    for pos, s in enumerate(sets):
        if pos == 0:
            db = "RELAXED_ORPHAN/100_random_{}.fasta".format(tags[1])
        else:
            db = "RELAXED_ORPHAN/100_random_{}.fasta".format(tags[0])
        for cnt, spec in enumerate(shared_species):
            if pos == 0:
                f = "{}_cluster.fasta".format(spec)
            else:
                f = "{}_orphan.fasta".format(spec)            
    		cmd = 'blastp -db {} -query {} -num_threads 30 -outfmt "6 qseqid sseqid evalue" -word_size 3 -out {}.tab'.format(db, f, f)
    		os.system(cmd)
    		print (cnt, "out of ", len(shared_species), "processed!")
    	print (tags[pos], "done!")
    print ("BLAST finished!")
def parseBlast(blast_out, acc2taxa_subject):
    results = defaultdict(list)
    for f_blast in blast_out:
        correct = a = re.findall("\d+", f_blast)[0]
        for line in open(f_blast, "r"):
            data = line.split()
            query = data[0]
            subject = data[1]
            e_val = float(data[2])
            results[correct].append((acc2taxa_subject[subject], e_val))
    blast_results = {}
    T = 0
    F = 0
    for t in results.keys():
        best_eVal = sorted(results[t], key=lambda tup: tup[1])
        if t in [el[0] for el in best_eVal[0:5]]:
            T += 1
        else:
            F +=1
        blast_results[t] = [el[0] for el in best_eVal[0:5]]
    print ("true", T, "false", F)
    return blast_results
if __name__ == '__main__':
	makeBlastDB()
    blast_out = {}
    acc2taxa_orphanDB = {}
    acc2taxa_clusterDB = {}
    dicts = [acc2taxa_orphanDB, acc2taxa_clusterDB]
    dbs = ["RELAXED_ORPHAN/100_random_orphan.fasta", "RELAXED_ORPHAN/100_random_cluster.fasta"]
    for pos, db in enumerate(dbs):
        for rec in SeqIO.parse(db, "fasta"):
            taxa = rec.description.split()[-1]
            acc = rec.id
            dicts[pos][acc] = taxa
    cluster_sets = glob.glob("RELAXED_ORPHAN/*_cluster.fasta.tab")
    orphan_sets = glob.glob("RELAXED_ORPHAN/*_orphan.fasta.tab")
    data_sets = ["orphanDBvsCluster", "clusterDBvsOrphan"]
    data = [cluster_sets, orphan_sets]
    for pos, d in enumerate(data):
        run = data_sets[pos]
        blast_out[run] = parseBlast(d, dicts[pos])
    pickle.dump(blast_out, open("blast_relaxed_tmp.p", "wb"))
    docs = list(blast_out["orphanDBvsCluster"].keys()) #same species were tested bith as clusters and orphans, and the same will be used later on in LSA based matching, therefore we save it
    pickle.dump([int(el) for el in docs], open("docs_relaxed.p", "wb"))
    res = assesOTU(blast_out["orphanDBvsCluster"], blast_out["clusterDBvsOrphan"], docs)
    pickle.dump(res, open("blast_SBH_VSM_results_relaxed.p", "wb")) #Blast results with SBH nad VSM methods on relaxed orphan dataset









