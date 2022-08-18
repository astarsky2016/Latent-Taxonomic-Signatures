#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle, random, plyvel
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis 
from ete3 import NCBITaxa
from itertools import islice, tee
from gensim import corpora, models
from scipy.signal import savgol_filter
ncbi = NCBITaxa()
results = pickle.load(open("LTS_TOX/500_results_kingdom_control.p", "rb"))
data = pickle.load(open("LTS_TOX/folded_vecs_500.p", "rb"))
control = data["control"]
print (results.keys())
print (results["sub500-"][1][1][0])
def calculatePrecision(values, classes):
    tags = ["SINGLE", "VOTE"]
    precisions = defaultdict(list)
    for v in values:
        for p, r in enumerate(v):  
            tag = tags[p]
            class_precisions = defaultdict(list)
            for kingdom in r[1][0].keys():
                tmp = r[1][0][kingdom]
                for g in classes:
                    prec = (tmp[g]["T"]/(tmp[g]["T"] + tmp[g]["F"]))*100
                    class_precisions[g].append(prec)
            CTEs = [] #mean
            error = [] #std   
            for g in classes:
                CTEs.append(np.mean(class_precisions[g]))
                error.append(np.std(class_precisions[g]))
            precisions[tag].append((CTEs, error))
    return precisions


# In[20]:


db = plyvel.DB('LTS_TOX/tax2prot_levDB')
subsets = data["500"]
virus_host = pickle.load(open("LTS_TOX/virus2kingdom.p", "rb"))
annotated_viruses = set(virus_host.keys())
vecs = defaultdict(list)
kingdom_names = ncbi.get_taxid_translator([2157, 2, 2759])
for tup in control:
    lineage = ncbi.get_lineage(tup[0])
    ranks = ncbi.get_rank(lineage)
    invert = dict([(v, k) for k, v in ranks.items()])
    if invert["superkingdom"] == 10239:
        exists = set(lineage).intersection(annotated_viruses)
        if exists:
            if len(exists) == 1:
                sk = list(virus_host[list(exists)[0]])[0]
                if sk == 2:
                    vecs["Viruses-bacterial host"].append([el[1] for el in tup[1]])
                if sk == 2759:
                    vecs["Viruses-eukaryotic host"].append([el[1] for el in tup[1]])      
    else:
        vecs[kingdom_names[invert["superkingdom"]]].append([el[1] for el in tup[1]])
for tup in subsets["sub500-"]:
    lineage = ncbi.get_lineage(tup[0])
    ranks = ncbi.get_rank(lineage)
    invert = dict([(v, k) for k, v in ranks.items()])
    if invert["superkingdom"] in [2759, 2157]:
        vecs[kingdom_names[invert["superkingdom"]]].append([el[1] for el in tup[1]])
print ("supplemented")
X_vecs = []
labels_true = []
target_names = []
for pos, k in enumerate(vecs.keys()):
    target_names.append(k)
    if len(vecs[k]) > 500:
        s = 500
    else:
        s = len(vecs[k])
        print (s, k)
    for v in random.sample(vecs[k], s):
        if len(v) == 400:
            X_vecs.append(np.array(v))
            labels_true.append(pos)
lda = LinearDiscriminantAnalysis(n_components=2)
X = np.array(X_vecs)
y = np.array(labels_true)
X_r = lda.fit(X, y).transform(X)
#This is for the benchmarking test and error bars
classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
ordered_keys = ["in50-", "in100-", "in500-", "sub50-", "sub100-", "sub500-"]
all_in = []
all_out = []
for k in ordered_keys:
    v = results[k]
    if "in" in k:
        all_in.append(v)
    else:
        all_out.append(v)
precision_in = calculatePrecision(all_in, classes)
precision_out = calculatePrecision(all_out, classes)   
#extarxt random animal and make signatures for random selections of proteins
matrix_organisms = pickle.load(open("LTS_TOX/species.documents", "rb"))
print (len(matrix_organisms))
kingdom_distro = defaultdict(list)
for o in matrix_organisms:
    lineage = ncbi.get_lineage(o)
    if 2759 in lineage:
        kingdom_distro["Eukaryota"].append(o)
    if 2157 in lineage:
        kingdom_distro["Archaea"].append(o)
    if 2 in lineage:
        kingdom_distro["Bacteria"].append(o)
    if 10239 in lineage:
        kingdom_distro["Viruses"].append(o)
for k, v in kingdom_distro.items():
    print (k, len(v))
selected_animal = random.choice(kingdom_distro["Eukaryota"])
dictionary = corpora.Dictionary.load('LTS_TOX/species.dictionary')
tfidf = models.TfidfModel.load("LTS_TOX/species.tfidf")
proteome = set([protein.decode("utf8") for protein in db.iterator(prefix=bytes("{}_".format(selected_animal), 'utf8'), include_key=False)])
print (len(proteome), "animal proteins total", selected_animal)
db.close()
total_proteome = len(proteome)
proteome_percents = []
p = 0.1
for i in range(10):
    proteome_percents.append((int(p*10), int(p*total_proteome)))
    p += 0.1
vector_representations = []
for partial, p in proteome_percents:
    proteome_cut = random.sample(proteome, p)
    texts = []
    for protein in proteome_cut:
        trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(protein, 3)))))
        texts += ["".join(el) for el in trigrams]
    vector_representations.append((partial, tfidf[dictionary.doc2bow(texts)]))
pickle.dump((vector_representations, selected_animal), open("selected_animal.vecs", "wb"))


# In[26]:


fig, axes = plt.subplots(ncols=2, nrows=3, constrained_layout=True, figsize=(12, 10))
all_ax = axes.ravel()
#plotting benchamrking tests
tags_method = ["SINGLE", "VOTE"]
tags_query = ["in", "out"]
cnt = 0
for pos, prec in enumerate([precision_in, precision_out]):
    for t in tags_method:
        means_50, std_50 = prec[t][0][0], prec[t][0][1]
        means_100, std_100 = prec[t][1][0], prec[t][1][1]
        means_500, std_500 = prec[t][2][0], prec[t][2][1]
        ind = np.arange(len(means_50))  # the x locations for the groups
        width = 0.20  # the width of the bars
        ax = all_ax[cnt]
        rects1 = ax.bar(ind - width/2, means_50, width, yerr=std_50,
                label='-{}50'.format(tags_query[pos]))
        rects2 = ax.bar(ind + width/2, means_100, width, yerr=std_100,
                label='-{}100'.format(tags_query[pos]))
        rects3 = ax.bar(ind + width*1.5, means_500, width, yerr=std_500,
                label='-{}500'.format(tags_query[pos]))
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Precision')
        ax.set_xlabel('Ranks')
        ax.set_title('taxonomy mapping results - {}'.format(t))
        ax.set_xticks(ind + width/2)
        ax.set_xticklabels(classes, rotation = 30)
        ax.grid(True)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        cnt += 1
#plotting LDA of vectors per kingdom
ax = all_ax[cnt]
colors = ["magenta", "blue", "turquoise", "green", "red"]
for color, i, target_name in zip(colors, [0, 1, 2, 3, 4], target_names):
    ax.scatter(X_r[y == i, 0], X_r[y == i, 1], alpha=0.8, color=color, label=target_name)
ax.plot([-14,10],[-6,8], 'black', linewidth=3, linestyle='dashed')
ax.plot([-14,10],[12,-11], 'black', linewidth=3, linestyle='dashed')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_title("Kingdoms of Life LTS signatures: 2-component LDA")
ax.set_ylabel('LD2')
ax.set_xlabel('LD1')
ax.grid(True)

cnt += 1
ax = all_ax[cnt]
increase = 0
step = 10
for v in vector_representations:
    p = v[0]
    xv = [el[0] for el in v[1]]
    yv = [el[1] + increase for el in v[1]]
    y_smooth = savgol_filter(yv, 999, 4)
    ax.plot(xv, y_smooth, label="{}%".format(step))
    step += 10
    increase += 0.005  
handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_title("Numida meleagris proteome coverage")
ax.set_ylabel('TF-IDF weighted frequency')
ax.set_xlabel('3-peptide motifs')
ax.grid(True)


# In[ ]:




