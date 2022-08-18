#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
classes = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
results_cluster = pickle.load(open("LTS_TOX/cluster_orphan_benchmarking.p", "rb"))
all_percents = defaultdict(dict)
for cat in results_cluster.keys():
    print (cat)
    for method in results_cluster[cat].keys():
        all_percents["{}_{}".format(cat, method)]["organism"] = results_cluster[cat][method][0]["T"]/sum(results_cluster[cat][method][0].values())*100
        for c in classes:
            all_percents["{}_{}".format(cat, method)][c] = results_cluster[cat][method][1][c]["T"]/sum(results_cluster[cat][method][1][c].values())*100
print (all_percents)


# In[3]:


import numpy as np
plt.figure(figsize=(10, 5))
plt.ylim(0,100)
cl = ["b", "g", "m", "c"]
col = 0
y_single = None
y_vote = None
for ds, v in all_percents.items():     
    if "cluster" in ds:
        ln = 'dashed'
        mark = "x"
    else:
        ln = 'solid'
        mark = "o"
    x = ["organism"] + classes
    y = [v[c] for c in x]
    if ds == "orphan_single":
        y_single = y[:]
    if ds == "orphan_votes":
        y_vote = y[:]
    plt.plot(x, y, color=cl[col], linewidth=2, linestyle=ln, label=ds, marker=mark, markersize=8)
    col += 1
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
x = ["organism"] + classes
plt.axhline(y=10.34,  xmin=1 , xmax=8, color='black', linestyle='--')
y = [10.34]*8
plt.fill_between(x, y, interpolate=True,  hatch='\\\\', facecolor='grey', alpha=.3)
plt.fill_between(x, y, y_single, color='blue', alpha=.5)
plt.fill_between(x, y_single, y_vote, color='green', alpha=.5)
plt.text(3, 4, 'RESTRICTED', fontsize = 15, bbox = dict(facecolor = 'red', alpha = 0.5))
plt.grid(linestyle='--')
plt.title("Cluster and orphan dataset taxonomy benchmarking")
plt.xlabel('Rank')
plt.ylabel('Precision in %')


# In[ ]:




