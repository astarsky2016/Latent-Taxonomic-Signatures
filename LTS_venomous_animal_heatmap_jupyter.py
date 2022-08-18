#!/usr/bin/env python
# coding: utf-8

# In[99]:


import pickle
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.signal import savgol_filter
from ete3 import NCBITaxa
ncbi = NCBITaxa()
coverages = pickle.load(open("LTS_TOX/coverages_toxprot.p", "rb"))


# In[167]:


import seaborn as sns
import pandas as pd
plt.style.use("seaborn")
y_x = defaultdict(list)
joint = []
venomous_taxa = set()
for annotation in coverages.keys():
    for dataset in coverages[annotation].keys():
        for tx in coverages[annotation][dataset].keys():
            if annotation == "venomous":
                venomous_taxa.add(tx)
            y_x[annotation].append(coverages[annotation][dataset][tx])
            if coverages[annotation][dataset][tx] > 100:
                print (coverages[annotation][dataset][tx])
            joint.append((annotation, coverages[annotation][dataset][tx], tx))
start = 10
above90 = set()
for el in joint:
    if el[0] == "mixed":
        if el[1] >= 0.9 and el[2] in venomous_taxa:
            above90.add(el[2])
print ("evo ih", len(above90))           
percentages = defaultdict(lambda: defaultdict(set))
datasets = set()
for el in joint:
    group = int(el[1]/start)
    percentages[group][el[0]].add(el[2])
    datasets.add(el[0])
datasets = ["venomous", "mixed", "non_venomous", "venomous|bacteria", "non_venomous|bacteria"]
annotations = list(percentages.keys())
columns = []
index = list(percentages.keys())
index.sort(reverse=True)
indices = ["{}%".format(l*10) for l in index]
x_axis = []
for el in datasets:
    if el == "mixed":
        x_axis.append('venomous|non_venomous')
    else:
        x_axis.append(el)
for a in annotations:
    columns.append(dict(zip(datasets,  [len(percentages[a][l]) for l in datasets])))
df = pd.DataFrame(columns, index=indices,columns=datasets)
print (df)
fig, axs = plt.subplots(ncols=2, figsize=(12, 6))
v = np.array(y_x["venomous"])
m = np.array(y_x["mixed"])
n = np.array(y_x["non_venomous"])
vb = np.array(y_x["venomous|bacteria"])
nb = np.array(y_x["non_venomous|bacteria"])
venomous_mean = np.mean(v)
mixed_mean = np.mean(m)
non_venomous_mean = np.mean(n)
venomous_bacteria_mean = np.mean(vb)
non_venomous_bacteria_mean = np.mean(nb)
venomous_std = np.std(v)
mixed_std = np.std(m)
non_venomous_std = np.std(n)
venomous_bacteria_std = np.std(vb)
non_venomous_bacteria_std = np.std(nb)
labels = ['venomous', 'venomous|non_venomous', 'non_venomous', "venomous|bacteria", "non_venomous|bacteria"]
x_pos = np.arange(len(labels))
CTEs = [venomous_mean, mixed_mean, non_venomous_mean, venomous_bacteria_mean, non_venomous_bacteria_mean]
error = [venomous_std, mixed_std, non_venomous_std, venomous_bacteria_std, non_venomous_bacteria_std]      
print (venomous_std)
# Build the plot
ax = axs[0]
ax.bar(x_pos, CTEs,
       yerr=error,
       align='center',
       alpha=0.5,
       ecolor='black',
       capsize=10)
ax.set_ylabel('Percentage of orthologs')
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, rotation=30)
ax.set_title('Overall ortholog coverage')
ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
hmap = sns.heatmap( df, cmap='coolwarm', linewidth = 1 , annot = True, ax=axs[1])
hmap.set_yticklabels(hmap.get_yticklabels(), rotation = 30, fontsize = 8)
hmap.set_xticklabels(x_axis, rotation = 30, fontsize = 8)
#set_xlabels
axs[1].set_title( "Animals sharing significant percentage of orthologs" )
fig.tight_layout()
plt.show()



plt.figure(figsize=(10, 5))
plt.ylim(0,100)
cl = ["r", "g", "b", "c", "m", "y", "k", "w"]
col = 0
for ds, v in y_x.items():
    if "non_" in ds:
        ln = 'dashed'
    else:
        ln = 'solid'
    x = list(range(len(v)))
    #v = savgol_filter(v, 9, 6)
    plt.plot(x, v, color=cl[col], linewidth=2, linestyle=ln, label=ds)
    col += 1
plt.title("BLAST interspecies homology")
plt.legend()
plt.grid(linestyle='--')


# In[ ]:




