#!/usr/bin/env python
# coding: utf-8

# In[530]:


import plotly.express as px
import plotly.offline as pyo
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.patches as mpatches

# Set notebook mode to work in offline
pyo.init_notebook_mode()
tsne = TSNE(n_components=2)
pca = PCA(n_components=2)
scaler=StandardScaler()
toxprot_data = pickle.load(open("LTS_TOX/toxprot_clustering.p", "rb"))
X_vecs = defaultdict(list)
labels_true = defaultdict(list)
annotations = {"venomous animals":1, "non_venomous animals":2,"bacteria":0}
for dataset in range(3):
    for k in toxprot_data.keys():
        if k == "non_venomous animals":
            taxa = toxprot_data[k][dataset].keys()
        else:
            taxa = toxprot_data[k].keys()
        for tx in taxa:
            #if tx in [9430, 7924]:
            #    continue
            if k == "non_venomous animals":
                folded = toxprot_data[k][dataset][tx]
            else:
                folded = toxprot_data[k][tx]
            v = [el[1] for el in folded]
            X_vecs[dataset].append(np.array(v))
            labels_true[dataset].append(annotations[k])
lda = LinearDiscriminantAnalysis(
    solver='svd', #{‘svd’, ‘lsqr’, ‘eigen’}, default=’svd’
    n_components=2, #int, default=None
)

fig, axes = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(12, 10))
all_ax = axes.ravel()
colour_scheme = ["blue", "cyan", "magenta"]
for d in range(3):
    X = np.array(X_vecs[d])
    y = np.array(labels_true[d])
    X_r = lda.fit(X, y).transform(X)
    venom_dots = np.array([X_r[p] for p, v in enumerate(labels_true[d]) if v == 1])
    non_venom_dots = np.array([X_r[p] for p, v in enumerate(labels_true[d]) if v == 2])
    bacteria_dots = np.array([X_r[p] for p, v in enumerate(labels_true[d]) if v == 0])
    all_ax[d].scatter(non_venom_dots[:, 0], non_venom_dots[:, 1], alpha=0.8, color=colour_scheme[d], label="Non-venomous {}".format(d+1))
    all_ax[d].scatter(bacteria_dots[:, 0], bacteria_dots[:, 1], alpha=0.8, color="green", label="Bacteria")
    all_ax[d].scatter(venom_dots[:, 0], venom_dots[:, 1], alpha=0.8, color="red", label="Venomous")
    all_ax[d].grid(True)
    all_ax[d].set_title("Dataset {}".format(d + 1))
    all_ax[d].set_xlabel("LD1")
    all_ax[d].set_ylabel("LD2")
    all_ax[d].legend()
fig.suptitle('Venomous vs Non_venomous animal LTS signatures', fontsize=16)
# Select the model and its parameters
# Fit transform the data
X_combined = []
labels_combined = []
labels_dataset = []
for d in range(3):
    vecs = X_vecs[d]
    labels = labels_true[d]
    if d == 0:
        bacteria = [vecs[p] for p, l in enumerate(labels) if l == 0]
        venomous = [vecs[p] for p, l in enumerate(labels) if l == 1]
        non_venomous = [vecs[p] for p, l in enumerate(labels) if l == 2]
        X_combined += bacteria
        X_combined += venomous
        X_combined += non_venomous
        labels_combined += [0]*len(bacteria)
        labels_combined += [1]*len(venomous)
        labels_combined += [2]*len(non_venomous)
        labels_dataset += [0]*len(bacteria)
        labels_dataset += [1]*len(venomous)
        labels_dataset += [2]*len(non_venomous)
    else:
        non_venomous = [vecs[p] for p, l in enumerate(labels) if l == 2]
        X_combined += non_venomous
        labels_combined += [2]*len(non_venomous)
        labels_dataset += [2 + d]*len(non_venomous)

X = np.array(X_combined)
y = np.array(labels_combined)
X_trans_lda=lda.fit_transform(X,y)
ld = np.array(labels_dataset)
# Print the results
print('*************** LDA Summary ***************')
print('Classes: ', lda.classes_)
print('Priors: ', lda.priors_)
print('Explained variance ratio: ', lda.explained_variance_ratio_)
x = X_trans_lda
estimator = KMeans(n_clusters=3)
y_kmeans = estimator.fit_predict(x)
#empty dictionaries

clusters_centroids=dict()
clusters_radii= dict()

'''looping over clusters and calculate Euclidian distance of 
each point within that cluster from its centroid and 
pick the maximum which is the radius of that cluster'''

for cluster in list([0, 1, 2]):

    clusters_centroids[cluster]=list(zip(estimator.cluster_centers_[:, 0],estimator.cluster_centers_[:,1]))[cluster]
    clusters_radii[cluster] = max([np.linalg.norm(np.subtract(i,clusters_centroids[cluster])) for i in zip(x[y_kmeans == cluster, 0],x[y_kmeans == cluster, 1])])
#Visualising the clusters and cluster circles
all_ax[3].scatter(x[ld == 2, 0], x[ld == 2, 1], c = 'blue', label = 'Non-venomous 1')
all_ax[3].scatter(x[ld == 3, 0], x[ld == 3, 1], c = 'cyan', label = 'Non-venomous 2')
all_ax[3].scatter(x[ld == 4, 0], x[ld == 4, 1], c = 'magenta', label = 'Non-venomous 3')
all_ax[3].scatter(x[ld == 0, 0], x[ld == 0, 1], c = 'green', label = 'Bacteria')
all_ax[3].scatter(x[ld == 1, 0], x[ld == 1, 1], c = 'red', label = 'Venomous')
#Visualising cluster circles
art = mpatches.Circle(clusters_centroids[0],clusters_radii[0], edgecolor='b',fill=False)
all_ax[3].add_patch(art)
art = mpatches.Circle(clusters_centroids[2],clusters_radii[2], edgecolor='r',fill=False)
all_ax[3].add_patch(art)
art = mpatches.Circle(clusters_centroids[1],clusters_radii[1], edgecolor='g',fill=False)
all_ax[3].add_patch(art)

#Plotting the centroids of the clusters
all_ax[3].scatter(estimator.cluster_centers_[:, 0], estimator.cluster_centers_[:,1], s = 50, c = 'black', label = 'Centroids')
all_ax[3].grid(True)
all_ax[3].set_title("All datasets clustered")
all_ax[3].set_xlabel("LD1")
all_ax[3].set_ylabel("LD2")
all_ax[3].legend()


# In[ ]:




