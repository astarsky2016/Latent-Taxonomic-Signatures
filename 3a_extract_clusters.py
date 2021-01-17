import glob, pickle, sqlite3, random, sys, os, wget #presuming you have wget installed both via pip and in your PATH
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#This is where we create relaxed orphan dataset based on NCBI Clusters
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    return None
def cluster_information():
	cluster_tags = ["CHL", "CLSC", "CLSM", "CLSN", "CLSP", "CLSZ", "MTH", "PCLA", "PHA", "PLN", "PTZ"]
	cluster_info = defaultdict(set)
	protein_clusters = []
	#Download CLuster info from NCBI ftp
	for tag in cluster_tags:
		protein_cluster = "{}_proteins.txt".format(tag)
		url = "https://ftp.ncbi.nih.gov/genomes/CLUSTERS/" + protein_cluster
		wget.download(url)
		protein_clusters.append(protein_cluster)
	for f in protein_clusters:
		with open(f, "r") as handle:
			next(handle)
			for line in handle:
				data = line.split()
				cluster = data[0]
				acc = data[1].split(".")[0]
				taxid = data[-2]
				cluster_info[taxid].add(acc)
	print (len(cluster_info))
	sorted_clusters = sorted(cluster_info.items(), key=lambda tup: len(tup[1]), reverse = True)
	pickle.dump(cluster_info, open("cluster_info.p", "wb"))
	return cluster_info
if __name__ == '__main__':
	path = os.getcwd()
	new_path = path +* "/RELAXED_ORPHAN"
	try:
		os.mkdir(path)
	except OSError:
		print ("Creation of the directory {} failed".format(new_path))
	else:
		print ("Successfully created the directory {}".format(new_path))
	refseqDB = create_connection("nr_species.db")
	c = refseqDB.cursor()
	taxons = set(["{}".format(el) for el in pickle.load(open("docs_orphan.p", "rb"))])
	cnt = 0
	cluster_info = cluster_information()
	for taxon in cluster_info.keys():
		if taxon in taxons:
			if len(cluster_info[taxon]) >= 100:
				c.execute("SELECT acc, fasta_seq FROM proteome WHERE taxon_id=?", (taxon,)) 
				records = dict([(el[0], el[1]) for el in c.fetchall()])
				cluster_candidates = set(records.keys()).intersection(cluster_info[taxon]) #these are Clusters representatives stored in a local Sqlite DB
				orphan_candidates = set(records.keys()).difference(cluster_info[taxon]) #these are relaxed orphans
				if len(orphan_candidates) >= 100:
					orphan_selected = random.sample(list(orphan_candidates), 100)
					clusters_selected = random.sample(list(cluster_candidates), 100)
					orphan_fastas = [SeqRecord(Seq(records[ac]), id=ac, description="taxa {}".format(taxon)) for ac in orphan_selected]
					cluster_fastas = [SeqRecord(Seq(records[ac]), id=ac, description="taxa {}".format(taxon)) for ac in clusters_selected]
					with open("RELAXED_ORPHAN/{}_cluster.fasta".format(taxon), "w") as proteome_out:
						SeqIO.write(cluster_fastas, proteome_out, "fasta")	
					with open("RELAXED_ORPHAN/{}_orphan.fasta".format(taxon), "w") as proteome_out:
						SeqIO.write(orphan_fastas, proteome_out, "fasta")	
					cnt += 1
	print (cnt, "taxa processed!")
	refseqDB.close()






