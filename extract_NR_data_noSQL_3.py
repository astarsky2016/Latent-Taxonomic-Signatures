"""This is the code that creates taxonomy benchmarking datasets for BLAST and LSA testing. Here, candidate taxa having in excess of 500 proteins are being selected,
the 500 proteins are being randomly excluded and datasets are being stored in a LevelDB instance.
"""
import glob, logging, pickle, random, shelve, sys, os, itertools, plyvel, random
from itertools import islice, tee
from ete3 import NCBITaxa
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from multiprocessing import Pool
def extractTrigrams(data):
	species = data[0]
	sequence_groups = data[1]
	results = {}
	for subset in sequence_groups.keys():
		texts = []
		for p in sequence_groups[subset]:
			trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
			texts += ["".join(el) for el in trigrams]
		results[subset] = (bytes(str(species), 'utf8'), bytes("|".join(texts), 'utf8'))
	return results
def prepareSpeciesTexts():
	db_text = plyvel.DB('./tax2text_500_levDB', create_if_missing=True)
	sub_dbs = {"main":db_text.prefixed_db(b'main-').write_batch(), "sub500":db_text.prefixed_db(b'sub500-').write_batch(), "sub100":db_text.prefixed_db(b'sub100-').write_batch(), 
	"sub50":db_text.prefixed_db(b'sub50-').write_batch(), "in500":db_text.prefixed_db(b'in500-').write_batch(), "in100":db_text.prefixed_db(b'in100-').write_batch(), "in50":db_text.prefixed_db(b'in50-').write_batch()}
	tags = ["main", "sub500", "sub100", "sub50", "in500", "in100", "in50"]
	blast_DB = open("main_500.fasta", "a")
	in_500 = open("in_500_query.fasta", "a")
	out_500 = open("out_500_query.fasta", "a")
	upper_limit = {2157: 15000, 10239: 3000, 2: 15000, 2759: 100000}
	lower_limit = {2157: 300, 10239: 10, 2: 300, 2759: 1000}
	db = plyvel.DB('./tax2prot_levDB')
	candidates = pickle.load(open("taxa_candidates_kingdoms_level.p", "rb"))
	docs_500 = []
	docs = []
	seq_counter = 0
	query_counter = 0
	for superkingdom in candidates.keys():
		print (superkingdom)
		batch = []
		for tx in candidates[superkingdom]:
			proteome = set()
			upper = upper_limit[superkingdom]
			lower = lower_limit[superkingdom]
			proteome = set([protein.decode("utf8") for protein in db.iterator(prefix=bytes("{}_".format(tx), 'utf8'), include_key=False)])
			size = len(proteome)
			if size >= lower:
				docs.append(tx)
				selection = list(proteome)
				subsample500 = []
				subsample100 = []
				subsample50 = []
				if size > upper:
					selection = random.sample(selection, upper)
				resized = len(selection)
				if resized >= 1000 and resized - 500 >= lower:
					docs_500.append(tx)
					random.shuffle(selection)
					subsample500 = selection[0:500]
					subsample100 = subsample500[0:100]
					subsample50 = subsample500[0:50]
					selection = selection[500:]
					random.shuffle(selection)
					in500 = selection[0:500]
					in100 = in500[0:100]
					in50 = in500[0:50]
					out_accs = ["{}_{}".format(tx, p) for p in range(query_counter, query_counter + len(subsample500))]
					in_accs = ["{}_{}".format(tx, p) for p in range(query_counter, query_counter + len(in500))]
					sequences_in = [SeqRecord(Seq(tup[0]), id=tup[1], description="BLAST in query") for tup in zip(in500, in_accs)]
					sequences_out = [SeqRecord(Seq(tup[0]), id=tup[1], description="BLAST out query") for tup in zip(subsample500, out_accs)]
					SeqIO.write(sequences_in, in_500, "fasta")
					SeqIO.write(sequences_out, out_500, "fasta")
					query_counter += 500
				accs = ["{}_{}".format(tx, p) for p in range(seq_counter, seq_counter + len(selection))]
				sequences = [SeqRecord(Seq(tup[0]), id=tup[1], description="BLAST DB") for tup in zip(selection, accs)]
				proteins_selected = {"main": selection, "sub500":subsample500, "sub100":subsample100, "sub50":subsample50, "in500":in500, "in100":in100, "in50":in50}
				batch.append((tx, proteins_selected))
				seq_counter += len(selection)
				SeqIO.write(sequences, blast_DB, "fasta")
			if len(batch) == 100:
				with Pool(15) as pool:
					results = pool.map(extractTrigrams, batch)
					for r in results:
						for subset in r.keys():
							if r[subset][1]:
								sub_dbs[subset].put(r[subset][0], r[subset][1])
				for tag in tags:
					sub_dbs[tag].write()
					sub_dbs[tag].clear()
				print (len(docs), len(docs_500))
				batch = []
		if batch:
			with Pool(15) as pool:
				results = pool.map(extractTrigrams, batch)
				for r in results:
					for subset in r.keys():
						if r[subset][1]:
							sub_dbs[subset].put(r[subset][0], r[subset][1])
			for tag in tags:
				sub_dbs[tag].write()
				sub_dbs[tag].clear()
	pickle.dump((docs, docs_500), open("species_removed500_control.documents", "wb"))
	db.close()
	db_text.close()
	blast_DB.close()
	in_500.close()
	out_500.close()
	print ("texts ready", len(docs), "species")
if __name__ == '__main__':
	prepareSpeciesTexts()

