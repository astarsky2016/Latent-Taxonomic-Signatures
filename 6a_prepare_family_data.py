import glob, pickle, sys, random, wget
from collections import defaultdict
from ete3 import NCBITaxa
from Bio import SeqIO
from itertools import islice, tee
from gensim import models, similarities, corpora
ncbi = NCBITaxa()
def build_pfam():
	#download Pfam proteomes data and pfamseq FASTA sequences
	path = os.getcwd()
	new_path = path + "/PFAM"
	os.mkdir(new_path)
	ftp = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/*.tsv.gz"
	wget.download(ftp, out=new_path)
	zipped = glob.glob(new_path + "/*.gz")
	for compressed in zipped:
		f = gzip.open(compressed, 'rb')
		bindata = f.read()
		f_name = compressed.replace(".gz", "")
		with open(f_name, 'w') as f_out:
			f_out.write(bindata)
	ftp_seq = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamseq.gz"
	wget.download(ftp, out=new_path)
	f = gzip.open(new_path + "/pfamseq.gz", 'rb')
	bindata = f.read()
	with open(new_path + "/pfamseq", 'w') as f_out:
		f_out.write(bindata)
	proteomes = glob.glob(new_path + "/*.tsv")
	taxRestrData = defaultdict(lambda: defaultdict(list))
	families2proteomes = defaultdict(lambda: defaultdict(list))
	pf2seq = {}
	print (len(proteomes))
	cnt = 0
	for p in proteomes:
		with open(p, "r") as proteome:
			taxon = int(os.path.basename(p).replace(".tsv", ""))
			for i in range(3):
				next(proteome)
			for line in proteome:
				data = line.split()
				acc = data[0]
				pfam = data[5]
				taxRestrData[pfam][taxon].append(acc)
				families2proteomes[taxon][pfam].append(acc)
				pf2seq[acc] = (pfam, taxon)
		cnt += 1
		print (cnt, "done!")
	bkp = {}
	bkp_fam = {}
	for pf in taxRestrData.keys():
		bkp[pf] = dict(taxRestrData[pf])
	for tx in families2proteomes.keys():
		bkp_fam[tx] = dict(families2proteomes[tx])
	pickle.dump(bkp_fam, open("taxa2pfam.p", "wb"))
	pickle.dump(bkp, open("pfam_taxa.p", "wb"))
	pickle.dump(pf2seq, open("pfam2seq.p", "wb"))
def restricted():
	#here we collect protein family sequences grouped under taxonomic ranks of genus and species 
	data = pickle.load(open("pfam_taxa.p", "rb"))
	matrix_data = set(pickle.load(open("docs.p", "rb")))
	tax_restricted = {"species":defaultdict(lambda: defaultdict(set)), "genus":defaultdict(lambda: defaultdict(set))}
	N = 100 #we require at least 100 sequences per family per taxonomic rank
	restricted_data = {}
	for pfam in data.keys():
		genus_family_map = {}
		matching_taxa = [species for species in data[pfam].keys() if species in matrix_data]
		for species in matching_taxa:
			if len(data[pfam][species]) >= N:
				try:
					lineage = ncbi.get_lineage(species)
					ranks = ncbi.get_rank(lineage)
					inverted = dict([(el[1], el[0]) for el in ranks.items()])
					g = inverted.get("genus", None)
					if ranks[species] == "species":
						tax_restricted["species"][pfam][species].update(set(data[pfam][species]))
						continue
					else:
						if g:
							tax_restricted["genus"][pfam][g].update(set(data[pfam][species]))
				except:
					continue
	pickle_restricted = {}
	for c_rank in tax_restricted.keys():
		tmp = {}
		for pf in tax_restricted[c_rank].keys():
			taxa = list(tax_restricted[c_rank][pf].keys())
			if len(taxa) >= 2:
				tmp[pf] = dict(tax_restricted[c_rank][pf])
		if tmp:
			pickle_restricted[c_rank] = tmp
	pickle.dump(pickle_restricted, open("pfam_restricted_data.p", "wb"))
def extract_mixedSequences():
	#here we mix all taxa sequences into a joined pool from we we construct mixed origin dataset
	experiment = {"species":defaultdict(dict), "genus":defaultdict(dict)}
	#data = pickle.load(open("pfam_restricted_data.p", "rb"))
	compare = pickle.load(open("pfam_hypothesis_5.p", "rb"))
	pf2seq = pickle.load(open("pfam2seq.p", "rb"))
	lsi = models.LsiModel.load("species_LSA.lsi")
	tfidf = models.TfidfModel.load('species_LSA.tfidf')
	dictionary = corpora.Dictionary.load('species_LSA.dict')
	pf_fasta = defaultdict(list)
	pfam = set()
	for rnk in compare.keys():
		for pf in compare[rnk].keys():
			pfam.add(pf)
	with open("pfamseq", "rU") as handle:
		for record in SeqIO.parse(handle, "fasta") :
			acc =record.id.split(".")[0]
			acc_grouping = pf2seq.get(acc, None)
			if acc_grouping:
				pf_fasta[acc_grouping[0]].append((str(record.seq).lower(), acc_grouping[1]))
	print ("fasta collected")
	for rnk in compare.keys():
		for pf in compare[rnk].keys():
			fasta_data = pf_fasta[pf]
			cmp_taxa = set(compare[rnk][pf].keys())
			groupings = set([el[1] for el in fasta_data])
			fastas = [el[0] for el in fasta_data]
			size = int(len(fastas)/len(cmp_taxa))
			if size > 100:
				size = 100
			for tx in cmp_taxa:
				texts = []
				for s in random.sample(fastas, size):
					trigrams = ["".join(el) for el in list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(s, 3)))))]
					texts += trigrams
				experiment[rnk][pf][tx] = lsi[tfidf[dictionary.doc2bow(texts)]]
	pickable = {}
	for rank in experiment.keys():
		pfams = dict(experiment[rank])
		pickable[rank] = pfams
	pickle.dump(pickable, open("pfam_hypothesis_5_mixed2.p", "wb"))
def extractSequences():
	#we prepare both species and genus sequence datasets
	#genus dataset contains sequences which were not included in species dataset - it is independent
	experiment = {"species":defaultdict(lambda: defaultdict(list)), "genus":defaultdict(lambda: defaultdict(list))}
	control_counter = {"species":defaultdict(lambda: defaultdict(int)), "genus":defaultdict(lambda: defaultdict(int))}
	data = pickle.load(open("pfam_restricted_data.p", "rb"))
	accepted = {}
	for c_rank in data.keys():
		for pf in data[c_rank].keys():
			for tx in data[c_rank][pf].keys():
				accs = list(data[c_rank][pf][tx])
				for acc in accs:
					accepted[acc] = {"pfam":pf, "taxon":tx, "rank":c_rank}
	print (len(accepted))
	cnt = 0
	with open("pfamseq", "rU") as handle:
		for record in SeqIO.parse(handle, "fasta") :
			acc =record.id.split(".")[0]
			is_accepted = accepted.get(acc, None)
			if is_accepted:
				threshold = control_counter[is_accepted["rank"]][is_accepted["pfam"]][is_accepted["taxon"]]
				if threshold < 100:
					s = str(record.seq).lower()
					trigrams = ["".join(el) for el in list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(s, 3)))))]
					experiment[is_accepted["rank"]][is_accepted["pfam"]][is_accepted["taxon"]] += trigrams
					control_counter[is_accepted["rank"]][is_accepted["pfam"]][is_accepted["taxon"]] += 1
	pickable = {}
	for rank in experiment.keys():
		pfams = {}
		for pf_family in experiment[rank].keys():
			taxons = {}
			for tx in experiment[rank][pf_family].keys():
				if control_counter[rank][pf_family][tx] >= 90:
					taxons[tx] = experiment[rank][pf_family][tx]
				else:
					print ("low number - skipping!")
			if len(taxons) >= 2:
				pfams[pf_family] = taxons
		if pfams:
			pickable[rank] = pfams
	pickle.dump(pickable, open("pfam_hypothesis_5.p", "wb"))
def getSequences():
	sequences = {}
	print ("getting")
	cnt = 0
	with open("pfamseq", "rU") as handle:
		for record in SeqIO.parse(handle, "fasta") :
			acc =record.id.split(".")[0]
			sequences[acc] = str(record.seq).lower()
			cnt += 1
	pickle.dump(sequences, open("pfam_sequences.p", "wb"))
if __name__ == '__main__':
	build_pfam()
	restricted()
	extractSequences()
	extract_mixedSequences()
	getSequences()

