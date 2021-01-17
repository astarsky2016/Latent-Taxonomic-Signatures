from Bio import SeqIO
from collections import defaultdict
from itertools import islice, tee
import random, pickle
from multiprocessing import Pool
def make_bates_mix(p):
    trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(p, 3)))))
    text = ["".join(el) for el in trigrams]
    shuffled_aa = [e for e in "".join(text)]
    random.shuffle(shuffled_aa)
    text_shuffle = ["".join(shuffled_aa[x:x+3]) for x in range(0,len(shuffled_aa),3)]
    alphabet = list(set(shuffled_aa))
    random_str = "".join([random.choice(alphabet) for _ in range(len(shuffled_aa))])
    text_random = [random_str[x:x+3] for x in range(0,len(random_str),3)]
    return {"real":text, "shuffle":text_shuffle, "random":text_random}

if __name__ == '__main__':
	frequencies = {"real":defaultdict(int), "shuffle":defaultdict(int), "random":defaultdict(int)}
	batch = []
	p = Pool(30)
	for rec in SeqIO.parse(open("500_proteome.fasta"), "fasta"):
		seq = str(rec.seq)
		batch.append(seq)
		if len(batch) == 100000:
			results = p.map(make_bates_mix, batch)
			for n_grams in results:
				for ds in n_grams.keys():
					for aa in n_grams[ds]:
						frequencies[ds][aa] += 1
			batch = []
			print ("done!")
	if batch:
		p = Pool(30)
		results = p.map(make_bates_mix, batch)
		for n_grams in results:
			for ds in n_grams.keys():
				for aa in n_grams[ds]:
					frequencies[ds][aa] += 1
		n_grams = make_bates_mix(batch)
		for ds in n_grams.keys():
			for aa in n_grams[ds]:
				frequencies[ds][aa] += 1
	frequencies_pickle = {"real":dict(frequencies["real"]), "shuffle":dict(frequencies["shuffle"]), "random":dict(frequencies["random"])}
	pickle.dump(frequencies_pickle, open("3_gram_frequencies.p", "wb"))
	p.close()
	p.join()
	n_grams = pickle.load(open("3_gram_frequencies.p", "rb"))
	for series in n_grams.keys():
		total = sum(n_grams[series].values())
		for k in n_grams[series].keys():
			n_grams[series][k] = n_grams[series][k]/total
	pickle.dump(n_grams, open("3_gram_frequencies_percent.p", "wb"))


