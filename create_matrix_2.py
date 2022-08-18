"""This is the code that creates LSA model based on collected species proteome datasets. For this purpose, these proteomes have to be collected (previous program does that) and Gensim - 
an open-source library for unsupervised topic modeling, document indexing, retrieval by similarity, and other natural language processing functionalities, using modern statistical machine learning; 
has to be installed. 
"""
import logging, pickle, sys, os, itertools, plyvel
from gensim import corpora, models, similarities
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
class MyCorpus(object):
	def __init__(self, dictionary, docs, db):
		self.dictionary = dictionary
		self.documents = docs
		self.cursor = db
	def __iter__(self):
		for species in self.documents:
			trigrams = self.cursor.get(bytes(str(species), "utf8"))
			texts = trigrams.decode("utf8").split("|")
			yield self.dictionary.doc2bow(texts)
class MyDictionary(object):
	def __init__(self, docs, db):
		self.documents = docs
		self.cursor = db
	def __iter__(self):
		for species in self.documents:
			trigrams = self.cursor.get(bytes(str(species), "utf8"))
			texts = trigrams.decode("utf8").split("|")
			yield texts
def matrix():
	db = plyvel.DB('./tax2text_levDB')
	docs = pickle.load(open("species.documents", "rb"))
	print ("starting")
	print (len(docs), "species collected in total")
	dictionary_mem = MyDictionary(docs, db)
	dictionary = corpora.Dictionary(s for s in dictionary_mem)
	dictionary.save("species.dictionary")
	dictionary = corpora.Dictionary.load("species.dictionary")
	print ("dictionary created")
	corpus = MyCorpus(dictionary, docs, db)
	corpora.MmCorpus.serialize('species.corpus', corpus)
	corpus = corpora.MmCorpus("species.corpus")
	print ("corpus created")
	tfidf = models.TfidfModel(corpus)
	corpus_tfidf = tfidf[corpus]
	tfidf.save("species.tfidf")
	print ("tfidf model created")
	lsi_model = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=400)
	lsi_model.save("species.lsi")
	corpus_lsi = lsi_model[corpus_tfidf]
	index = similarities.MatrixSimilarity(corpus_lsi)
	index.save("species.index")
	print ("finished matrix creation")
	db.close()
if __name__ == '__main__':
	matrix()
