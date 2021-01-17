import sqlite3, glob, logging, pickle, random, shelve, sys, os, re, shutil, gzip
from Bio import SeqIO
from ete3 import NCBITaxa
from collections import defaultdict
from gensim import corpora, models, similarities
from itertools import islice, tee
from gensim.test.utils import datapath
import urllib.request as request
from contextlib import closing
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
ncbi = NCBITaxa()
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
    return None
def create_table(c, create_table_sql):
    try:
        c.execute(create_table_sql)
    except:
        print("error")
def insertToTaxonomy(batch, c):
	task = 'INSERT INTO acc2taxa VALUES (?,?)'
    c.executemany(task1, batch)
def insertToDB(batch, c):
    task = 'INSERT INTO proteome VALUES (?,?,?,?)'
    c.executemany(task, batch)
def extractBatch(sequence_batch, taxa_cursor, refseq_cursor, counter):
	batch = []
	cnt = counter
	for seq_record in sequence_batch:
		accs = [seq_record.id.split(".")[0]]
		desc = seq_record.description
		the_rest = desc.split("]")
		if len(the_rest) >= 2:
			for part in the_rest[1:]:
				full = part.split()
				if full:
					multi = full[0].strip().replace('\x01', '')
					if multi:
						accs.append(multi.split(".")[0])
		for acc in accs:
			taxa_cursor.execute("SELECT taxa FROM acc2taxa WHERE acc=?", (acc,)) 
			rows = taxa_cursor.fetchall()
			for tx in rows:
				species_id = tx[0]
				sequence = str(seq_record.seq).lower()
				cnt += 1
				batch.append((cnt, species_id, acc, sequence))
	insertToDB(batch, refseq_cursor)
	return cnt
def makeDocs():
	#this is where we tokenize proteins into 3-peptides and store them into a database for later use
	sequencesDB = create_connection("nr_species.db")
	cur = sequencesDB.cursor()
	conn = create_connection("nr_text.db")
	sql_create_documents_table = """CREATE TABLE IF NOT EXISTS texts (
									taxon_id integer PRIMARY KEY,
	                                trigrams text NOT NULL
	                            );"""
	store = conn.cursor()
	create_table(store, sql_create_documents_table)
	conn.commit()
	cur.execute("select DISTINCT taxon_id FROM proteome")
	all_taxa = [el[0] for el in cur.fetchall()]
	cnt = 0
	docs = []
	for id in all_taxa:
		cur.execute("SELECT fasta_seq FROM proteome WHERE taxon_id=?", (id,)) 
		fastas = list(set([el[0] for el in cur.fetchall()]))
		lenght = len(fastas)
		sample_size = 0
		if lenght >= 100:
			if lenght > 5000:
				sample_size = 5000
			else:
				sample_size = lenght
		if sample_size:
			choice = random.sample(fastas, sample_size)
			text = ""
			for s in choice:
				trigrams = list(zip(*(islice(seq, index, None) for index, seq in enumerate(tee(s, 3)))))
				trigram_protein = "|".join(["".join(el) for el in trigrams]) + "|"
				text += trigram_protein
			store.execute("INSERT INTO texts VALUES (?, ?)", (id, text[0:-1]))
			conn.commit()
			docs.append(id)
			cnt += 1
			print(cnt, "processed")
	print ("texts extracted for {} species".format(len(docs)))
	store.execute("CREATE INDEX taxon_idx ON texts(taxon_id)")
	conn.commit()
	pickle.dump(docs, open("docs.p", "wb"))
	store.close()
	cur.close()
	return (docs)
class dictionary_stream(object):
	def __init__(self, db, docs):
		self.store = db
		self.docs = docs
	def __iter__(self):
		for id in self.docs:
			self.store.execute("SELECT trigrams FROM texts WHERE taxon_id=?", (id,))
			yield self.store.fetchall()[0][0].split("|")
class MyCorpus(object):
	def __init__(self, dict, docs, db):
		self.dictionary = dict
		self.stored_random = db
		self.documents = docs
	def __iter__(self):
		for t in self.documents:
			self.stored_random.execute("SELECT trigrams FROM texts WHERE taxon_id=?", (t,))
			text = self.stored_random.fetchall()[0][0].split("|")
			yield self.dictionary.doc2bow(text)
def populateDB_Taxonomy():
	#mapping taxonmy identifiers to protein acc
	acc2taxa = create_connection("acc2taxId.db")
	acc2taxa_cursor = acc2taxa.cursor()
	sql_create_acc2taxa_table = """CREATE TABLE IF NOT EXISTS acc2taxa (
	                                acc text PRIMARY KEY,
	                                taxa integer NOT NULL
	                            );"""
	create_table(acc2taxa_cursor, sql_create_acc2taxa_table)
	acc2taxa.commit()
	batch = []
	source = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
	with closing(request.urlopen(source)) as r:
		with open('prot.accession2taxid.gz', 'wb') as f:
			shutil.copyfileobj(r, f)
	f = gzip.open('prot.accession2taxid.gz', 'rb')
	file_content = f.read()
	f_line = 0
	for line in file_content:
		if f_line == 0:
			f_line = 1
			continue
		decoded = line.decode().strip()
		data = decoded.split()
		acc = data[0]
		taxId = data[2]
		batch.append((acc, taxId))
		if len(batch) == 50000:
			insertToDB(batch, acc2taxa_cursor)
			acc2taxa.commit()
			batch = []
	if batch:
		insertToDB(batch, acc2taxa_cursor)
		acc2taxa.commit()
		print ("Taxonomy mapped to accession")	
	print ("indexing")
	acc2taxa.execute("CREATE INDEX acc_idx ON acc2taxa(acc)")
	acc2taxa.commit()
	acc2taxa.close()

def populateDB_nr():
	#storing nr protein sequences locally
	acc2taxa = create_connection("acc2taxId.db")
	refseqDB = create_connection("nr_species.db")
	acc2taxa_cursor = acc2taxa.cursor()
	refseqDB_cursor = refseqDB.cursor()
	sql_create_proteins_table = """CREATE TABLE IF NOT EXISTS proteome (
									record_id integer PRIMARY KEY,
	                                taxon_id integer NOT NULL,
	                                acc text NOT NULL,
	                                fasta_seq text NOT NULL
	                            );"""
	create_table(refseqDB_cursor, sql_create_proteins_table)
	refseqDB.commit()
	batch = []
	source = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
	with closing(request.urlopen(source)) as r:
		with open('nr.gz', 'wb') as f_in:
			shutil.copyfileobj(r, f_in)
		f = gzip.open('nr.gz', 'rb')
		bindata = f.read()
		with open('nr', 'w') as f_out:
			f_out.write(bindata)
	for seq_record in SeqIO.parse("nr", "fasta"):
		batch.append(seq_record)
		if len(batch) == 50000:
			rec_counter = extractBatch(batch, acc2taxa_cursor, refseqDB_cursor, rec_counter)
			refseqDB.commit()
			print (rec_counter, "done!")
			batch = []
	if batch:
		rec_counter = extractBatch(batch, acc2taxa_cursor, refseqDB_cursor, rec_counter)
		refseqDB.commit()
		print (rec_counter, "final!")
		batch = []
	print ("nr sequences stored in a local sqlite database")
	print ("indexing")
	refseqDB_cursor.execute("CREATE INDEX taxon_idx ON proteome(taxon_id)")
	refseqDB.commit()
	refseqDB.close()
	acc2taxa.close()
if __name__ == '__main__':
	populateDB_Taxonomy()
	populateDB_nr()
	docs = makeDocs()
	texts_conn = create_connection("nr_text.db")
	texts = texts_conn.cursor()
	print ("Building LSA species model with all taxa >= 100 protein sequences in nr database")
	mem_friendly = dictionary_stream(texts, docs)
	dictionary = corpora.Dictionary(ngrams for ngrams in mem_friendly)
	dictionary.save("species_LSA.dict".format(d))
	dictionary = corpora.Dictionary.load("species_LSA.dict")
	corpus = MyCorpus(dictionary, docs, texts)
	corpora.MmCorpus.serialize("species_LSA.corp", corpus)
	corpus = corpora.MmCorpus("species_LSA.corp")
	tfidf = models.TfidfModel(corpus)
	tfidf.save("species_LSA.tfidf")
	tfidf = models.TfidfModel.load("species_LSA.tfidf")
	corpus_tfidf = tfidf[corpus]
	lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=400)
	lsi.save("species_LSA.lsi")
	index = similarities.MatrixSimilarity(lsi[corpus_tfidf])
	index.save("species_LSA.index")
	print ("LSA species model built with {} species!".format(len(docs)))
	texts.close()
	cursor_sequences.close()

