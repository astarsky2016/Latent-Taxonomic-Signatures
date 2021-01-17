import sqlite3, pickle, random
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def create_connection(db_file):
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return None
def over_1000(cur):
    taxons = defaultdict(int)
    cur.execute("SELECT taxon_id FROM proteome")
    for row in cur.fetchall():
        taxons[row[0]] += 1
    above = [el[0] for el in taxons.items() if el[1] >= 1000]
    print (len(above))
    pickle.dump(above, open("docs_above1000.p", "wb"))

if __name__ == '__main__':
    docs = pickle.load(open("docs_above1000.p", "rb"))
    conn = create_connection("nr_species.db")
    cur = conn.cursor()
    with open("500_proteome.fasta", "a") as proteome_out:
        cnt = 0
        for doc in docs:
            cur.execute("SELECT acc, fasta_seq FROM proteome WHERE taxon_id=?", (doc,)) 
            proteome = [el for el in cur.fetchall()]
            random_500 = random.sample(proteome, 500)
            proteome_fastas = [SeqRecord(Seq(rec[1]), id=rec[0], description="taxa {}".format(doc)) for rec in random_500]
            SeqIO.write(proteome_fastas, proteome_out, "fasta") 
            cnt += 1
    conn.close()
    cmd = "makeblastdb -in 500_proteome.fasta -dbtype prot"
    os.system(cmd)
    print ("Blast database created")
   
