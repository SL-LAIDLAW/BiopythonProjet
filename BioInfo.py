#!/usr/bin/env python
# -*- coding: utf-8 -*-
	## essentiel afin de ne pas encontrer des problèmes liés aux accents avec python < 3


### Initialisation
#Packages utilisés
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# Les arguments passés à notre script python à partir du .sh
search_query = str(sys.argv[1])
format_arg = str(sys.argv[2])
Entrez.email = str(sys.argv[3])





### Obtention des données
# Lancons la rêquete
requete_Brut = Entrez.esearch(db=format_arg, term=search_query, retmax=5)
requete_Parsed = Entrez.read(requete_Brut)
Id_correspondants = requete_Parsed["IdList"]
requete_Brut.close()


# Obtenons nos sequences à partir des Id qu'on vient de trouver
sequences_Brut = Entrez.efetch(db=format_arg, id=Id_correspondants, rettype="gb")
sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")


# Exportons ces sequences sous format fasta
multifasta_var = []
for chaqueSequence in sequences_Parsed:
	multifasta_var.append(chaqueSequence)
SeqIO.write(multifasta_var, "sequences.fasta","fasta")
sequences_Brut.close()




### Analyse des données
# Lisons le contenu de notre fichier fasta
with open('sequences.fasta', 'r') as sequence_fasta:
    multifasta_var=sequence_fasta.read()


# Envoyons ces donnés pour un BLAST chez NCBI
Blast_xml = NCBIWWW.qblast('blastn', 'nt', multifasta_var)
blast_File = open("BLAST.xml", "w")
blast_File.write(Blast_xml.read())
sequence_fasta.close()
blast_File.close()
Blast_xml.close()


# Analyse du fichier BLAST
with open('BLAST.xml', 'r') as blast_File:
    Blast_xml=blast_File.read()














