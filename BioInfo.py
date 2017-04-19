#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Initialisation
#Packages utilisés
import sys
import re
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# Les arguments passés à notre script python à partir du .sh
default_search = str(sys.argv[1])
format_arg = str(sys.argv[2])
Entrez.email = str(sys.argv[3])
skip_humans = str(sys.argv[4])
skip_predictions = str(sys.argv[5])




## Obtention des données
if default_search == "True":
	# Si on choisi default
	Id_correspondants = "12583668,110559298,12583670,157278206,1046887222"


if default_search == "False":
	search_query = raw_input("Votre Recherche?.. ")
	# Lancons la rêquete
	requete_Brut = Entrez.esearch(db=format_arg, term=search_query, retmax=10)
	requete_Parsed = Entrez.read(requete_Brut)
	Id_list = requete_Parsed("IdList")

	print(requete_Parsed["IdList"])
	whichID = raw_input("Quel ID à choisir?.. ")
	Id_correspondants = whichID
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
Blast_xml.close()
blast_File.close()


#Mise en place de la dictionaire {HSP bit-score : Identifiant}
HSPdict = {}
for record in NCBIXML.parse(open("BLAST.xml")) :
    if record.alignments :
        for align in record.alignments :
            for hsp in align.hsps :
                newdictValue = {align.hit_def : hsp.bits}
                HSPdict.update(newdictValue)




### Triage des données
# Suprimons pairs de donnés qui concerne les Homo sapiens si demandé
if skip_humans == "True":
	regex = r"^.*\b(homo sapiens)\b.*$"
	for key in HSPdict.keys():
		HomoSapien_Present = bool(re.search(regex, key.lower()))
		if HomoSapien_Present is True:
			del HSPdict[key]


# Suprimons pairs de donnés qui sont des prédictions si demandé
if skip_predictions == "True":
	regex = r"^.*\b(predicted)\b.*$"
	for key in HSPdict.keys():
		predicted_Present = bool(re.search(regex, key.lower()))
		if predicted_Present is True:
			del HSPdict[key]		


# Obtenons le top 10 HSP bit-score de notre recherche
HPS_list = []
for value in HSPdict.values():
	HPS_list.append(value)
HPS_list = sorted(HPS_list, reverse=True)
top_HPS = HPS_list[:10]


# Echanger keys et values afin d'extraire l'identifiant à partir du HSP bit-score
HSPdict = {y:x for x,y in HSPdict.items()}


#Obtenons l'identifiant à partir de HSP bit-score trouvé ci-haut
for score in top_HPS:
	print score , "\t", HSPdict[score]
