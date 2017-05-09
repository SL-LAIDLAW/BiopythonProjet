#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Initialisation
#Packages utilisés
import os
import sys
import subprocess
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

from Bio.Align.Applications import ClustalOmegaCommandline
from contextlib import redirect_stdout
from Bio import AlignIO
from Bio import Phylo
import pylab
import re
import io


# Les arguments passés à notre script python à partir du .sh
default_search = str(sys.argv[1])
Entrez.email = str(sys.argv[2])
nom_utilisateur = str(sys.argv[3])


# Creeons notre dossier de travail
directory = "runFiles"
if not os.path.exists(directory):
    os.makedirs(directory)

# Creeons notre dossier de résultat
LaTexdir = "LaTexResults"
if not os.path.exists(LaTexdir):
    os.makedirs(LaTexdir)



## Obtention des données
if default_search == "True":
	# Si on choisi default
	Id_correspondant = "12583668"
	nom_gene_interet = "DEC2"


if default_search == "False":
	print(">> Vous avez choisi de rechercher une séquence")
	search_query = input("Votre Recherche?.. ")
	# Lancons la rêquete
	choixRetstart = 10
	choixRetmax = 10
	while True:
		requete_Brut = Entrez.esearch(db="Nucleotide", term=search_query, retmax=choixRetmax, retstart=choixRetstart)
		requete_Parsed = Entrez.read(requete_Brut)
		for myId in requete_Parsed["IdList"]:
			myId_name = Entrez.esummary(db="Nucleotide", id=myId, rettype="gb")
			for record in Entrez.parse(myId_name):
				Texte_a_ecrire = '{:>15}    {:.57}'.format(" " + str(myId), str(record['Title']))
				print("%s" % Texte_a_ecrire)

		Id_correspondant = input("\nQuel ID à choisir? (tappez 'plus' pour plus de valeurs, ou 'retour' pour changer l'objet de recherche).. ")
		if Id_correspondant == "plus":
			choixRetstart = int(choixRetstart) + 10
		elif Id_correspondant == "retour":
			search_query = input("Votre Recherche?.. ")
		else:
			nom_gene_interet = input("Nom du gène?.. ")
			break
	requete_Brut.close()


# Cherchons le séquence correspondant à notre Id
multifasta_var = []
print(">> Obtention des données...")
sequences_Brut = Entrez.efetch(db="Nucleotide", id=Id_correspondant, rettype="gb")
sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")

for chaqueSequence in sequences_Parsed:
	multifasta_var.append(chaqueSequence)
SeqIO.write(multifasta_var, str("/".join([directory,"gene_interet.fasta"])),"fasta")
print(">> Données enregistrés...")
sequences_Brut.close()

# Exportons séquence du gène d'intêret sous format fasta
for record in multifasta_var:
	Organisme_interet = record.description.replace("PREDICTED: ","")
	regex = r"^((?:\S+\s+){1}\S+).*"
	Organisme_interet = re.findall(regex, Organisme_interet)[0]
	Organisme_gene_interet = Organisme_interet
	print("Organisme du gene initiale est : " + Organisme_interet)



# Envoyons ces donnés pour un BLAST chez NCBI
print(">> BLASTons notre gène d'interet au NCBI...")
for record in SeqIO.parse(str("/".join([directory,"gene_interet.fasta"])),"fasta"):
	Blast_xml = NCBIWWW.qblast('blastn', 'nt', record.seq)
	blast_File = open(str("/".join([directory,"gene_interet_BLAST.xml"])), "w")
	blast_File.write(Blast_xml.read())
	Blast_xml.close()
	blast_File.close()


### Analyse des données
#Mise en place de la dictionaire {HSP bit-score : Identifiant}
HSPdict = {}
orgIDict = {}
OrganismeDict = {}
regexID = r"(?<=gi\|)(\d*)(?=\|)"
regexOrganisme = r"^((?:\S+\s+){1}\S+).*"

for record in NCBIXML.parse(open(str("/".join([directory,"gene_interet_BLAST.xml"])))) :
    if record.alignments :
        for align in record.alignments :
            for hsp in align.hsps:
            	Id_correspondant = re.findall(regexID, align.hit_id)[0]
            	newdictValue = {Id_correspondant : hsp.bits}
            	Organisme = align.hit_def.replace("PREDICTED: ","")
            	Organisme = re.findall(regexOrganisme, Organisme)[0]
            	orgvalue = {hsp.bits : Organisme}
            	orgID = { Organisme : Id_correspondant}
            	orgIDict.update(orgID)
            	OrganismeDict.update(orgvalue)
            	HSPdict.update(newdictValue)

### Triage des données
# Creer et trier une liste à partir des scores du dictionnaire
HPS_trie = []
HPS_brut = []
HPS_organismes = []
for value in HSPdict.values():
	HPS_brut.append(value)

HPS_brut = sorted(HPS_brut, reverse=True)


# Creeons une liste identique sans duplicats ou résultats de l'espèce d'interet
for value in HPS_brut:
	if str(OrganismeDict[value]) in Organisme_interet:
		pass
	else:
		if str(OrganismeDict[value]) in HPS_organismes:
			pass

		else:
			HPS_organismes.append(str(OrganismeDict[value]))
			HPS_trie.append(value)

top_HPS = HPS_trie[:5]

# Inversons le dictionnaire pour pouvoir utiliser le [key] avec les « values »
HSPdict = {y:x for x,y in HSPdict.items()}

# Ecrivons les titres dans le fichier résultats
results_File = open(str("/".join([directory,"results.txt"])), "w")
Texte_a_ecrire = '{:^9}     {:^25} {:^21}'.format("Bit-Score", "Espèce","Identifiant")
results_File.write("%s\n" % Texte_a_ecrire)
Texte_a_ecrire = '{:_<9}     {:_<25} {:_<21}'.format("", "","")
results_File.write("%s\n" % Texte_a_ecrire)
results_File.write("\n>> Résultats pour " + Organisme_interet + " :\n")

# Ecrivons ces premiers résultats dans le fichier résultats
TableaugeneInteret = ""
for score in top_HPS:
	Texte_a_ecrire = '{:>9}      {:25} {:21}'.format(str(float(score)), str(OrganismeDict[score]),str(HSPdict[score]))
	results_File.write("%s\n" % Texte_a_ecrire)
	TableaugeneInteret = TableaugeneInteret + (str(float(score)) + " & " + OrganismeDict[score] + " & " + HSPdict[score] + "\\\\\n")
results_File.write("" + "\n")


# Mettons tous les identifiants ensemble pour que le NCBI les comprends 
for value in HPS_trie[:4]:
	Id_correspondant = ",".join([str(HSPdict[value]),Id_correspondant])

# Exportons séquence du gène d'intêret sous format fasta
sequences_Brut = Entrez.efetch(db="Nucleotide", id=Id_correspondant, rettype="gb")
print(">> Données enregistrés...")
sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")

multifasta_var = []
for chaqueSequence in sequences_Parsed:
	multifasta_var.append(chaqueSequence)

# Creeons le base de donnés qui servira ensuite pour le BLAST en local
SeqIO.write(multifasta_var, str("/".join([directory,"base_de_donnes.fasta"])),"fasta")
sequences_Brut.close()

cmd = "makeblastdb -in " + str("/".join([directory,"base_de_donnes.fasta"])) + " -parse_seqids -dbtype nucl"
subprocess.check_output(
	cmd,
	stderr=subprocess.STDOUT,
	shell=True)


# Separons les identifiants pour que mon script les comprends 
top_HPS_ids = Id_correspondant.split(",")

myIds = []
print("")
for value in HPS_trie[0:5]:
	# print(str(value) + "\t" + str(HSPdict[value]) + "\t" +  str(OrganismeDict[value]))
	myIds.append(HSPdict[value])

# On inverse de nouveau, c'est à dire qu'on revient aux clés/valeurs qu'on avait avant
HSPdict = {y:x for x,y in HSPdict.items()}


# Commoncons le compteur pour pouvoir suivre le nombre de boucles qui ont lieu
compteur = 0
Tableaugene = {}
ListeOrganismesInteret = {}

# Boucle pour blaster, et ecrire les résultats pour chacun de nos résultats du blast précedent
for myId in myIds:
	# Declarons variables, et disons ce qu'on fait
	compteur +=1
	Organisme_interet = str(OrganismeDict[HSPdict[myId]])
	newdictValue = {compteur : Organisme_interet}
	ListeOrganismesInteret.update(newdictValue)
	message = "\n> Début du boucle nb." + str(compteur) + " sur l'espèce :"
	message = str("echo \"$(tput setaf 6)" + message + "$(tput sgr 0)\"")
	os.system(message)
	
	print(str(myId) + " : " + Organisme_interet)
	
	# Cherchons le séquence correspondant à notre Id
	sequences_Brut = Entrez.efetch(db="Nucleotide", id=myId, rettype="gb")
	print(">> Données enregistrés...")
	sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")


	# Exportons séquence du gène d'intêret sous format fasta
	multifasta_var = []
	exportNameSeq = str(directory) + "/" + str(compteur) + "_" + Organisme_interet + ".fasta"

	for chaqueSequence in sequences_Parsed:
		multifasta_var.append(chaqueSequence)
	SeqIO.write(multifasta_var, exportNameSeq,"fasta")
	sequences_Brut.close()


	### Analyse des données
	print(">> Analyse BLAST des séquences...")
	
	# create unique export name
	exportNameBLAST = str(directory) + "/" + str(compteur) + "_" + Organisme_interet + "_BLAST.xml"

	# BLASTons en local
	Blast_xml = NcbiblastnCommandline(query=exportNameSeq,db=str("/".join([directory,"base_de_donnes.fasta"])), outfmt=5, out=exportNameBLAST)
	os.system(str(Blast_xml))


	#Mise en place des dictionaires
	microHSPdict = {}
	for record in NCBIXML.parse(open(exportNameBLAST)) :
		if record.alignments :
			for align in record.alignments :
				for hsp in align.hsps:
					newdictValue = {align.hit_id : hsp.bits}
					Organisme = align.hit_def.replace("PREDICTED: ","")
					Organisme = re.findall(regexOrganisme, Organisme)[0]
					orgvalue = {hsp.bits : Organisme}
					OrganismeDict.update(orgvalue)
					microHSPdict.update(newdictValue)	            		


	### Triage des données
	HPS_trie = []
	HPS_brut = []
	HPS_organismes = []
	for value in microHSPdict.values():
		HPS_brut.append(value)

	HPS_brut = sorted(HPS_brut, reverse=True)
	for value in HPS_brut:
		if str(OrganismeDict[value]) in Organisme_interet:
			pass
		else:
			HPS_trie.append(value)

	# Echanger keys et values afin d'extraire l'identifiant à partir du HSP bit-score
	microHSPdict = {y:x for x,y in microHSPdict.items()}

	# Obtenons l'identifiant à partir de HSP bit-score trouvé ci-haut
	HSPresults = []
	results_File = open(str("/".join([directory,"results.txt"])), "a")
	results_File.write(">> Résultats pour " + Organisme_interet + " :\n")

	accencion2Id = {}
	Tableaugene_temp = ""
	regexID = r"(?<=[a-z]{3}\|)(.*)(?=\|)"
	for score in HPS_trie:
		microID = orgIDict[OrganismeDict[score]]
		try:
			accencionID = re.findall(regexID, microHSPdict[score])[0]
		except:
			accencionID = microHSPdict[score]
		newdictValue = {accencionID : microID}
		accencion2Id.update(newdictValue)
		Texte_a_ecrire = '{:>9}      {:25} {:21}'.format(" " + str(float(score)), str(OrganismeDict[score]),str(microID))
		results_File.write("%s\n" % Texte_a_ecrire)
		Tableaugene_temp = str(Tableaugene_temp) + (str(float(score)) + " & " + OrganismeDict[score] + " & " + microID + "\\\\\n")
	results_File.write("" + "\n")

	Tableaugene.update({compteur : Tableaugene_temp}) 
	print("> fin du boucle nb." + str(compteur))


print("\n>> Resultats Enregistrés sous \"resultats.txt\"")


### Creation de l'arbre phylogenique
# Creeons le fichier newark pour l'arbre
multifasta = "base_de_donnes.fasta"
guidetree_accession = str("/".join([directory,"guidetree_accession.txt"]))
clustalomega_cline = ClustalOmegaCommandline(infile=multifasta, outfile="base_de_donnes_alignes.fasta", guidetree_out="guidetree_accession.txt",verbose=True, auto=True, force=True)
cmd = "cd runFiles\n" + str(clustalomega_cline)
cmd = str(cmd)
os.system(cmd)


# Remplacons le numero ascencion avec le nom de l'espèce
guidetree_binomiale = str("/".join([directory,"guidetree_binomiale.txt"]))
with open(guidetree_accession, 'r') as file :
  filedata = file.read()

orgIDict = {y:x for x,y in orgIDict.items()}
for key in accencion2Id.keys():
	especeName = orgIDict[accencion2Id[key]]
	especeName = especeName.replace(" ", "_")
	filedata = filedata.replace(key, especeName)

with open(guidetree_binomiale, 'w+') as file :
  file.write(filedata)


# Creer à partir de guidetree binomiale l'arbre en ASCII
arbre_out = io.StringIO()
arbre_in = Phylo.read(guidetree_binomiale, "newick")

Phylo.draw_ascii(arbre_in,file=arbre_out,column_width=90)
arbre_ascii = str(arbre_out.getvalue())
arbre_out.close()





### Creons le fichier LaTex pour afficher les résultats
with open("Resultat_Template.tex", 'r') as file :
  resultat_latex = file.read()

# Creer le header
resultat_latex = resultat_latex.replace("[NomUtilisateur*]", nom_utilisateur)
resultat_latex = resultat_latex.replace("[emailUtilisateur*]", Entrez.email)
resultat_latex = resultat_latex.replace("[geneInteret*]", nom_gene_interet)

# Creer les Tableaux dans le fichier LaTex
resultat_latex = resultat_latex.replace("[TableaugeneInteret*]", str(TableaugeneInteret))
resultat_latex = resultat_latex.replace("[especegeneInteret*]", Organisme_gene_interet)

resultat_latex = resultat_latex.replace("[Tableaugene1*]",str(Tableaugene[1]))
resultat_latex = resultat_latex.replace("[especegene1*]",ListeOrganismesInteret[1])

resultat_latex = resultat_latex.replace("[Tableaugene2*]",str(Tableaugene[2]))
resultat_latex = resultat_latex.replace("[especegene2*]",ListeOrganismesInteret[2])

resultat_latex = resultat_latex.replace("[Tableaugene3*]",str(Tableaugene[3]))
resultat_latex = resultat_latex.replace("[especegene3*]",ListeOrganismesInteret[3])

resultat_latex = resultat_latex.replace("[Tableaugene4*]",str(Tableaugene[4]))
resultat_latex = resultat_latex.replace("[especegene4*]",ListeOrganismesInteret[4])

resultat_latex = resultat_latex.replace("[monarbrephylo*]", arbre_ascii)
with open("LaTexResults/Resultats.tex", 'w+') as file :
  file.write(resultat_latex)

os.system("cd LaTexResults\npdflatex --shell-escape --file-line-error --synctex=1 Resultats.tex")
