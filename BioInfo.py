#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Initialisation
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from contextlib import redirect_stdout
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import subprocess
import pylab
import os
import sys
import re
import io

# Les arguments passés à notre script python à partir du .sh
default_search = str(sys.argv[1])
Entrez.email = str(sys.argv[2])
nom_utilisateur = str(sys.argv[3])
skip_ncbi = str(sys.argv[4])

# Creeons notre dossier de travail
directory = "runFiles"
if not os.path.exists(directory):
    os.makedirs(directory)

# Creeons notre dossier de résultat
LaTexdir = "LaTexResults"
if not os.path.exists(LaTexdir):
    os.makedirs(LaTexdir)



### Obtention des données :
# Si on choisi de poursuivre le gène du projet
if default_search == "True":
	Identifiant = "12583668"
	nom_gène_intérêt = "DEC2"

# Si on veut choisir notre gène à partir de notre moteur de recherche
if default_search == "False":
	print(">> Vous avez choisi de rechercher une séquence")
	objet_récherche = input("Votre Recherche?.. ")
	numéro_RetStart = 10
	numéro_RetMax = 10
	while True:
		# Lancons la rêquete
		requête_Brut = Entrez.esearch(db="Nucleotide", term=objet_récherche, retmax=numéro_RetMax, retstart=numéro_RetStart)
		requête_Parsed = Entrez.read(requête_Brut)

		# Montrer le nom et identifiant de chacun des résultats de recherche
		for mon_Identifiant in requête_Parsed["IdList"]:
			mon_relevé = Entrez.esummary(db="Nucleotide", id=mon_Identifiant, rettype="gb")
			for record in Entrez.parse(mon_relevé):
				texte_à_écrire = '{:>15}    {:.57}'.format(" " + str(mon_Identifiant), str(record['Title']))
				print("%s" % texte_à_écrire)

		# Obtenons identifiant du gène souhaité / ou afficher plus de résultats / ou changer l'objet de recherche
		Identifiant = input("\nQuel ID à choisir? (tappez 'plus' pour plus de résultats, ou 'retour' pour \
		changer l'objet de recherche).. ")
		
		if Identifiant == "plus":
			numéro_RetStart = int(numéro_RetStart) + 10

		elif Identifiant == "retour":
			objet_récherche = input("Votre Recherche?.. ")

		else:
			nom_gène_intérêt = input("Nom du gène?.. ")
			break
	requête_Brut.close()


# Cherchons et enregistrer la séquence correspondant à notre Id
multifasta_var = []
print(">> Obtention des données...")
sequences_Brut = Entrez.efetch(db="Nucleotide", id=Identifiant, rettype="gb")
sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")

for chaqueSequence in sequences_Parsed:
	multifasta_var.append(chaqueSequence)
SeqIO.write(multifasta_var, str("/".join([directory,"gene_interet.fasta"])),"fasta")
print(">> Données enregistrés...")
sequences_Brut.close()


# Determinons le nom de l'espèce pour chaque « record » in fasta avec un expression régulier
for record in multifasta_var:
	espèce_initiale = record.description.replace("PREDICTED: ","")
	regex = r"^((?:\S+\s+){1}\S+).*"
	espèce_initiale = re.findall(regex, espèce_initiale)[0]
	espèce_gene_interet = espèce_initiale
	print("Espèce du gene initiale est : " + espèce_initiale)


# Envoyons le fasta de notre gène d'interet pour un BLAST chez NCBI
if skip_ncbi == "False":
	print(">> BLASTons notre gène d'interet au NCBI...")
	for record in SeqIO.parse(str("/".join([directory,"gene_interet.fasta"])),"fasta"):
		rétour_BLAST = NCBIWWW.qblast('blastn', 'nt', record.seq)
		fichier_BLAST = open(str("/".join([directory,"gene_interet_BLAST.xml"])), "w")
		fichier_BLAST.write(rétour_BLAST.read())
		rétour_BLAST.close()
		fichier_BLAST.close()


### Analyse des données
# Mise en place des dictionaires
Dictionnaire_ID_HSP = {}
Dictionnaire_Espèce_ID = {}
Dictionnaire_HSP_Espèce = {}
regexID = r"(?<=gi\|)(\d*)(?=\|)"
regexEspèce = r"^((?:\S+\s+){1}\S+).*"

for record in NCBIXML.parse(open(str("/".join([directory,"gene_interet_BLAST.xml"])))) :
    if record.alignments :
        for align in record.alignments :
            for hsp in align.hsps:
				# Determinons le nom de l'espèce
            	Nom_espèce = align.hit_def.replace("PREDICTED: ","")
            	Nom_espèce = re.findall(regexEspèce, Nom_espèce)[0]

				# Determinons l'identifiant
            	Identifiant = re.findall(regexID, align.hit_id)[0]
            	
            	# Ajoutons nom de l'espèce et identifiant à dictionnaire
            	Espèce_ID__entry = { Nom_espèce : Identifiant}
            	Dictionnaire_Espèce_ID.update(Espèce_ID__entry)

            	# Ajoutons bit-score et nom de l'espèce à dictionnaire
            	HSP_Espèce__entry = {hsp.bits : Nom_espèce}
            	Dictionnaire_HSP_Espèce.update(HSP_Espèce__entry)

            	# Ajoutons identifiant et bit-score à dictionnaire
            	ID_HSP__entry = {Identifiant : hsp.bits}
            	Dictionnaire_ID_HSP.update(ID_HSP__entry)






### Triage des données
HSP_liste = []
nom_espèces = []
HSP_filtrés = []

# Créer et trier une liste à partir des scores du dictionnaire
for value in Dictionnaire_ID_HSP.values():
	HSP_liste.append(value)

HSP_liste = sorted(HSP_liste, reverse=True)


# Creeons une liste identique sans duplicats ou résultats de l'espèce d'interet
for value in HSP_liste:
	if str(Dictionnaire_HSP_Espèce[value]) in espèce_initiale:
		pass
	else:
		if str(Dictionnaire_HSP_Espèce[value]) in nom_espèces:
			pass

		else:
			# Ajouter le nom de l'espèce à la liste nom_espèce
			nom_espèces.append(str(Dictionnaire_HSP_Espèce[value]))
			# Ajouter le score à la liste des scores HSP filtrés
			HSP_filtrés.append(value)


# Inversons le dictionnaire pour pouvoir utiliser le [key] avec les « values »
Dictionnaire_ID_HSP = {y:x for x,y in Dictionnaire_ID_HSP.items()}


# Ecrivons les titres dans le fichier texte : résultats.txt
fichier_résultats = open(str("/".join([directory,"résultats.txt"])), "w")
texte_à_écrire = '{:^9}     {:^25} {:^21}'.format("Bit-Score", "Espèce","Identifiant")
fichier_résultats.write("%s\n" % texte_à_écrire)
texte_à_écrire = '{:_<9}     {:_<25} {:_<21}'.format("", "","")
fichier_résultats.write("%s\n" % texte_à_écrire)
fichier_résultats.write("\n>> Résultats pour " + espèce_initiale + " :\n")


# Ecrivons ces premiers résultats dans le fichier résultats
Latex_table__initial = ""
for score in HSP_filtrés[:5]:
	texte_à_écrire = '{:>9}      {:25} {:21}'.format( str(float(score)), str(Dictionnaire_HSP_Espèce[score]), str(Dictionnaire_ID_HSP[score]))
	fichier_résultats.write("%s\n" % texte_à_écrire)
	Latex_table__initial = Latex_table__initial + ( str(float(score)) + " & " + Dictionnaire_HSP_Espèce[score] + " & " + Dictionnaire_ID_HSP[score] + "\\\\\n")
fichier_résultats.write("" + "\n")


# Mettons tous les identifiants ensemble pour que le NCBI les comprends 
for value in HSP_filtrés[:4]:
	Identifiant = ",".join([str(Dictionnaire_ID_HSP[value]),Identifiant])


# Exportons séquence du gène d'intêret sous format fasta
sequences_Brut = Entrez.efetch(db="Nucleotide", id=Identifiant, rettype="gb")
print(">> Données enregistrés...")
sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")

multifasta_var = []
for chaqueSequence in sequences_Parsed:
	multifasta_var.append(chaqueSequence)

# Sauvgardons le fichier multifasta et creeons le base de donnés qui servira ensuite pour le BLAST en local
SeqIO.write(multifasta_var, str("/".join([directory,"base_de_donnes.fasta"])),"fasta")
sequences_Brut.close()
cmd = "makeblastdb -in " + str("/".join([directory,"base_de_donnes.fasta"])) + " -parse_seqids -dbtype nucl"
try:
	subprocess.check_output(
		cmd,
		stderr=subprocess.STDOUT,
		shell=True)
except:
	os.system(cmd)


# Créons la liste « mes_identifiants » qui contiendra les séquences les plus proche à notre gène d'interet
mes_Identifiants = []
print("")
for value in HSP_filtrés[0:5]:
	mes_Identifiants.append(Dictionnaire_ID_HSP[value])

# On inverse de nouveau, c'est à dire qu'on revient aux clés/valeurs qu'on avait avant
Dictionnaire_ID_HSP = {y:x for x,y in Dictionnaire_ID_HSP.items()}


# Commencons le compteur pour pouvoir suivre le nombre de boucles qui ont lieu
compteur = 0
Latex_table = {}
ListeOrganismesInteret = {}

# Boucle pour blaster, et ecrire les résultats pour chacun de nos résultats du blast précedent
for mon_Identifiant in mes_Identifiants:
	# Initiation
	compteur +=1
	espèce_boucle = str(Dictionnaire_HSP_Espèce[Dictionnaire_ID_HSP[mon_Identifiant]])
	newdictValue = {compteur : espèce_boucle}
	ListeOrganismesInteret.update(newdictValue)
	message = "\n> Début du boucle nb." + str(compteur) + " sur l'espèce :"
	message = str("echo \"$(tput setaf 6)" + message + "$(tput sgr 0)\"")
	os.system(message)
	print(str(mon_Identifiant) + " : " + espèce_boucle)
	

	# Cherchons le séquence correspondant à notre Id
	sequences_Brut = Entrez.efetch(db="Nucleotide", id=mon_Identifiant, rettype="gb")
	print(">> Données enregistrés...")
	sequences_Parsed = SeqIO.parse(sequences_Brut,"gb")


	# Exportons séquence du gène d'intêret sous format fasta
	multifasta_var = []
	exportNameSeq = str(directory) + "/" + str(compteur) + "_" + espèce_boucle + ".fasta"

	for chaqueSequence in sequences_Parsed:
		multifasta_var.append(chaqueSequence)
	SeqIO.write(multifasta_var, exportNameSeq,"fasta")
	sequences_Brut.close()


	### Analyse des données
	print(">> Analyse BLAST des séquences...")
	
	# create unique export name
	exportNameBLAST = str(directory) + "/" + str(compteur) + "_" + espèce_boucle + "_BLAST.xml"

	# BLASTons en local
	rétour_BLAST = NcbiblastnCommandline(query=exportNameSeq,db=str("/".join([directory,"base_de_donnes.fasta"])), outfmt=5, out=exportNameBLAST)
	os.system(str(rétour_BLAST))


	# Mise en place des dictionaires
	boucle_Dictionnaire_ID_HSP = {}
	for record in NCBIXML.parse(open(exportNameBLAST)) :
		if record.alignments :
			for align in record.alignments :
				for hsp in align.hsps:
					# Determinons le nom de l'espèce
					Nom_espèce = align.hit_def.replace("PREDICTED: ","")
					Nom_espèce = re.findall(regexEspèce, Nom_espèce)[0]
					
					# Ajoutons bit-score et nom de l'espèce à dictionnaire
					HSP_Espèce__entry = {hsp.bits : Nom_espèce}
					Dictionnaire_HSP_Espèce.update(HSP_Espèce__entry)
					
					# Ajoutons identifiant et bit-score à dictionnaire
					ID_HSP__entry = {align.hit_id : hsp.bits}
					boucle_Dictionnaire_ID_HSP.update(ID_HSP__entry)	            		


	### Triage des données
	HSP_filtrés = []
	HSP_liste = []
	nom_espèces = []

	# Créer et trier une liste à partir des scores du dictionnaire
	for value in boucle_Dictionnaire_ID_HSP.values():
		HSP_liste.append(value)

	HSP_liste = sorted(HSP_liste, reverse=True)


	# Creeons une liste identique sans résultats de notre l'espèce d'interet de cette boucle
	for value in HSP_liste:
		if str(Dictionnaire_HSP_Espèce[value]) in espèce_boucle:
			pass
		else:
			HSP_filtrés.append(value)


	# Inversons le dictionnaire afin d'extraire l'identifiant à partir du HSP bit-score
	boucle_Dictionnaire_ID_HSP = {y:x for x,y in boucle_Dictionnaire_ID_HSP.items()}

	# Ecrivons les résultats pour ce boucle dans le fichier texte : résultats.txt
	HSPresults = []
	fichier_résultats = open(str("/".join([directory,"résultats.txt"])), "a")
	fichier_résultats.write(">> Résultats pour " + espèce_boucle + " :\n")

	# Obtenons l'identifiant à partir de HSP bit-score
	Dictionnaire_Accencion_ID = {}
	Latex_table__boucle = ""
	regexID = r"(?<=[a-z]{3}\|)(.*)(?=\|)"

	for score in HSP_filtrés:
		# Obtenons l'identifiant de l'espèce à partir du bit-score
		Identifiant_score = Dictionnaire_Espèce_ID[Dictionnaire_HSP_Espèce[score]]

		# Creeons une dictionnaire pour convertir numéro accencion en identifiant
		accencionID = re.findall(regexID, boucle_Dictionnaire_ID_HSP[score])[0]
		Accencion_ID__entry = {accencionID : Identifiant_score}
		Dictionnaire_Accencion_ID.update(Accencion_ID__entry)

		# Ecrivons résultats dans le fichier résultats.txt
		texte_à_écrire = '{:>9}      {:25} {:21}'.format(" " + str(float(score)), str(Dictionnaire_HSP_Espèce[score]),str(Identifiant_score))
		fichier_résultats.write("%s\n" % texte_à_écrire)

		# De même, ecrivons les résultats dans un dictionnaire qui mettra ces résultats dans le fichier LaTeX
		Latex_table__boucle = str(Latex_table__boucle) + (str(float(score)) + " & " + Dictionnaire_HSP_Espèce[score] + " & " + Identifiant_score + "\\\\\n")
	fichier_résultats.write("" + "\n")

	Latex_table.update({compteur : Latex_table__boucle}) 
	print("> fin du boucle N° " + str(compteur))

print("\n>> Resultats Enregistrés sous \"resultats.txt\"")





### Création de l'arbre phylogénétique
# Créeons le fichier newark pour l'arbre
multifasta = "base_de_donnes.fasta"
guidetree_accession = str("/".join([directory,"guidetree_accession.txt"]))
clustalomega_cline = ClustalOmegaCommandline(infile=multifasta, outfile="base_de_donnes_alignes.fasta", guidetree_out="guidetree_accession.txt",verbose=True, auto=True, force=True)
cmd = "cd runFiles\n" + str(clustalomega_cline)
cmd = str(cmd)

# Lancer le shell avec subprocess pour cacher le stdout pour une interface plus propre, si possible
print(">> Lancons le ClustalOmega sur nos séquences...")
try:
	subprocess.check_output(
		cmd,
		stderr=subprocess.STDOUT,
		shell=True)
except:
	os.system(cmd)


# Remplacons le numero accencion avec le nom de l'espèce dans le guidetree_binomiale
guidetree_binomiale = str("/".join([directory,"guidetree_binomiale.txt"]))
with open(guidetree_accession, 'r') as file :
  filedata = file.read()

Dictionnaire_Espèce_ID = {y:x for x,y in Dictionnaire_Espèce_ID.items()}
for key in Dictionnaire_Accencion_ID.keys():
	especeName = Dictionnaire_Espèce_ID[Dictionnaire_Accencion_ID[key]]
	especeName = especeName.replace(" ", "_")
	filedata = filedata.replace(key, especeName)

with open(guidetree_binomiale, 'w+') as file :
  file.write(filedata)


# Créons l'arbre à partir du guidetree_binomiale
arbre_in = Phylo.read(guidetree_binomiale, "newick")
arbre_type = ""

# Essayer de produire le beau arbre en png
print(">> Générons notre arbre phylogénétique")
try:
	Phylo.draw(arbre_in,do_show=None)
	figure = plt.gcf()
	figure.set_size_inches(11, 8.5)
	plt.savefig(LaTexdir + "/monarbre.png", dpi = 100)
	arbre_type = "image"

# Si pas possible, produit arbre ASCII
except:
	arbre_out = io.StringIO()
	Phylo.draw_ascii(arbre_in,file=arbre_out,column_width=90)
	arbre_ascii = str(arbre_out.getvalue())
	arbre_out.close()
	arbre_type = "ASCII"





### Creons le fichier LaTex pour afficher les résultats
with open("Resultat_Template.tex", 'r') as file :
  resultat_latex = file.read()

# Creer le titre
resultat_latex = resultat_latex.replace("[NomUtilisateur*]", nom_utilisateur)
resultat_latex = resultat_latex.replace("[emailUtilisateur*]", Entrez.email)
resultat_latex = resultat_latex.replace("[geneInteret*]", nom_gène_intérêt)

# Creer les Tableaux dans le fichier LaTex grace à notre dictionnaire Latex_table
resultat_latex = resultat_latex.replace("[TableaugeneInteret*]", str(Latex_table__initial))
resultat_latex = resultat_latex.replace("[especegeneInteret*]", espèce_gene_interet)

resultat_latex = resultat_latex.replace("[Tableaugene1*]",str(Latex_table[1]))
resultat_latex = resultat_latex.replace("[especegene1*]",ListeOrganismesInteret[1])

resultat_latex = resultat_latex.replace("[Tableaugene2*]",str(Latex_table[2]))
resultat_latex = resultat_latex.replace("[especegene2*]",ListeOrganismesInteret[2])

resultat_latex = resultat_latex.replace("[Tableaugene3*]",str(Latex_table[3]))
resultat_latex = resultat_latex.replace("[especegene3*]",ListeOrganismesInteret[3])

resultat_latex = resultat_latex.replace("[Tableaugene4*]",str(Latex_table[4]))
resultat_latex = resultat_latex.replace("[especegene4*]",ListeOrganismesInteret[4])

if arbre_type == "image":
	resultat_latex = resultat_latex.replace("[monarbrephylo*]", "\\includegraphics[width=\columnwidth,keepaspectratio]{monarbre.png}")
else:
	resultat_latex = resultat_latex.replace("[monarbrephylo*]", "\\begin{verbatim}" + arbre_ascii + "\end{verbatim}")

# Écrivons ces modification à notre Fichier Resultats
with open("LaTexResults/Resultats.tex", 'w+') as file :
  file.write(resultat_latex)

# Essayons de compiler le .tex
try:
	print(">> Compilons le fichier resultats en LaTeX..")
	try:
		cmd = "cd LaTexResults\npdflatex --shell-escape --file-line-error --synctex=1 Resultats.tex"
		subprocess.check_output(
			cmd,
			#stderr=subprocess.STDOUT,
			shell=True)

	except:
		os.system(cmd)
	print(">> Resultats Enregistrés sous \"LaTexResults/Resultats.tex\"")
except:
	print("Erreur en compilant le fichier 'LaTexResults/Resultats.tex' ; veuillez-voir 'runFiles/resultats.txt'")
