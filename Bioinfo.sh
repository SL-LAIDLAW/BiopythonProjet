#!/bin/bash

# Header
 	echo "******************************************************************************"
	echo " "
	echo "$(tput setaf 0)$(tput setab 0)***$(tput sgr 0)"
	echo "$(tput setaf 0)$(tput setab 0)***$(tput sgr 0)   $(tput setaf 0)$(tput bold)HLIN404 Projet 2017                             $(tput setaf 0)$(tput bold)Introduction généraliste$(tput sgr 0)"
	echo "$(tput setaf 0)$(tput setab 0)   $(tput sgr 0)   $(tput setaf 0)$(tput bold)                                                        à l'informatique$(tput sgr 0)"
	echo "$(tput setaf 0)$(tput setab 0)   $(tput sgr 0)   $(tput setaf 5)$(tput bold)Sean L. LAIDLAW                               Application à la génomique$(tput sgr 0)"
	echo "$(tput setaf 0)$(tput setab 0)***$(tput sgr 0)"
	echo " "


# Obtenir terme de recherche
	echo "$(tput setaf 5)####$(tput sgr 0) Votre recherche? (tappez 'd' pour default, c'est à dire, ma requete du projet)"
	# read search_term
		while true; do
		read search_term
		if [ "$search_term" = "d" ]
			then
				search_term="dec2 mRNA for bhLH protein complete cds "
		elif [ "$search_term" = "" ]
			then
				echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, veuillez entrer des termes de recherche$(tput sgr 0)"
				continue
		fi
		break
	done
echo " "
echo " "


# Obtenir nature de séquence
	echo "$(tput setaf 5)####$(tput sgr 0)Quel est la nature de votre séquence, nucleotidique ou proteique (n/p)"
		while true; do
			read format
			# transforme la case de notre variable 'format'
			format=$(echo "$format" | tr '[:upper:]' '[:lower:]')
			if [ "$format" = "n" ]
				then
				format="nucleotide"
			elif [ "$format" = "p" ]
				then
				format="Protein"
			elif [ "$format" != "n" ] || [ "$format" != "p" ]
			then
				echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, repondez par \"n\" ou \"p\"$(tput sgr 0)"
				continue
		fi
		break
	done


echo " "
echo " "



# Obtenir addresse mail de l'utilisateur
	echo "$(tput setaf 5)####$(tput sgr 0) Êtes-vous Sean Laidlaw? (y/n)"
		while true; do
			read identify
			if [ "$identify" = "n" ] || [ "$identify" = "N" ]
				then
					echo "$(tput setaf 6)######$(tput sgr 0) Quelle est votre addresse email ?"
					read email

			elif [ "$identify" = "y" ] || [ "$identify" = "Y" ]
				then email="sean.laidlaw@etu.umonptellier.fr"

			elif [ "$identify" != "y" ] || [ "$identify" != "n" ]
				then
					echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, repondez par \"y\" ou \"n\"$(tput sgr 0)"
					continue
		fi
		break
	done


echo "Lancons le .py : cette opération est suseptible de prendre quelques minutes..."
python Bioinfo.py "$search_term" "$format" "$email"