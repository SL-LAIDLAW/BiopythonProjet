#!/bin/bash


# Header
	echo " "
	echo " "
	echo "$(tput setaf 0)$(tput setab 0)***************************************************************************$(tput sgr 0)"
	# echo " - - - - - - - - - - - - - - HLIN404 Projet 2017 - - - - - - - - - - - - - "
	echo "$(tput setaf 0)$(tput setab 0) - - - - - - - - - - - - - - $(tput setaf 7)HLIN404 Projet 2017$(tput sgr 0)$(tput setaf 0)$(tput setab 0) - - - - - - - - - - - - - $(tput sgr 0)"
	# echo " - - - - - - - - - - - - - - - - S.L. LAIDLAW  - - - - - - - - - - - - - - "
	echo "$(tput setaf 0)$(tput setab 0) - - - - - - - - - - - - - - - - $(tput sgr 0)$(tput setaf 5)$(tput setab 0)$(tput bold)S.L. LAIDLAW$(tput sgr 0)$(tput setaf 0)$(tput setab 0)  - - - - - - - - - - - - - - $(tput sgr 0)"
	echo "$(tput setaf 0)$(tput setab 0)***************************************************************************$(tput sgr 0)"
	echo " "
	echo " "
	echo " "

# Obtenir terme de recherche
	echo "#### Terme de recherche? (tappez 'd' pour default, c'est à dire, ma requete du projet)"
	# read search_term
		while true; do
		read search_term
		if [ "$search_term" = "d" ]
			then
				search_term="dec2 mRNA for bhLH protein complete cds "
		elif [ "$search_term" = "" ]
			then
				echo "$(tput setaf 1)response n'est pas valide, veuillez entrer des termes de recherche svp$(tput sgr 0)"
				continue
		fi
		break
	done
echo " "
echo " "


# Obtenir nature de séquence
	echo "#### Quel est la nature de votre séquence, nucleotidique ou proteique (n/p)"
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
				echo "$(tput setaf 1)Erreur : response n'est pas valide, repondez par \"n\" ou \"p\" svp$(tput sgr 0)"
				continue
		fi
		break
	done


echo " "
echo " "



# Obtenir addresse mail de l'utilisateur
	echo "#### Etes-vous Sean Laidlaw? (y/n)"
		while true; do
			read identify
			if [ "$identify" = "n" ] || [ "$identify" = "N" ]
				then
					echo "###### Quel est votre addresse mail ?"
					read email

			elif [ "$identify" = "y" ] || [ "$identify" = "Y" ]
				then email="sean.laidlaw@etu.umonptellier.fr"

			elif [ "$identify" != "y" ] || [ "$identify" != "n" ]
				then
					echo "$(tput setaf 1)Erreur : response n'est pas valide, repondez par \"y\" ou \"n\" svp$(tput sgr 0)"
					continue
		fi
		break
	done


python Bioinfo.py "$search_term" "$format" "$email"