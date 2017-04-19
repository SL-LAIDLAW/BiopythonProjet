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
	echo "$(tput setaf 5)####$(tput sgr 0) Recherche du projet ou personalisé? (projet/custom)"
	# read default_search
		while true; do
		read default_search
		if [ "$default_search" = "projet" ]
			then
				default_search="True"
		elif [ "$default_search" = "custom" ]
			then
				default_search="False"
		else
			echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, veuillez entrer \"projet\" ou \"custom\"$(tput sgr 0)"
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


echo " "
echo " "


# Garde resultats des Humains?
	echo "$(tput setaf 5)####$(tput sgr 0) Est-ce qu'on garde les résultats des humains (y/n)"
		while true; do
			read skip_humans
			if [ "$skip_humans" = "n" ] || [ "$skip_humans" = "N" ]
				then skip_humans="True"
			elif [ "$skip_humans" = "y" ] || [ "$skip_humans" = "Y" ]
				then skip_humans="False"
			elif [ "$skip_humans" != "y" ] || [ "$skip_humans" != "n" ]
				then
					echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, repondez par \"y\" ou \"n\"$(tput sgr 0)"
					continue
		fi
		break
	done


echo " "
echo " "


# Garde resultats des Humains?
	echo "$(tput setaf 5)####$(tput sgr 0) Est-ce qu'on garde les prédictions? (y/n)"
		while true; do
			read skip_predictions
			if [ "$skip_predictions" = "n" ] || [ "$skip_predictions" = "N" ]
				then skip_predictions="True"
			elif [ "$skip_predictions" = "y" ] || [ "$skip_predictions" = "Y" ]
				then skip_predictions="False"
			elif [ "$skip_predictions" != "y" ] || [ "$skip_predictions" != "n" ]
				then
					echo "$(tput setaf 1)Erreur : la réponse n'est pas valide, repondez par \"y\" ou \"n\"$(tput sgr 0)"
					continue
		fi
		break
	done


echo "Lancons le .py..."

python Bioinfo.py "$default_search" "$format" "$email" "$skip_humans" "$skip_predictions"