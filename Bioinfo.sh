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


# Obtenir l'objet de recherche
	echo "$(tput setaf 5)####$(tput sgr 0) Recherche du projet ou personalisé ? (projet/custom)"
	# read default_search
		while true; do
		read default_search
		default_search=$(echo "$default_search" | tr '[:upper:]' '[:lower:]')
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



# Obtenir addresse mail de l'utilisateur
	echo "$(tput setaf 5)####$(tput sgr 0) Quelle est votre addresse email ?"
	read email


echo " "
echo " "


echo "Lancons le .py..."

python3 BioInfo.py "$default_search" "$email" "$skip_predictions" "$blast_type"