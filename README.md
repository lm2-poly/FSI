# MEC6514 - Interactions fluide-structure

## Professeur : Frédérick Gosselin

Vous trouverez ci-joint du contenu complémentaire au cours donné en classe, avec des rappels sur la matière couverte, des détails supplémentaires et des codes vous permettant d'approfondir et d'obtenir les résultats vus en classe. Une partie de la documentation est écrite en anglais, et le reste en français.

Plusieurs portions ont été couvertes, et pour chacune d'elle, vous trouverez un Jupyter Notebook, qu'il est possible :
- De consulter en ligne, en cliquant sur le lien associé ci-dessous (quelques bogues d'affichage des équations sont présents, pour une meilleure lisibilité, se réferer aux options suivantes) ;
- D'exécuter localement au moyen d'une installation Anaconda et avec les librairies précisées dans le fichier ``` requirements_python3-10-4.txt ``` (cf. la procédure ci-dessous) ;
- D'exécuter en ligne, sans installation nécessaire, en cliquant sur [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lm2-poly/FSI/HEAD)

### Analyse dimensionnelle 

[Masse ajoutée - Cylindre](https://github.com/lm2-poly/FSI/blob/main/Chapitre-2_Analyse-dimensionnelle/1_Masse_ajoutee/Cylindre/Masse_ajoutee.ipynb)

[Masse ajoutée - Disque et verre](https://github.com/lm2-poly/FSI/blob/main/Chapitre-2_Analyse-dimensionnelle/1_Masse_ajoutee/Disque_Verre/added_mass.ipynb)

[Analyse de stabilité - Système masse-ressort](https://github.com/lm2-poly/FSI/blob/main/Chapitre-2_Analyse-dimensionnelle/2_Analyse_stabilite/Analyse_stabilite.ipynb)

### Vibrations induites par vortex

[VIV - Cylindre](https://github.com/lm2-poly/FSI/blob/main/Chapitre-3_Vibrations-induites-par-vortex-(VIV)/VIV.ipynb)

### Ecoulements transverses sur structures élancées

[Flottement d'une aile](https://github.com/lm2-poly/FSI/blob/main/Chapitre-4_Ecoulements-transverses-sur-structures-elancees/Flottement/Flottement.ipynb)

### Ecoulement axial

[Modes d'un tuyau](https://github.com/lm2-poly/FSI/blob/main/Chapitre-5_Ecoulement-axial/Tuyau/Modes/Modes.ipynb)

[Vitesse critique d'un tuyau](https://github.com/lm2-poly/FSI/blob/main/Chapitre-5_Ecoulement-axial/Tuyau/Vitesse-critique_flottement/Vitesse_critique.ipynb)

### Exécution locale des Jupyter Notebooks
- Installer Anaconda au moyen du lien suivant : https://www.anaconda.com/products/distribution ;
- Ouvrir Anaconda Prompt à partir du menu Démarrer ;
- Se rendre dans le dossier de votre choix avec le chemin vers celui-ci : ``` cd Utilisateurs/Utilisateur/Documents/FSI ``` ;
- Y déposer le dossier téléchargé depuis Github ;
- Créer un environnement avec les librairies nécessaires : ``` conda create --name name --file requirements.txt ``` ;
- Activer l'environnement : ``` conda activate name ``` (et pour le désactiver ``` conda deactivate ```) ;
- Exécuter Jupyter Notebook : ``` jupyter notebook ``` ;
- Depuis l'onglet ouvert dans le navigateur internet, choisissez le fichier à ouvrir.
