# MLA_Projet
Ce repository  sert de plateforme commune pour les projets de Modèles de Localisation et Applications (MLA).
<u>Attention : </u> Il faut me contacter afin d'être en collaborateur sur le Git (Val Bc sur messenger)
## Environnement
julia version 1.8.5
Bibliothèques :
- JuMP
- Gurobi

## Utilisation
- Installation des bibliothèques 
```bash
julia ./requirement.jl
```

- Run
```bash
julia ./src/main.jl
```
## Structure du projet
- Le dossier `data/`contient toutes les instances.
- `requirement.jl` est un fichier d'installation et de mise à jour des différentes bibliothèques nécessaires.
- Le dossier `src/`contient l'ensemble du code du projet.
- Le fichier `src/main.jl` est le fichier principal du projet qui contient les différentes utilisations des fonctions.
- Le fichier `src/readData.jl` gère le parsing des données.
- Les fichiers `src/groupeX.jl` contiennent le code fourni par le groupe X.

## Consignes
- Chaque groupe n'utilise et ne modifie que son fichier Julia. Si vous avez besoin de créer d'autres fichiers, soyez attentifs et assurez-vous d'une notation claire.
- En cas de modification du fichier d'un autre groupe, assurez-vous de les en informer au préalable.
- Le fichier `main.jl` doit rester propre et léger.
- Ne pas push d'autres branches et restez sur la branche `dev`
- Maintenez à jour le `requirement.jl` si vous utilisez des bibliothèques spéciales.
- Assurez-vous de mettre des commits précis décrivant les changements apportés.