# mosquitoes
An experiment in trying to find the genetic factors that contribute to the anthropophily (preference for feeding on humans) of mosquito species.


In short, we try to test the following propositions:

- There're certain genetic characteristics that control the anthropophily of a mosquito species.
- Species with large genetic distance and similar anthropophily probably have higher than average identity for the genes that affect anthropophily.
- Species with small genetic distance and differing anthropophily probably have lower than average identity for the genes that affect anthropophily.

The project is for now based on proteome analysis.

Essentially, the approach is:
1. Get the proteome of each mosquito species for which we can.
2. For each pair of species, get each protein in their proteomes and find the closest homolog from the other species along with the similarity score. This gives us a correspondence table with columns: $(species_1, species_2, protein_1, protein_2, score)$
3. For each row in the corrospondence table (pair of homologous proteins) we calculate a relevance score: _relevance = (genes_similarity/species_similarity)_<sup>anthropophily_similarity</sup>, where _anthropophily_similarity_ is in the range $[-1,1]$ (-1 for high dissimilarity and 1 for high similarity).
4. Each protein is assigned the average relevance score from all homologous pairs it is part of.
5. Get the top N proteins by relevance.

## Prerequisites:
You have to have python3 installed, with version >= 3.8.

You have to have the **diamond** executable in the your PATH.

You can download it from: https://github.com/bbuchfink/diamond/releases and decompress it.
If you're using bash, you can make it accessible by: 
```bash
$ mkdir ~/bin
$ mv path/to/diamond ~/bin
$ chmod +x ~/bin/diamond

# then add this to your ~/.bashrc (or .zshrc if you're using zsh)
export PATH=$PATH:~/bin

# Make sure that the changes are sourced
$ source ~/.bashrc
# or source ~/.zshrc (respectively)
```

You also need the **clustal omega** executable in the your PATH. On Ubuntu:
```bash
sudo apt-get install clustalo
```

## Setup

```bash
$ python -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```