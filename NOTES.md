# Notes


## Propositions

- There're certain genetic characteristics that control the antropophily of a mosquito specie.
- Species with large genetic distance and similar anthropophily probably have high identity for the genes that affect antropophily.
- Species with small genetic distance and differing anthropophily probably have low identity for the genes that affect antropophily.



## Approach ideas scratchpad
- Get species in taxonomic subtree of a taxon
- Download protein (?) sequences (result: set of fasta files)
- Write a script that produces a gene correspondence table: (specie1, specie2, s1_gene_id, s2_gene_id, identity)
- Calculate "relevance score" for each gene, something like identity/<average identity for other genes of the same species>. The score should find genes that are different in species that are similar.
- Get scores for antropophily for species.
- Find correlations between most relevant genes and antropophily scores.
- Can we use rank lowering technique such as Singular Value Decomposition (similar to Latent Semantic Analysis) to find the most important genes.
- We have too many combinations of species (51 species with available proteome => 1200+ possible pairs). We have to compare only a subset comprising of the pairs that are most similar.
	- Find what is the minimal random subset of the proteome that would indicate reliably the overall identity.
	- Ignore species without good data for antropophily.
	- Compare all pairs using this minimal indicative subset.
	- Get the N best ones and analyse them using their whole proteomes.
- From the correspondence table make a relevance table. The relevance score can be something like (id(Prot_A,i, Prot_B,j)/id(Genome_A, Genome_B))^(anthrop(A)*anthrop(b)), where anthrop(x) gives the anthropophily of the species as 1 for high and -1 for low.