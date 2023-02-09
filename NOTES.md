# Notes


## Propositions

- There're certain genetic characteristics that control the anthropophily of a mosquito species.
- Species with large genetic distance and similar anthropophily probably have higher than average identity for the genes that affect anthropophily.
- Species with small genetic distance and differing anthropophily probably have lower than average identity for the genes that affect anthropophily.



## Approach ideas scratchpad
- Get species in taxonomic subtree of a taxon
- Download protein (?) sequences (result: set of fasta files)
- Write a script that produces a gene correspondence table: (specie1, specie2, s1_gene_id, s2_gene_id, identity)
- Calculate "relevance score" for each gene, something like identity/<average identity for other genes of the same species>. The score should find genes that are different in species that are similar.
- Get scores for anthropophily for species.
- Find correlations between most relevant genes and anthropophily scores.
- Can we use rank lowering technique such as Singular Value Decomposition (similar to Latent Semantic Analysis) to find the most important genes.
- We have too many combinations of species (51 species with available proteome => 1200+ possible pairs). We have to compare only a subset comprising of the pairs that are most similar.
	- Find what is the minimal random subset of the proteome that would indicate reliably the overall identity.
	- Ignore species without good data for antropophily.
	- Compare all pairs using this minimal indicative subset.
	- Get the N best ones and analyse them using their whole proteomes.
- From the correspondence table make a relevance table. The relevance score can be something like (id(Prot_A,i, Prot_B,j)/id(Genome_A, Genome_B))^(anthrop(A)*anthrop(b)), where anthrop(x) gives the anthropophily of the species as 1 for high and -1 for low.


## TODO:

- Check if community detection is reasonable.
	- e.g. get some known families and check if they are correctly identified by the community detection algorithm.
    - if we have annotations for the genes, can we use them to check if annotations within a community are homogeneus.

- Enrichment analysis.
	- Check if mapping genes to e.g. Drosophila and doing the enrichment analysis there is a good idea. Maybe because it has more annotations.

- Get top candidate proteins and then do further analysis for them on genetic level.
	- Also, we can include more species outside the set of species with available proteome.

- Use the same technique on a well-studied phenotype in another family of species and see if it is able to detect the relevant genes.

- Take each protein, find the set of its orthologs in all the species, calculate a relevancy score for each set, this score is assigned to the protein from which the set originated. Rank the proteins for each species and check if orthologs of the same protein are near the top in most of them. (this is a way to remove the necessity of community detection, which may be introducing errors)

- Separate species into groups depending on anthropophily (negative, ambivalent and positive) and construct consensus sequences from the high-ranking proteins. You can predict the anthropophily of new species by checking to which consensus sequences they align better. This can be used to validate the conjecture without modifying mosquiotoes genetically. 