Recreating the Clustering algorithm experiments from the article: "The what, The from, and The to: The Migration Games in Deduplicated Systems"
Roei Kisous, Ariel Kolikant, Abhinav Duggal, Sarai Sheinvald, Gala Yadgar

The purpose of this is to explain how to recreate the experiments done with the HC algorithm.
The algorithm creates a bit appearance matrix based on the system snapshots provided, where each file represents a row and each column represents a block, 1 is set if the file contains the block, 0 otherwise.
With this matrix we generate our dissimilarity matrix, where each row and column correspond to a file, and the values represent the Jaccard distance between them. As the matrix is symmetric, only the lower triangular is used for the algorithm; the upper triangular is used to restore the original values if multiple runs are performed.
Using the parameters set, these value are changed in accordance with the explanation in the article. The algorithm performs hierarchical clustering and outputs the results into results files. 

For more information, please read 'HCexperiments.pdf'.