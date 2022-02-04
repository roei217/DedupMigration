Recreating the algorithms' experiments from the article: "The what, The from, and The to: The Migration Games in Deduplicated Systems"
Roei Kisous, Ariel Kolikant, Abhinav Duggal, Sarai Sheinvald, Gala Yadgar

There are 4 folders in this project:
1. Clustering
The purpose of this is to explain how to recreate the experiments done with the HC algorithm.
The algorithm creates a bit appearance matrix based on the system snapshots provided, where each file represents a row and each column represents a block, 1 is set if the file contains the block, 0 otherwise.
With this matrix we generate our dissimilarity matrix, where each row and column correspond to a file, and the values represent the Jaccard distance between them. As the matrix is symmetric, only the lower triangular is used for the algorithm; the upper triangular is used to restore the original values if multiple runs are performed.
Using the parameters set, these value are changed in accordance with the explanation in the article. The algorithm performs hierarchical clustering and outputs the results into results files. 

For more information, please read 'Clustering/HCexperiments.pdf'.

2. ILP
The purpose of this is to explain how to recreate the experiments done with the ILP algorithm.
The algorithm creates an ILP model representing legal migration plans and utilizes the Gurobi ILP solver to find the model's optimal solution. The solution is then translated to a migration plan.The constraints of the model include contraints to maintain legal migrations and constraints to limit traffic usage and to achieve load balancing. In a case where no legal solution exists, none is returned. In a case where the optimal solution has not been reached within the given time, the best solution so far is given.

For more information, please read 'ILP/ILPExperiments.pdf'.

3. Greedy
The purpose of this is to explain how to recreate the experiments done with the Greedy algorithm.
The algorithm works in iterations, choosing the best file to migrate for every iteration. The defenition of best changes in accordance to the two types of phases; 1) the balancing phase attempts to achieve load balancing and prefers moving files from the biggest volume to the smallest volume 2) the optimization phase prefers to migrate files resulting in the most data removal. Greedy performs 5 iteration of each phase.

For more information, please read 'Greedy/GreedyExperiments.pdf'.

4. CostCalculator
For evaluating the nature of a specific migration plan, we implemented a calculator that outputs the costs and requirements.
From this plan, the calculator provides the traffic consumption, score achieved, and deletion obtained while also providing the initial and final size of each volume with specific traffic and deletion. Each volume is iterated over and its incoming traffic and resulting deletion is calculated. The final traffic is equal to the sum of all traffic calculated, and the final deletion is equal to the sum of all deletions minus the sum of all traffic.
We added a results cache to save time, so that in case a migration plan needs to be calculated again (which we might not know about), it is not calculated again and instead the results are retrieved from the cache. The migration plan's hash is inserted as a key into the cache database with the corresponding results. The hash is recalculated and searched in the database when trying to find this plan again. It is returned if it is found, otherwise it is recalculated. The results will be fetched if the results are in the db specified in the -cache_path tag and the -no_cache tag was not used.
Our file-based lock allows the calculator to run on multiple threads concurrently without user interaction. There is no restriction on the number of threads the user may run as long as he does not delete the txt file created under the folder in which the database is stored, this is the file lock.
Furthermore, to test our plan, we use the calculator both on the sampled system (for real-time evaluation) and on the original system (for static evaluation), as shown below.

For more information, please read 'CostCalculator/CostCalc.pdf'.
