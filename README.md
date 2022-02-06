Paper
------
To read more about this please check our *paper*(https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous, Ariel Kolikant, Gala Yadgar (Technion - Israel Institute of Technology);
           Abhinav Duggal (Dell Technologies);
           Sarai Sheinvald (ORT Braude College of Engineering).
		   
Clustering
------
Clustering is a well-known technique for grouping objects based on their similarity. 
It is fast and effective, and is directly applicable to our domain: files are similar if they share a large portion of their blocks. 
Our approach is thus to create clusters of similar files and to assign each cluster to a volume, remapping those files that were assigned to a volume different from their original location.

For more information, please read 'Clustering/HCexperiments.pdf'.

ILP
------
The purpose of this is to explain how to recreate the experiments done with the ILP algorithm.
The algorithm creates an ILP model representing legal migration plans and utilizes the Gurobi ILP solver to find the model's optimal solution. The solution is then translated to a migration plan.The constraints of the model include contraints to maintain legal migrations and constraints to limit traffic usage and to achieve load balancing. In a case where no legal solution exists, none is returned. In a case where the optimal solution has not been reached within the given time, the best solution so far is given.

For more information, please read 'ILP/ILPExperiments.pdf'.

Greedy
------
The purpose of this is to explain how to recreate the experiments done with the Greedy algorithm.
The algorithm works in iterations, choosing the best file to migrate for every iteration. The defenition of best changes in accordance to the two types of phases; 1) the balancing phase attempts to achieve load balancing and prefers moving files from the biggest volume to the smallest volume 2) the optimization phase prefers to migrate files resulting in the most data removal. Greedy performs 5 iteration of each phase.

For more information, please read 'Greedy/GreedyExperiments.pdf'.

CostCalculator
------
For evaluating the nature of a specific migration plan, we implemented a calculator that outputs the costs and requirements.
Those costs dictate the nature of a migration plan.

For more information, please read 'CostCalculator/CostCalc.pdf'.


