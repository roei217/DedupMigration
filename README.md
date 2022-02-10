Paper
------
To read more about this please check our *paper*(https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous and Ariel Kolikant, Technion - Israel Institute of Technology; 
		 Abhinav Duggal, DELL EMC; 
		 Sarai Sheinvald, ORT Braude College of Engineering; 
		 Gala Yadgar, Technion - Israel Institute of Technology
		   
Clustering
------
Clustering is a well-known technique for grouping objects based on their similarity. 
It is fast and effective, and is directly applicable to our domain: files are similar if they share a large portion of their blocks. 
Our approach is thus to create clusters of similar files and to assign each cluster to a volume, remapping those files that were assigned to a volume different from their original location.

For more information, please read 'Clustering/Clustering.pdf'.

ILP
------
The ILP algorithm constructs an ILP model for migration plans from the input deduplicated system, solves the ILP model using the Gurobi ILP solver, then translates the solution into a migration plan. The ILP solver finds the optimal solution for the model unless it runs out of time. If a timeout occurs, the best solution found before the timeout is returned. Although the optimal solution of the model is not the optimal solution of the migration plan, it is shown to be close. 

For more information, please read 'ILP/ILP.pdf'.

Greedy
------
The greedy algorithm works in iterations, in which the "best" file migration is selected for each iteration. Iteration type influences the method for calculating the best migration. During balancing iterations, the best migration will be the one that brings the system closer to balance; during optimization phases, the best migration will be the one that reduces the total system size without breaking the load balance constraints.

For more information, please read 'Greedy/Greedy.pdf'.

CostCalculator
------
For evaluating the nature of a specific migration plan, we implemented a calculator that outputs the costs and requirements.
Those costs dictate the nature of a migration plan.

For more information, please read 'CostCalculator/CostCalc.pdf'.


