Recreating the ILP algorithm experiments from the article: "The what, the from, and the to: Migration Games in Deduplicated Systems"
Roei Kisous, Ariel Kolikant, Abhinav Duggal, Sarai Sheinvald, Gala Yadgar

The purpose of this is to explain how to recreate the experiments done with the ILP algorithm.
The algorithm creates an ILP model representing legal migration plans and utilizes the Gurobi ILP solver to find the model's optimal solution. The solution is then translated to a migration plan.The constraints of the model include contraints to maintain legal migrations and constraints to limit traffic usage and to achieve load balancing. In a case where no legal solution exists, none is returned. In a case where the optimal solution has not been reached within the given time, the best solution so far is given.

For more information, please read 'ILPExperiments.pdf'.
