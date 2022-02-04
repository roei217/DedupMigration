Recreating the Greedy algorithm experiments from the article: "The what, the from, and the to: Migration Games in Deduplicated Systems"
Roei Kisous, Ariel Kolikant, Abhinav Duggal, Sarai Sheinvald, Gala Yadgar

The purpose of this is to explain how to recreate the experiments done with the Greedy algorithm.
The algorithm works in iterations, choosing the best file to migrate for every iteration. The defenition of best changes in accordance to the two types of phases; 1) the balancing phase attempts to achieve load balancing and prefers moving files from the biggest volume to the smallest volume 2) the optimization phase prefers to migrate files resulting in the most data removal. Greedy performs 5 iteration of each phase.

For more information, please read 'GreedyExperiments.pdf'.
