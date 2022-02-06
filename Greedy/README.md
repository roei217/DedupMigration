# Greedy

General Information
-------------------
Greedy is a migration plan creator for deduplicated storage systems.

Running
-------

1. Install and configure libboost

1. Complie the code
   > make

2. Run Greedy
   > ./GreedyLoadBalancerUnited {volumelist} {output} {summaryFile} {timelimit} {Traffic} {margin}
		
    1)	Volumelist: volume List section
    2)	Output: the path to the where the csv migration plan is to be written.
    3)	summaryFile: path to the summary file.
    4)	TimeLimit: after how many seconds should the algorithm stop
    5)	Traffic: what is the maximum amount of traffic allowed for the algorithm (0-100, representing percentage)
    6)	Margin: the error margin allowed from the desired balance. If this margin is 1 or more, the code would run without load balancing. (0-1, representing fraction)

Examples
------------------
Find a migration plan for examples/basicLoadBalance_k1_3Vols, saving the summary in summary.csv, the migration plan cannot use more than 20% of the original size as traffic and must be within margin of 0.3 from the desired load balance and cannot use more than 6 houres (21600 seconds). save the plan as outputFile.csv.

  > ./GreedyLoadBalancerUnited basicLoadBalance_k1_3Vols outputFile.csv summaryFile.csv 21600 20 0.3
  
Paper
------
To read more about Greedy please check our [paper](https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous, Ariel Kolikant, Gala Yadgar (Technion - Israel Institute of Technology);
           Abhinav Duggal (Dell Technologies);
           Sarai Sheinvald (ORT Braude College of Engineering).

Email : sarielko@campus.technion.ac.il

For more information, please read 'Greedy.pdf'.
