# ILP

General Information
-------------------
ILP is a migration plan creator for deduplicated storage systems.

Running
-------
1. install and configure gurobi and libboost

2. Complie the code
   > ./compile_Thanos-united.sh

3. Run ILP
   > ./Thanos-United {volumeList} {summaryFile} { Traffic } {margin} {innerFilter} {outDir} {timeLimit} {seed} {num of threads}
		
    1)	Volumelist: volume List section
    2)	summaryFile: path to the summary file.
    3)	Traffic: what is the maximum amount of traffic allowed for the algorithm (0-100, representing percentage)
    4)	Margin: the error margin allowed from the desired balance. If this margin is 1 or more, the code would run without load balancing. (0-1, representing fraction)
    5)	innerFilter: a huge relaxation of the load balancing constraint. It is recommended to be left at 0 unless significant time improvements are needed (with accuracy reduction being accepted)
    6)	outDir: a directory where the output files would be stored
    7)	timeLimit: how many seconds is the gurobi solver allowed to run
    8)	seed: the random seed given to the gurobi runner, different models can find better/faster solutions under different seeds
    9)	num of threads: how many threads could gurobi use. Recommended to use just a bit bellow the server’s maximum amount of CPU’s

Examples
------------------
Find a migration plan for examples/basicLoadBalance_k1_3Vols, saving the summary in summary.csv, the migration plan cannot use more than 20% of the original size as traffic and must be within margin of 0.25 from the desired load balance and cannot run mmore than 6 hours (21600 seconds). use all blocks without inner filter, random seed of 10 and 10 threads, saving the output into the outDir folder:
   > ./Thanos-United examples/basicLoadBalance_k1_3Vols summary.csv 20 0.25 0 outDir/ 21600 10 10

Paper
------
To read more about ILP please check our [paper](https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous and Ariel Kolikant, Technion - Israel Institute of Technology; 
		 Abhinav Duggal, DELL EMC; 
		 Sarai Sheinvald, ORT Braude College of Engineering; 
		 Gala Yadgar, Technion - Israel Institute of Technology

Email : sarielko@campus.technion.ac.il

For more information, please read 'ILP.pdf'.
