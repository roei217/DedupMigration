# Clustering

General Information
-------------------
Clustering is a migration plan creator for deduplicated storage systems.

Running
-------

1. Complie the code
   > make

2. Run clustering instance
   > hc -workloads workloads_list -fps num_fps -WT WTs_list -seed seeds_list -gap gaps_list -output_path_prefix path 
		*-lb* *-lb_sizes sizes_list* *-eps eps_value* *-clusters num_clusters* 
		
		2.1. workloads_list - list of full paths to the volumes
		2.2. num_fps - number of min hash fingerprints (input 'all' for all fps)
		2.3. WTs_list - list of WTs
		2.4. seeds_list - list of seeds
		2.5. gaps_list - list of gaps
		2.6. path - relative path + prefix of the result files
		Optional parameters:
			2.7. lb - use -lb for load balancing, otherwise it will not be considered
			2.8. sizes_list -  a list of clusters' requested sizes, must sum to 100 (default is even distribution)
			2.9. eps_value - % to add to system's initial size at every iteration, default is 5%
			2.10. num_clusters - number of clusters (final volumes), default is same number of input volumes

Make sure -output_path_prefix directory actually exists.

Examples
------------------
Find a migration plan for *example/vol1.csv example/vol2.csv* using *10,000* blocks, WT *1*, 
seed *0*, gap *1%* with *load balance* and the output prefix is *results/example1*.
The other parameters are default: sizes_list is *50%, 50%*, num_clusters is *2*, eps_value is *5%*.
   > ./hc -workloads example/vol1.csv example/vol2.csv -fps 10000 -WT 1 -lb -seed 0 -gap 1 -output_path_prefix results/example1
  
Find a migration plan for *example/vol1.csv example/vol2.csv* using *all* the blocks, two different WT *0,1*, 
two different seeds *0, 10*, two different gaps *0%, 1%*, eps is *5%* with *load balance* and the output prefix is *results/example2*.
The other parameters are default: sizes_list is *50%, 50%*, num_clusters is *2*.
   > ./hc -workloads example/vol1.csv example/vol2.csv -fps all -WT 0 1 -lb -seed 0 10 -gap 0 1 -eps 5 -output_path_prefix results/example2
  
Find a migration plan for *example/vol1.csv example/vol2.csv* using *all* the blocks, two different WT *0,1*, 
two different seeds *0, 10*, two different gaps *0%, 1%*, eps is *5%*, sizes_list is *70%, 30%* with *load balance* and the output prefix is *results/example3*.
The other parameters are default: num_clusters is *2*.
   > ./hc -workloads example/vol1.csv example/vol2.csv -fps all -WT 0 1 -lb -seed 0 10 -gap 0 1 -lb_sizes 30 70 -output_path_prefix results/example3
  
  
Find a migration plan for *example/vol1.csv example/vol2.csv* using *all* the blocks, two different WT *0,1*, 
two different seeds *0, 10*, two different gaps *0%, 1%*, eps is *5%*, with *no load balance* and the output prefix is *results/example4*.
The other parameters are default: sizes_list is *50%, 50%* but irrelevant, num_clusters is *2*.
   > ./hc -workloads example/vol1.csv example/vol2.csv example/empty.csv -fps all -WT 0 1 -seed 0 10 -gap 0 1 -output_path_prefix results/example4
  
Paper
------
To read more about Clustering please check our [paper](https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous and Ariel Kolikant, Technion - Israel Institute of Technology; 
		 Abhinav Duggal, DELL EMC; 
		 Sarai Sheinvald, ORT Braude College of Engineering; 
		 Gala Yadgar, Technion - Israel Institute of Technology

Email : kisous[dot]roei[at]campus[dot]technion[dot]ac[dot]il

For more information, please read 'Clustering.pdf'.