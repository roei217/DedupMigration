# Calculator

General Information
-------------------
Calculator results a migration plan's cost for deduplicated storage systems.

Dependencies
------------
1. sqlite3 is required to use sqlite3 database.

Running
-------

1. Complie the code
   > make

2. Run calculator instance
   > calc -file file_path *-output output_name* -cache_path db_path -no_cache -lb_sizes sizes_list
		
		2.1. file - the migration plan's file
		2.2. cache_path - the full path to a sqlite3 database (.db)
		Optional parameters:
			2.3. no_cache - do not use the cache (useful for small volumes).
			2.4. lb_sizes -  a list of the requested sizes for each final volume, default is set to perfect load balancing
			2.5. output - the name of the output file, default is random

Note that if more than 20 volumes are used, the cache will not be used regardless of the -no_cache option.

Examples
------------------
Calculate a migration plan's cost for the plan in *example/example_migration_sample.csv* with cache path *example/example.db*
and the output file is *example/output_sampled.csv*.
   > ./calc -file example/example_migration_sample.csv -output example/output_sampled.csv -cache_path example/example.db
  
Calculate a migration plan's cost for the plan in *example/example_migration_original.csv* with cache path *example/example.db* but do not use it
and the output file is *example/output_original.csv*.
   > ./calc -file example/example_migration_original.csv -output example/output_original.csv -cache_path example/example.db -no_cache
  
Paper
------
To read more about Clustering please check our *paper*(https://www.usenix.org/conference/fast22/presentation/kisous):
   "The what, The from, and The to: The Migration Games in Deduplicated Systems". In 20th USENIX Conference on File and Storage Technologies (FAST 22), 2022.

Authors: Roei Kisous and Ariel Kolikant, Technion - Israel Institute of Technology; 
		 Abhinav Duggal, DELL EMC; 
		 Sarai Sheinvald, ORT Braude College of Engineering; 
		 Gala Yadgar, Technion - Israel Institute of Technology

Email : kisous[dot]roei[at]campus[dot]technion[dot]ac[dot]il

For more information, please read 'CostCalc.pdf'.