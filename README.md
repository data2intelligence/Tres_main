Example code of run Tres analysis and reproduce major results  

# Stage 1: download single-cell datasets  
Please CD into the src folder and run "./download.py" first to download all data files. There are 16 single-cell datasets.  
  
So, you may run "./run.py inx 16" where inx is a number between 0 and 16 to compute the correlation results for each dataset.  
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 16". You may re-write this file for your local HPC.  

# Stage 2: merge computation scores into signature files  
After computing the correlation results for each dataset, you can run "./run.py" to generate the merged Tres signature.

# Stage 3: predict clinical response of immune checkpoint blockade, CAR T, and adoptive cell transfer  
After creating the Tres signature, you can run "./predict.py" to generate therapy response prediction results.

