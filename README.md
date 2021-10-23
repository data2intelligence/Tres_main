Example code of run Tres analysis and reproduce major results  
  
Please CD into the src folder and run "./download.py" first to download all data files. There are 17 single-cell datasets.  
  
So, you may run "./run.py inx 17" where inx is a number between 0 and 17 to compute the correlation results for each dataset.  
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 17". You may re-write this file for your local HPC.  

After computing the correlation results for each dataset, you can run "./run.py" to generate the merged Tres signature.  