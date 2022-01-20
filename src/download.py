#!/usr/bin/env python

import os, pathlib, tarfile
import urllib.request

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
src_path = os.path.join(base_path, 'src')

for url in ['https://hpc.nih.gov/~Jiang_Lab/Tres/sc_cohorts.tar.gz', 'https://tres.data2intelligence.net/download/plots/evaluation.pickle.gz/']:
    f = os.path.basename(url.rstrip('/'))
    
    out = os.path.join(data_path, f)
    urllib.request.urlretrieve(url, out)
    
    if url.find('.tar.gz') > 0:
        my_tar = tarfile.open(out)
        my_tar.extractall(data_path)
        my_tar.close()
        
        os.remove(out)
