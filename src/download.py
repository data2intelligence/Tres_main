#!/usr/bin/env python

import os, pathlib, tarfile
import urllib.request

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
src_path = os.path.join(base_path, 'src')

url = 'https://hpc.nih.gov/~Jiang_Lab/Tres'

for f in ['sc_cohorts.tar', 'validation.tar']:
    out = os.path.join(data_path, f + '.gz')
    
    urllib.request.urlretrieve(os.path.join(url, f + '.gz'), out)

    my_tar = tarfile.open(out)
    my_tar.extractall(data_path)
    my_tar.close()
    
    os.remove(out)
