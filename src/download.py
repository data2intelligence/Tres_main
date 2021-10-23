#!/usr/bin/env python

import os, pathlib, tarfile
import urllib.request

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
src_path = os.path.join(base_path, 'src')

url = 'https://hpc.nih.gov/~Jiang_Lab/Tres'

for f in ['sc_cohorts.tar']:
    urllib.request.urlretrieve(os.path.join(url, f + '.gz'), os.path.join(data_path, f + '.gz'))

# expand TCGA and GTEx cohorts
f = os.path.join(data_path, 'sc_cohorts.tar.gz')

my_tar = tarfile.open(f)
my_tar.extractall(data_path)
my_tar.close()

os.remove(f)
