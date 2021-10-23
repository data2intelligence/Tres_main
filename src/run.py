#!/usr/bin/env python

import os, sys, pathlib

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')


# Cancer Set
dataset_Tumor_CD8 = [
    ['Glioblastoma.GSE163108.10x', 'CD8+_recurrent', None],
    
    ['Uveal.GSE139829', 'Cytotoxic CD8', None],
    
    ['Colorectal.GSE108989', 'TTC', None],
    
    ['Colorectal.GSE146771', 'T_CD8 T cell', None],
    ['Colorectal.GSE146771.10x', 'T_CD8 T cell', None],
    
    ['Head_Neck.GSE103322', 'Tumor_T cell', None],
    
    ['Liver.GSE125449', 'T cell', None],
    
    ['Liver.GSE98638', 'TTC', None],
       
    ['Liver.GSE140228', 'Tumor_CD8', None],
    ['Liver.GSE140228.10x', 'Tumor_CD8', None],
        
    ['NSCLC.GSE99254', 'TTC', None],
    ['Breast.GSE114725', 'TUMOR_T:CD8+EM', None],
]

for location in ['Core', 'Middle', 'Edge']:
    dataset_Tumor_CD8.append(['NSCLC.E-MTAB-6149', '%s_T CD8' % location, location])

for ttime in ['Naive', 'Post']:
    dataset_Tumor_CD8.append(['Melanoma.GSE115978', '%s_T_CD8' % ttime, ttime])

for sub_type in ['ER+', 'HER2+', 'TNBC']:
    dataset_Tumor_CD8.append(['Breast.GSE176078', '%s_T cells CD8+' % sub_type, sub_type])

for sub_type in ['ABC', 'GCB']:
    dataset_Tumor_CD8.append(['DLBCL.GSE182434', 'DLBCL-%s_T cells CD8' % sub_type, sub_type])



# custom single cell cohorts
datasets = [    
    # small datasets first
    'Melanoma.GSE115978',
    
    'Head_Neck.GSE103322',
    'Glioblastoma.GSE131928',
    
    'Colorectal.GSE108989',
    
    'Liver.GSE98638',
    'Liver.GSE140228',
    'Colorectal.GSE146771', 
    
    'NSCLC.GSE99254',
    'Breast.GSE114725',
    'Liver.GSE125449', # 10x small dataset
    
    ################################################
    # 10x region, lots of cells
    'Breast.GSE176078',
    'Glioblastoma.GSE163108.10x',
    'Liver.GSE140228.10x',
    'Colorectal.GSE146771.10x',
    'NSCLC.E-MTAB-6149',
    'Uveal.GSE139829',
    
    'DLBCL.GSE182434',
]

def main():
    input_path = os.path.join(data_path, 'sc_cohorts')
    output_path = os.path.join(data_path, 'output')
    
    if len(sys.argv) == 3:    
        inx, Nnode = int(sys.argv[1]), int(sys.argv[2])
        assert Nnode == len(datasets)
        
        dataset = datasets[inx]
        
        output = os.path.join(output_path, dataset)
        os.system(' '.join(['Tres.py -n 1', '-i', os.path.join(input_path, dataset + '.pickle.gz'), '-o', output]))
    
    elif len(sys.argv) == 1:
        print('merge')
    
    else:
        sys.stderr.write('Error parameter input.\n')
        return 1
    
    return 0

if __name__ == '__main__': main()
