#!/usr/bin/env python
import os, sys, pathlib, pandas
from sklearn.metrics import roc_auc_score

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)
data_path = os.path.join(base_path, 'data')
input_path = os.path.join(data_path, 'sc_cohorts')
output_path = os.path.join(data_path, 'output')
    

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



def median_merge(output_path):
    output = os.path.join(output_path, 'merge')

    lst = []
    
    pivots = ['TGFB1', 'TRAIL', 'PGE2']
        
    for pivot in pivots:
        merge = []
        
        # extract CD8 T cells
        for dataset, cell_pivot, sub_title in dataset_Tumor_CD8:
            result = pandas.read_csv(os.path.join(output_path, dataset), sep='\t', index_col=0)
            
            flag = [(v.find('t.' + cell_pivot) == 0 and v.split('.')[-1] == pivot) for v in result.columns]
            assert sum(flag) > 0
            
            result = result.loc[:, flag]
            assert result.columns.value_counts().max() == 1 # if sample ID extraction is correct, must have no redundancy
            
            if sub_title is not None: dataset += ('.' + sub_title)        
            
            result = result.loc[result.isnull().mean(axis=1) < 0.5]
            result = result.median(axis=1)
            result.name = dataset
                
            merge.append(result)
        
        merge = pandas.concat(merge, axis=1, join='outer', sort=False)
        assert merge.columns.value_counts().max() == 1
        
        merge.to_csv(output + '.' + pivot, sep='\t', index_label=False)
        lst.append(merge)
    
    # create median signature across three immuno-suppressive signals    
    merge = pandas.concat(lst, axis=1, join='inner')
    
    cnt_map = merge.columns.value_counts()
    assert sum(cnt_map != len(pivots)) == 0
    
    merge = merge.groupby(merge.columns, axis=1).median()
    merge.to_csv(output + '.Median', sep='\t', index_label=False)




def ROC_AUC_set(arr, logFC):
    common = arr.index.intersection(logFC.index)
    arr = arr.loc[common]
    logFC = logFC.loc[common]
        
    return roc_auc_score(logFC > 0, arr)


def compute_signature_AUC():
    qthres = 0.05
    logFC_thres = 0.5
    null_thres = 0.5
    AUC_thres = 0.7
    
    signature = pandas.read_excel(os.path.join(data_path, 'signature', 'Tpersistance.Krishna2020.xlsx'), index_col=0)
    signature = signature.loc[signature['FDR'] < qthres, 'logFC']
    signature = signature.loc[signature.abs() > logFC_thres]    
    
    print('Pos, Neg = %d, %d' % (sum(signature > 0), sum(signature < 0)))
    
    output = os.path.join(output_path, 'merge.Median')
    
    result = pandas.read_csv(output, sep='\t', index_col=0)
    result = result.loc[result.isnull().mean(axis=1) < null_thres]
    
    AUC = result.apply(lambda v: ROC_AUC_set(v.dropna(), signature)).sort_values(ascending=False)
    AUC.name = 'AUC'
    AUC.to_csv(output + '.AUC', sep='\t', index_label=False)
    
    # create a median signature across all dataset for prediction purpose later
    signature = result.loc[:, AUC.index[AUC > AUC_thres]].dropna().median(axis=1)
    signature.name = 'Tres'
    signature.to_csv(output + '.signature', sep='\t', index_label=False)




def main():
    if len(sys.argv) == 3:    
        inx, Nnode = int(sys.argv[1]), int(sys.argv[2])
        assert Nnode == len(datasets)
        
        dataset = datasets[inx]
        
        output = os.path.join(output_path, dataset)
        os.system(' '.join(['Tres.py -n 1', '-i', os.path.join(input_path, dataset + '.pickle.gz'), '-o', output]))
    
    elif len(sys.argv) == 1:
        median_merge(output_path)
        compute_signature_AUC()
    else:
        sys.stderr.write('Error parameter input.\n')
        return 1
    
    return 0

if __name__ == '__main__': main()
