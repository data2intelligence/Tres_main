#!/usr/bin/env python
import os, sys, pathlib, pandas
from sklearn.metrics import roc_auc_score

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)
data_path = os.path.join(base_path, 'data')
input_path = os.path.join(data_path, 'sc_cohorts')
output_path = os.path.join(data_path, 'output')
    

# Cancer Set
dataset_Tumor_CD8 = []

# custom single cell cohorts
datasets = []

platform_map = {}
cancer_map = {}


def load_single_cell_datasets():
    """
    Fill in dataset information from the catalog
    """
    
    catalog = pandas.ExcelFile(os.path.join(data_path, 'sc_cohorts', 'catalog.xlsx'), engine='openpyxl')
    
    for category in catalog.sheet_names:
        data = pandas.read_excel(catalog, category, index_col=0)
        data = data.loc[~data['Cancer'].isnull()] # in case any null rows
        
        for title, info in data.iterrows():            
            f = os.path.join(data_path, 'sc_cohorts', title + '.pickle.gz')
            assert os.path.exists(f)
            
            datasets.append([title, f])        
            platform_map[title] = info['Platform']
            cancer_map[title] = info['Cancer']
            
            # all na fields as non-existing
            info = info.dropna()
                    
            Subcohort = None
            if 'Subcohort' in info.index:
                Subcohort = info['Subcohort'].split(',')
                for i, v in enumerate(Subcohort):
                    v = Subcohort[i] = v.strip()
                    platform_map[title + '.' + v] = info['Platform']
            
            if 'CD8' in info.index:
                pivot = info['CD8']
                
                if Subcohort is None:
                    dataset_Tumor_CD8.append([title, pivot, None])
                else:
                    for v in Subcohort: dataset_Tumor_CD8.append([title, v + '_' + pivot, v])
    
    catalog.close()


def strip_cancer_type_list(lst):
    return [v.split('_')[0].split('.')[0] if v.find('Liver') < 0 else v.split('_')[0].split('.', 1)[1] for v in lst]



def median_merge(output_path):
    qthres = 0.05
    frac_thres = 1e-3
    
    output = os.path.join(output_path, 'merge')
    
    lst = []
    
    pivots = ['TGFB1', 'TRAIL', 'PGE2']
        
    for pivot in pivots:
        merge = []
        
        # extract CD8 T cells
        for dataset, cell_pivot, sub_title in dataset_Tumor_CD8:
            result = pandas.read_csv(os.path.join(output_path, dataset), sep='\t', index_col=0)
            
            # focus on current pivot
            flag = [(v.split('.')[-1] == pivot) for v in result.columns]
            assert sum(flag) > 0
            
            result = result.loc[:, flag]
            result = result.loc[result.isnull().mean(axis=1) < 1]
            
            # extract t-values
            flag = [(v.find('t.' + cell_pivot) == 0) for v in result.columns]
            assert sum(flag) > 0
            
            result_t = result.loc[:, flag]
            
            # strip t and pivot
            result_t.columns = ['.'.join(v.split('.')[1:-1]) for v in result_t.columns]
            assert result_t.columns.value_counts().max() == 1 # if sample ID extraction is correct, must have no redundancy
            
            # extract q-values
            flag = [(v.find('q.' + cell_pivot) == 0) for v in result.columns]
            assert sum(flag) > 0
            
            result_q = result.loc[:, flag]
            
            # strip t and pivot
            result_q.columns = ['.'.join(v.split('.')[1:-1]) for v in result_q.columns]
            assert result_q.columns.value_counts().max() == 1 # if sample ID extraction is correct, must have no redundancy
            
            flag = (result_q < qthres).mean() > frac_thres
            
            if flag.sum() == 0:
                # nothing to include
                continue
            
            result = result_t.loc[:, flag]
            
            cancer = cancer_map[dataset]
            if sub_title is not None: cancer += ('.' + sub_title)        
            
            result.columns = cancer + '_' + dataset + '_' + result.columns
            #result = result.loc[result.isnull().mean(axis=1) < 0.5]
            #result = result.median(axis=1)
            #result.name = cancer + '_' + dataset
            
            merge.append(result)
        
        merge = pandas.concat(merge, axis=1, join='outer', sort=False)
        assert merge.columns.value_counts().max() == 1
        
        merge.to_csv(output + '.' + pivot, sep='\t', index_label=False)
        lst.append(merge)
    
    # create an average score as Tres score
    common_col = None
    
    for result in lst:
        print(result.shape)    
        # get single cells that are profiled by all signals
        if common_col is None:
            common_col = result.columns
        else:
            common_col = common_col.intersection(result.columns)
    
    print(common_col.shape)
    for i, result in enumerate(lst): lst[i] = result.loc[:, common_col]

    # create median signature across three immuno-suppressive signals    
    merge = pandas.concat(lst, axis=1, join='inner')
    
    cnt_map = merge.columns.value_counts()
    assert sum(cnt_map != len(pivots)) == 0
    
    merge = merge.groupby(merge.columns, axis=1).median()
    print(merge.shape)
    
    merge.to_csv(output + '.Median', sep='\t', index_label=False)

    # create cancer-level merge
    #merge = pandas.read_csv(output + '.Median', sep='\t', index_col=0)
    merge = merge.loc[merge.isnull().mean(axis=1) < 0.5]
    print(merge.shape)
    
    # don't merge liver cancer
    flag = strip_cancer_type_list(merge.columns)
    
    merge = merge.groupby(flag, axis=1).median()
    merge.to_csv(output, sep='\t', index_label=False)
    
    # create a median signature across all dataset for prediction purpose later
    merge = merge.dropna().median(axis=1)
    merge.name = 'Tres'
    merge.to_csv(output + '.signature', sep='\t', index_label=False)


    




def ROC_AUC_set(arr, logFC):
    common = arr.index.intersection(logFC.index)
    arr = arr.loc[common]
    logFC = logFC.loc[common]
        
    return roc_auc_score(logFC > 0, arr)


def compute_signature_AUC():
    qthres = 0.05
    logFC_thres = 0.5
    null_thres = 0.5
    
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




def main():
    load_single_cell_datasets()
    
    if len(sys.argv) == 3:
        inx, Nnode = int(sys.argv[1]), int(sys.argv[2])
        assert Nnode == len(datasets)
        
        dataset, f = datasets[inx]
        
        output = os.path.join(output_path, dataset)
        os.system(' '.join(['Tres.py -n 1', '-i', f, '-o', output]))
    
    elif len(sys.argv) == 1:
        median_merge(output_path)
        compute_signature_AUC()
    else:
        sys.stderr.write('Error parameter input.\n')
        return 1
    
    return 0

if __name__ == '__main__': main()
