#!/usr/bin/env python

import os, numpy, pathlib
import pandas
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
import scipy.stats as stats
import pickle

from sklearn.metrics import roc_auc_score

from lifelines import CoxPHFitter
from lifelines.fitters.kaplan_meier_fitter import KaplanMeierFitter

from matplotlib.backends.backend_pdf import PdfPages

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
output_path = os.path.join(data_path, 'output')

font_size = 30
figure_width = 7

plt.rcParams.update({'font.size': font_size})

# plot functions
def boxplot_one(handle, arr, i, col, flag_dot=True, marker='o'):
    bp = handle.boxplot(arr, widths=0.6, showfliers=False, patch_artist=True, positions=[i])
    
    for median in bp['medians']: median.set(color=col, linewidth=5)
    
    for box in bp['boxes']:
        box.set_facecolor((0,0,0,0))
        #box.set_alpha(alpha)
    
    if flag_dot:
        x = numpy.random.normal(i, 0.1, size=arr.shape[0])
        handle.plot(x, arr, color=col, marker=marker, linestyle='none', markersize=5)


def ROC_plot(data_lst, output=None, flag_sort=False):
    fig = plt.figure(figsize=(figure_width, figure_width), frameon=False)
    
    sort_lst = []
    
    for title, flag, arr in data_lst:
        fpr, tpr, _ = metrics.roc_curve(flag, arr)
        sort_lst.append([title, metrics.auc(fpr, tpr), fpr, tpr])
    
    # rank by AUC in descending order
    if flag_sort: sort_lst.sort(key = lambda v: v[1], reverse=True)
    
    for title, AUC, fpr, tpr in sort_lst:
        plt.plot(fpr, tpr, label= '%s %.2f' % (title, AUC))
    
    plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    
    #plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tick_params(pad=10)
    plt.legend(frameon=False, fontsize=font_size)
    
    if output is not None:
        fig.savefig(output + '.pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        return fig




def main():
    # create a prediction signature with Tres and Tpersistance
    merge = []
    
    signature = pandas.read_csv(os.path.join(output_path, 'merge.signature'), sep='\t', index_col=0)['Tres']
    merge.append(signature)

    signature = pandas.read_excel(os.path.join(data_path, 'signature', 'Tpersistance.Krishna2020.xlsx'), index_col=0)
    signature = signature.loc[:, 'logFC']
    signature.name = 'T Persistance'
    merge.append(signature)
    
    signature = pandas.read_csv(os.path.join(data_path, 'signature', 'GSE23321.diff.gz'), sep='\t', index_col=0)['Tscm']
    merge.append(signature)
    
    signature = pandas.concat(merge, axis=1, join='outer')
    
    fin = open(os.path.join(data_path, 'evaluation.pickle.gz'), 'rb') 
    data_lst = pickle.load(fin)
    
    for treatment, timepoint, cohort, cancer, response, data, compare_lst in data_lst:
        title = '%s_%s_%s_%s' % (treatment, timepoint, cohort, cancer)
        print(title)
        
        output = os.path.join(output_path, title)
        
        correlation = signature.apply(lambda v: data.corrwith(v))
        
        common = correlation.index.intersection(response.index)
        response = response.loc[common]
        correlation = correlation.loc[common]
        
        if type(response) == pandas.Series:
            fig = plt.figure(figsize = (0.6*figure_width, figure_width), frameon=False)
            
            arr_in = correlation.loc[~response, 'Tres']
            arr_out = correlation.loc[response, 'Tres']
            
            boxplot_one(plt, arr_in, 0, 'blue', flag_dot = True, marker='o')
            boxplot_one(plt, arr_out, 1, 'red', flag_dot = True, marker='v')
            
            plt.ylabel('Correlation')
            
            if treatment.find('COVID19') == 0:
                labels = ['Severe', 'Mild']
            else:
                labels = ['Non-responder', 'Responder']
            
            plt.xticks([0, 1], labels, rotation=30, ha='right')
            
            plt.tick_params(pad=10)
            plt.axhline(0, ls='--', color='grey')
            
            z, p = stats.ranksums(arr_in, arr_out)
            print('z =',z, 'P =',p)    
            plt.title('P = %.1e' % p, fontsize=font_size)
            
            fig.savefig(output + '.boxplot.pdf', bbox_inches='tight', transparent=True)
            plt.close(fig)
            
            plot_lst = []
            
            for title, arr in correlation.items():
                plot_lst.append([title, response, arr])
            
            ROC_plot(plot_lst, output + '.ROC')
            
        else:
            out = output + '.%s' % response.columns[0].split()[0]
            
            # part 1: Cox-PH Wald test for all metrics
            cf = CoxPHFitter(penalizer=1e-3)
            
            metrics = []
            
            for title, arr in correlation.items():
                data = pandas.concat([response, arr], axis=1, join='inner')
                cf.fit(data, data.columns[0], event_col=data.columns[1])
                metrics.append(cf.summary.loc[arr.name])
            
            metrics = pandas.concat(metrics, axis=1, join='inner')    
            
            fig = plt.figure(figsize = (figure_width, figure_width), frameon=False)
            plt.bar(metrics.columns, metrics.loc['z'])
            plt.ylabel('Risk z-score')
            plt.xticks(rotation=30, ha='right')
            fig.savefig(out + '.bar.pdf', bbox_inches='tight', transparent=True)
            plt.close(fig)
            
            z, p = metrics.loc['z', 'Tres'], metrics.loc['p', 'Tres']
            print('z =',z, 'P =',p)
            
            # part 2: KM plot for Tres
            kmf = KaplanMeierFitter()
    
            fig = plt.figure(figsize = (figure_width, figure_width), frameon=False)
            
            flag = correlation.loc[:, 'Tres'] > 0
            
            kmf.fit(response.iloc[:,0].loc[flag], response.iloc[:,1].loc[flag], label= 'Tres > 0 (n=%d)' % sum(flag))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='-')
        
            kmf.fit(response.iloc[:,0].loc[~flag], response.iloc[:,1].loc[~flag], label= 'Tres < 0 (n=%d)' % sum(~flag))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='--')
            
            plt.ylabel('Fraction')
            plt.xlabel(response.columns[0])
            plt.title('P = %.1e' % p, fontsize=font_size)
            plt.legend(frameon=False)
            
            fig.savefig(out + '.kmplot.pdf', bbox_inches='tight', transparent=True)
            plt.close(fig)
    
    fin.close()
    
    return 0

if __name__ == '__main__': main()
