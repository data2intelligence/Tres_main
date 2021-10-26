#!/usr/bin/env python

import os, numpy, pathlib
import pandas
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
import scipy.stats as stats

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



############################################################
# one function for each dataset
def predict_ICB_Caushi2021_response(signature):
    output = os.path.join(output_path, 'Validation.ICB_Caushi2021')
    
    expression = pandas.read_csv(os.path.join(data_path, 'validation', 'GSE176021.tumor_MANA.gz'), sep='\t', index_col=0)
    expression = expression.subtract(expression.mean(axis=1), axis=0)
    
    flag = pandas.Series([v.split('.')[0] == 'MPR' for v in expression.columns], index=expression.columns, name='Response')
    correlation = signature.apply(lambda v: expression.corrwith(v))
    
    data_lst = []
        
    fig = plt.figure(figsize = (0.5*figure_width, figure_width), frameon=False)
    
    for title, arr in correlation.items():
        in_arr, out_arr = arr.loc[flag], arr.loc[~flag]    
        z, p = stats.ranksums(in_arr, out_arr)
        auc = roc_auc_score(flag, arr)
            
        incnt, outcnt = sum(flag), sum(~flag)
        
        print(title, auc, z, p, incnt, outcnt)
        
        data_lst.append([title, flag, arr])
            
        if title == 'Tres':
            boxplot_one(plt, in_arr, 0, 'red', marker='o')
            boxplot_one(plt, out_arr, 1, 'blue', marker='v')
                
            plt.text(0, in_arr.min(), incnt, horizontalalignment='center', fontsize=30)
            plt.text(1, out_arr.min(), outcnt, horizontalalignment='center', fontsize=30)
                    
            plt.title('P = %.1e' % p, fontsize=font_size)
    
    plt.xticks([0, 1], ['Responder', 'Non'], rotation=30, ha='right')
    
    plt.ylabel('Correlation')
    plt.tick_params(pad=10)
        
    fig.savefig(output + '.boxplot.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)
        
    ROC_plot(data_lst, output + '.ROC')



def predict_ICB_SadeFeldman2018_response(signature):
    output = os.path.join(output_path, 'Validation.ICB_SadeFeldman2018')
    
    expression = os.path.join(data_path, 'validation', 'Melanoma.GSE120575.merge.pickle.gz')
    expression = pandas.read_pickle(expression)
    expression = expression.subtract(expression.mean(axis=1), axis=0)
    
    included = ['Exhausted T_CD8', 'Lymphocyte exhausted/cell cycle']
    flag = [('Exhausted.' + v.split('.',1)[1] if v.split('.')[0] in included else v) for v in expression.columns]
    expression = expression.groupby(flag, axis=1).median()
    
    correlation = signature.apply(lambda v: expression.corrwith(v))
    
    flag = [v.split('.')[0] + '.' + v.split('.')[1].split('_')[0] for v in correlation.index]
    correlation_group = correlation.groupby(flag)
    
    for cell_type in ['Exhausted.Pre', 'Exhausted.Post']:
        fig = plt.figure(figsize = (0.5*figure_width, figure_width), frameon=False)
            
        correlation = correlation_group.get_group(cell_type)
        flag = pandas.Series([v.split('-')[-1] == 'R' for v in correlation.index], index=correlation.index, name='Response')
            
        data_lst = []
            
        for title, arr in correlation.items():
            in_arr, out_arr = arr.loc[flag], arr.loc[~flag]
            z, p = stats.ranksums(in_arr, out_arr)
            auc = roc_auc_score(flag, arr)
            incnt, outcnt = sum(flag), sum(~flag)
        
            print(cell_type, title, auc, z, p, incnt, outcnt)
            data_lst.append([title, flag, arr])
                
            if title == 'Tres':
                boxplot_one(plt, in_arr, 0, 'red', marker='o')
                boxplot_one(plt, out_arr, 1, 'blue', marker='v')
                    
                plt.text(0, in_arr.max(), incnt, horizontalalignment='center', fontsize=30)
                plt.text(1, out_arr.max(), outcnt, horizontalalignment='center', fontsize=30)
                
                plt.title('P = %.1e' % p, fontsize=font_size)

        plt.xticks([0, 1], ['Responder', 'Non'], rotation=30, ha='right')
        plt.ylabel('Correlation')
        plt.tick_params(pad=10)
            
        out = output + '.' + cell_type
        fig.savefig(out + '.boxplot.pdf', bbox_inches='tight', transparent=True)
        plt.close(fig)
            
        ROC_plot(data_lst, out + '.ROC')



def predict_CD19CAR_Fraietta2018_response(signature):
    fprefix = os.path.join(data_path, 'validation', 'CD19CAR_Fraietta2018')
    output = os.path.join(output_path, 'Validation.CD19CAR_Fraietta2018')
    
    responder = pandas.read_csv(fprefix + '.responders', sep='\t', index_col=0)
    responder = responder.index[responder['Best Overall Response'].apply(lambda v: v in ['CR', 'PRTD'])]
    
    # Infusion product transcriptomics
    expression = pandas.read_csv(fprefix + '.CAR.self_subtract.gz', sep='\t', index_col=0)
    correlation = signature.apply(lambda v: expression.corrwith(v))
        
    flag = pandas.Series([v in responder for v in correlation.index], index=correlation.index, name='Responder')
        
    data_lst = []
        
    for title, arr in correlation.items():
        in_arr, out_arr = arr.loc[flag], arr.loc[~flag]    
        z, p = stats.ranksums(in_arr, out_arr)
        auc = roc_auc_score(flag, arr)
            
        incnt, outcnt = sum(flag), sum(~flag)
        print(title, auc, z, p, incnt, outcnt)
            
        data_lst.append([title, flag, arr])
            
        if title == 'Tres':
            fig = plt.figure(figsize = (0.5*figure_width, figure_width), frameon=False)
            boxplot_one(plt, in_arr, 0, 'red', marker='o')
            boxplot_one(plt, out_arr, 1, 'blue', marker='v')
                
            plt.text(0, in_arr.min(), incnt, horizontalalignment='center', fontsize=30)
            plt.text(1, out_arr.min(), outcnt, horizontalalignment='center', fontsize=30)
            
            plt.title('P = %.1e' % p, fontsize=font_size)
            plt.xticks([0, 1], ['Responder', 'Non'], rotation=30, ha='right')
            plt.ylabel('Correlation')
            plt.tick_params(pad=10)
    
            fig.savefig(output + '.boxplot.pdf', bbox_inches='tight', transparent=True)
            plt.close(fig)
                
    ROC_plot(data_lst, output + '.ROC')




def predict_CD19CAR_Chen2021_response(signature):
    margin = 0.05
    metric_type = 'Tres'
    
    fprefix = os.path.join(data_path, 'validation', 'CD19CAR_Chen2021')
    output = os.path.join(output_path, 'Validation.CD19CAR_Chen2021')
    
    # Predict CAR-T response through pre-manufecture    
    clinical = pandas.read_csv(fprefix + '.clinical', sep='\t', index_col=0)
    data = pandas.read_csv(fprefix + '.self_subtract.gz', sep='\t', index_col=0)
    data_groups = data.groupby([v.split('.')[-1] for v in data.columns], axis=1)
    
    # always add some penalty in case perfect separation    
    cf = CoxPHFitter(penalizer=1e-3)
    kmf = KaplanMeierFitter()
    
    pdf = PdfPages(output + '.example.pdf')
        
    for cell_type, data in data_groups:
        data.columns = [v.split('.')[0] for v in data.columns]
        correlation = signature.apply(lambda v: data.corrwith(v))
        
        arr = correlation[metric_type]
        
        data = pandas.concat([clinical[['BCA_MONTHS', 'BCELL_RECOVERY', 'AGE']], arr], axis=1, join='inner')
        cf.fit(data, data.columns[0], event_col=data.columns[1])
                
        # plot here
        if cell_type == 'Average':        
            fig = plt.figure(figsize = (figure_width, figure_width), frameon=False)
            
            flag = (data.iloc[:,-1] > 0)
            kmf.fit(data.iloc[:,0].loc[flag], data.iloc[:,1].loc[flag], label='%s > 0 (n=%d)' % (metric_type.split()[0], sum(flag)))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='solid', color='red')
                    
            kmf.fit(data.iloc[:,0].loc[~flag], data.iloc[:,1].loc[~flag], label='%s < 0 (n=%d)' % (metric_type.split()[0], sum(~flag)))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='dashdot', color='blue')
    
            p = cf.summary.loc[metric_type, 'p']
            p = p/2 # convert two-sided to one-sided p-values
            plt.text(data.iloc[:,0].median(), 0.9, 'P = %.1e' % p)
                    
            plt.title(cell_type, fontsize=font_size)
            plt.xlabel('Bcell aplasia survival (month)')
            plt.ylabel('Fraction')
            plt.legend(frameon=False)
                    
            plt.tick_params(pad=10)
            plt.ylim([0 - margin,1 + margin])
                    
            pdf.savefig(fig, bbox_inches='tight', transparent=True)
            plt.close(fig)
    
    pdf.close()




def predict_ACT_Lauss2017_response(signature):
    margin = 0.05
    metric_type = 'Tres'
    
    fprefix = os.path.join(data_path, 'validation', 'Lauss2017')
    output = os.path.join(output_path, 'Validation.ACT_Lauss2017')
    
    expression = pandas.read_csv(fprefix + '.self_subtract.gz', sep='\t', index_col=0)
    correlation = signature.apply(lambda v: expression.corrwith(v))
    
    CTL = expression.loc[['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1']].median()
    CTL.name = 'CTL'
    
    # add a little penalty to fix the perfection separation here
    cf = CoxPHFitter(penalizer=1e-3)
    kmf = KaplanMeierFitter()
    

    pdf = PdfPages(output + '.pdf')
    
    for surv_type in ['OS', 'PFS']:
        clinical = pandas.read_csv(fprefix + '.' + surv_type, sep='\t', index_col=0)
    
        arr = correlation[metric_type]
                
        # separate CTL level analysis
        for CTL_type, CTL_flag in [
            ['CTL > 0', (CTL > 0)],
            ['CTL < 0', (CTL < 0)],
            ]:
                
            # only look at CTL rich samples
            data = pandas.concat([clinical, arr.loc[CTL_flag]], axis=1, join='inner')
            cf.fit(data, data.columns[0], event_col=data.columns[1])
            
            
            fig = plt.figure(figsize = (figure_width, figure_width), frameon=False)
                    
            flag = (data.iloc[:,-1] > 0)
            kmf.fit(data.iloc[:,0].loc[flag], data.iloc[:,1].loc[flag], label='%s > 0 (n=%d)' % (metric_type.split()[0], sum(flag)))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='solid', color='red')

            kmf.fit(data.iloc[:,0].loc[~flag], data.iloc[:,1].loc[~flag], label='%s < 0 (n=%d)' % (metric_type.split()[0], sum(~flag)))
            kmf.plot(ci_show=False, show_censors=True, linewidth=2, ls='dashdot', color='blue')
                
            p = cf.summary.loc[metric_type, 'p']
            p /= 2 # convert two-sided p-value to one-side pvalue
            
            plt.title('%s : P = %.1e' % (CTL_type, p), fontsize=font_size)
                    
            if surv_type == 'OS':
                plt.xlabel('Overall survival (month)')
            elif surv_type == 'PFS':
                plt.xlabel('Prog-free survival (month)')
            else:
                plt.xlabel(surv_type + ' (month)')
                    
            plt.ylabel('Fraction')
            plt.legend(frameon=False)
            
            plt.tick_params(pad=10)
            plt.ylim([0 - margin,1 + margin])
                    
            pdf.savefig(fig, bbox_inches='tight', transparent=True)
            plt.close(fig)
    
    pdf.close()



def main():
    # create a prediction signature with Tres and Tpersistance
    merge = []
    
    signature = pandas.read_csv(os.path.join(output_path, 'merge.Median.signature'), sep='\t', index_col=0)['Tres']
    merge.append(signature)

    signature = pandas.read_excel(os.path.join(data_path, 'signature', 'Tpersistance.Krishna2020.xlsx'), index_col=0)
    signature = signature.loc[:, 'logFC']
    signature.name = 'T Persistance'
    merge.append(signature)    

    signature = pandas.concat(merge, axis=1, join='outer')
    
    # validation and performance comparison by ROC curves
    predict_ICB_Caushi2021_response(signature)
    predict_ICB_SadeFeldman2018_response(signature)
    predict_CD19CAR_Fraietta2018_response(signature)
    
    # survival analysis using pre-manufacture data
    predict_CD19CAR_Chen2021_response(signature)
    predict_ACT_Lauss2017_response(signature)
    
    return 0

if __name__ == '__main__': main()
