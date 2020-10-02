import pandas as pd
import interpret_kraken2_output
import docopt

def run_interpret_kraken2_output(kraken_inhandle):
    kraken_results = interpret_kraken2_output.read_in_parsed_kraken_result(kraken_inhandle)
    kraken_interpretation = interpret_kraken2_output.make_output_dataframe(kraken_results)
    contam = interpret_kraken2_output.identify_contaminated_samples(kraken_results)
    # print(contam.index)
    return contam

def read_sistr(sistr_inhandle):
    # with open(sistr_inhandle) as fi:
    sistr_results = pd.read_csv(sistr_inhandle, sep = '\t')
    sistr_results = sistr_results[['genome', 'qc_messages', 'qc_status', 'serovar', 'serovar_antigen', 'serovar_cgmlst']]
    # print(sistr_results[['genome']])
    sistr_results = sistr_results.set_index('genome')
    return sistr_results

def read_assembly_stats(assembly_stats_inhandle):
    assembly_stats_results = pd.read_csv(assembly_stats_inhandle, sep = '\t')
    ## since just catted together all the outputs, there is a header for each one.
    assembly_stats_results = assembly_stats_results.loc[assembly_stats_results['filename'] != 'filename']
    assembly_stats_results = assembly_stats_results.assign(sample_name = [x.split('/')[-1].rstrip('_contigs.fa') for x in assembly_stats_results['filename']])
    assembly_stats_results = assembly_stats_results[['sample_name', 'total_length', 'number', 'N50']]
    assembly_stats_results = assembly_stats_results.set_index('sample_name')
    return assembly_stats_results

def make_master_df(contam_samples, sistr_samples, assembly_stats_samples):
    # all_sample_names = list(set(contam_samples).union(set(sistr_samples), set(assembly_stats_samples)))
    # print(len(all_samples))
    # print(all_samples)
    ## not sure how genome snuck in there, from sistr results probably
    if 'genome' in all_sample_names:
        all_sample_names.remove('genome')
    all_samples = pd.DataFrame({'samples':all_sample_names})
    # print(all_samples)
    return all_samples

def main(kraken_inhandle, sistr_inhandle, assembly_stats_inhandle):
    '''
    1. read in
        a. sistr results
        b. assembly stats results
        c. kraken2 parsed results
    2. filter/pre-process a, b & c
    3. combine a, b, c
    4. make decision about each one
    '''
    contam = run_interpret_kraken2_output(kraken_inhandle)
    sistr_results = read_sistr(sistr_inhandle)
    assembly_stats_results = read_assembly_stats(assembly_stats_inhandle)
    # print(contam[[0]])
    # all_samples = make_master_df(list(contam.index), sistr_results[['genome']], assembly_stats_results['sample_name'])
    # print(contam)



root_dir = '/Users/flashton/Dropbox/GordonGroup/ben_kumwenda_genomes'
kraken_inhandle = f'{root_dir}/kraken2/results/2020.09.30/2020.10.01.parsed_results.v2.txt'
sistr_inhandle = f'{root_dir}/sistr/results/2020.09.25/2020.09.25.ben_k_sistr_results.tsv'
assembly_stats_inhandle = f'{root_dir}/qc/results/2020.09.28/2020.09.28.assembly_stats.ben_k.tsv'

if __name__ == '__main__':
    main(kraken_inhandle, sistr_inhandle, assembly_stats_inhandle)