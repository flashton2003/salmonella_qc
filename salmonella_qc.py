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
    # print(sistr_results)
    sistr_results = sistr_results.loc[sistr_results['cgmlst_ST'] != 'cgmlst_ST']
    sistr_results = sistr_results[['genome', 'qc_messages', 'qc_status', 'serovar', 'serovar_antigen', 'serovar_cgmlst']]
    sistr_results = sistr_results.assign(genome = [x.strip('_contigs.fa') for x in sistr_results['genome']])
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

def read_mlst(mlst_inhandle):
    mlst_results = pd.read_csv(mlst_inhandle, sep = '\t', header = None)
    
    ## since just catted together all the outputs, there is a header for each one.
    mlst_results = mlst_results.assign(sample_name = [x.rstrip('_contigs.fa') for x in mlst_results[0]])
    # print(mlst_results)
    mlst_results = mlst_results[['sample_name', 2]]
    mlst_results = mlst_results.set_index('sample_name')
    # print(mlst_results)
    return mlst_results 

def take_index_union(contam, sistr_results, assembly_stats_results, mlst_results):
    sample_union = pd.Index.union(contam.index, sistr_results.index)
    sample_union = sample_union.union(assembly_stats_results.index)
    sample_union = sample_union.union(mlst_results.index)
    all_samples = pd.DataFrame({'samples':list(sample_union)})
    all_samples = all_samples.set_index('samples')
    return all_samples

def combine(all_samples, contam, sistr_results, assembly_stats_results, mlst_results, merged_outhandle):
    merge = pd.merge(all_samples, contam, how = 'outer', left_index = True, right_index = True)
    merge = pd.merge(merge, sistr_results, how = 'outer', left_index = True, right_index = True)
    merge = pd.merge(merge, assembly_stats_results, how = 'outer', left_index = True, right_index = True)
    merge = pd.merge(merge, mlst_results, how = 'outer', left_index = True, right_index = True)
    # print(merge)
    merge.to_csv(merged_outhandle, sep = '\t')

def main(kraken_inhandle, sistr_inhandle, assembly_stats_inhandle, mlst_inhandle, merged_outhandle):
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
    mlst_results = read_mlst(mlst_inhandle)
    all_samples = take_index_union(contam, sistr_results, assembly_stats_results, mlst_results)
    combine(all_samples, contam, sistr_results, assembly_stats_results, mlst_results, merged_outhandle)
    

root_dir = '/Users/flashton/Dropbox/GordonGroup/ben_kumwenda_genomes'
# kraken_inhandle = f'{root_dir}/kraken2/results/2020.09.30/2020.10.01.parsed_results.v2.txt'
# sistr_inhandle = f'{root_dir}/sistr/results/2020.09.25/2020.09.25.ben_k_sistr_results.tsv'
# assembly_stats_inhandle = f'{root_dir}/qc/results/2020.09.28/2020.09.28.assembly_stats.ben_k.tsv'
# merged_outhandle = f'{root_dir}/qc/results/2020.10.02/2020.10.02.ben_k_merged_qc.tsv'
kraken_inhandle = f'{root_dir}/kraken2/results/2020.10.08/2020.10.08.parsed_results_Feasy_Ent.txt'
sistr_inhandle = f'{root_dir}/sistr/results/2020.10.08/2020.10.08.feasey_ent_sistr_res.tsv'
assembly_stats_inhandle = f'{root_dir}/qc/results/2020.10.08/2020.10.08.assembly_stats_feasey_ent.tsv'
mlst_inhandle = f'{root_dir}/mlst/results/2020.10.08/2020.10.08.mlst_feasey_ent.tsv'
merged_outhandle = f'{root_dir}/qc/results/2020.10.08/2020.10.02.feasey_ent_merged_qc.tsv'



if __name__ == '__main__':
    main(kraken_inhandle, sistr_inhandle, assembly_stats_inhandle, mlst_inhandle, merged_outhandle)