import docopt

class Sample():


def main():
    '''
    1. read in
        a. sistr results
        b. assembly stats results
        c. kraken2 parsed results
    2. filter/pre-process a, b & c
    3. combine a, b, c
    4. make decision about each one
    '''
    pass

root_dir = '/Users/flashton/Dropbox/GordonGroup/ben_kumwenda_genomes'
kraken_inhandle = f'{root_dir}/kraken2/results/2020.09.30/2020.10.01.parsed_results.txt'
sistr_inhandle = f'{root_dir}/sistr/results/2020.09.25/2020.09.25.ben_k_sistr_results.tsv'
assembly_stats_inhandle = f'{root_dir}/qc/results/2020.09.28/2020.09.28.assembly_stats.ben_k.tsv'

if __name__ == '__main__':
    main()