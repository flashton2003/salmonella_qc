import string
import pandas as pd

class KrakenResult():

    def __init__(self, split_line):
        self.sample_name = split_line[0]
        self.percent_reads_assigned = float(split_line[1])
        self.number_reads_rooted_here = int(split_line[2])
        self.number_reads_assigned_here = int(split_line[3])
        self.taxonomic = split_line[4]
        self.rank = self.full_rank.rstrip(string.digits)
        self.ncbi_id = int(split_line[5])
        self.name = split_line[6].strip()

def read_in_parsed_kraken_result(kraken_inhandle):
    # output_dict = {}
    # with open(kraken_inhandle) as fi:
        # for line in fi.readlines():
            # print(line.strip())
            # split_line = line.strip().split('\t')
            # kraken_result = KrakenResult(split_line)
    kraken_results = pd.read_csv(kraken_inhandle, sep = '\t', header = None)
    


def main(kraken_inhandle):
    read_in_parsed_kraken_result(kraken_inhandle)

root_dir = '/Users/flashton/Dropbox/GordonGroup/ben_kumwenda_genomes'
kraken_inhandle = f'{root_dir}/kraken2/results/2020.09.30/2020.10.01.parsed_results.txt'

if __name__ == '__main__':
    main(kraken_inhandle)