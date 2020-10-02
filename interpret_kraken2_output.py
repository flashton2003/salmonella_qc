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
    kraken_results = pd.read_csv(kraken_inhandle, sep = '\t', header = None)
    return kraken_results

def make_output_dataframe(kraken_results):
    kraken_interpretation = pd.DataFrame()
    kraken_interpretation['sample_name'] = kraken_results[0].unique().tolist()
    return kraken_interpretation

def get_count_column(contam):
    output_list = []
    counting_dict = {}
    for name in contam[0]:
        if name in counting_dict:
            counting_dict[name] += 1
            output_list.append(counting_dict[name])
        else:
            counting_dict[name] = 1
            output_list.append(counting_dict[name])
    return output_list

def identify_contaminated_samples(kraken_results):
    contam = kraken_results.loc[kraken_results[6] != 'Salmonella enterica']
    contam = contam.assign(contam_string = [f'{x}; {y}' for x, y in zip(contam[6], contam[1])])
    contam = contam[[0, 'contam_string']]
    ## need to add a column index, in addition to the "row" index which will be the sample name
    count_column = get_count_column(contam)
    contam = contam.assign(count_column = count_column)
    contam = contam.pivot(index = 0, columns = 'count_column', values = 'contam_string')
    return contam

def main(kraken_inhandle):
    kraken_results = read_in_parsed_kraken_result(kraken_inhandle)
    kraken_interpretation = make_output_dataframe(kraken_results)
    contam = identify_contaminated_samples(kraken_results)
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        # print(contam)

root_dir = '/Users/flashton/Dropbox/GordonGroup/ben_kumwenda_genomes'
kraken_inhandle = f'{root_dir}/kraken2/results/2020.09.30/2020.10.01.parsed_results.v2.txt'

if __name__ == '__main__':
    main(kraken_inhandle)

