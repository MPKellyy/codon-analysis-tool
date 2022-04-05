import matplotlib
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
from tqdm import tqdm
import seaborn as sns


# Function used to identify and keep track of incoming codons
# Input: sequence, codon frequency dictionary, upstream position -1 of last nucleotide in site A
# NOTE: ASSUMES DNA DATA IS 5-TO-3 PRIME CODING STRANDS
def track_codons(inputseq, site_map, codon_limit, offset):
    codonsread = 0  # Keeps track of codons read
    length = len(inputseq)  # Saving length of current sequence

    # If offset is larger than sequence size, skip
    if (offset >= length-1) or (offset+1 >= length-1) or (offset+2 >= length-1):
        return codonsread

    end = offset  # Starting index for codon frame read
    start = end + 2  # Ending index for codon frame read

    # Reversing sequence for easier read logic
    inputseq = inputseq[::-1]

    # While the frame has three nucleotides to read
    for site in range(0, codon_limit):
        # Ensures reading frame doesn't run-off DNA fragment (ie accessing 3rd codon when only 6 nucleotides were given)
        if (end >= length) or (end+1 >= length) or (start >= length):
            break

        # Saving the current codon in the reading frame
        codon = inputseq[start] + inputseq[end+1] + inputseq[end]

        # Try incrementing number of times seen of current codon in frequency dictionary
        try:
            site_map[site][codon] = site_map[site][codon] + 1
        # If an exception is thrown, current codon was not seen yet, add it to dictionary
        except:
            site_map[site][codon] = 1

        # Shift reading frame up three nucleotides
        start += 3
        end += 3

        # Update number of codons read
        codonsread += 1

    # Return number of codons read
    return codonsread


# Function used to convert codon frequency dictionary to a dataframe
# Input: Name of read, length of read, codon frequency dictionary
def create_df(readname, total_codons, site_map):
    # Codon to amino acid table
    amino_acid_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
        'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W',
    }

    # Creating an empty list to hold all entries in the dataframe
    data = []
    codon_index = 0

    sites = ["A", "P", "E"]  # Keeps track of which site is being read

    for site in range(0, len(site_map)):
        # Iterating through all codons saved in codon map (aka frequency dictionary)
        for codon in site_map[site]:
            # Creating a list to store current entry in dataframe
            amino_acid = ""
            try:
                amino_acid = amino_acid_table[codon]
            except:
                amino_acid = None

            # Setting which type read was
            type = None
            if "WT" in readname:
                type = "Wild"
            elif "K" in readname:
                type = "Mutant"

            # Cols: Name of read codon came from, codon sequence, times seen, frequency in read, contains an ambiguous read
            codondata = [readname, codon, amino_acid, site_map[site][codon], sites[site], "N" in codon, type]

            # Adding entry to dataframe list
            data.append(codondata.copy())

    # Convert dataframe list into a dataframe, removing unnecessary index column from dataframe that is set on default
    df = pd.DataFrame(data, columns=['read_name', 'codon', 'amino_acid', 'times_seen', 'site', 'contains_ambiguous', 'type'])
    # Print and return dataframe
    df = df.set_index('read_name')
    return df


def read_file(filename):
    # Dictionary used to keep track of seen codons (frequency dictionary)
    # codonmap = {}
    site_map = [{}, {}, {}]  # Index 0: A, Index 1: P, Index 2: E
    codoncounter = 0
    count = 0

    # Iterates through all entries in a fastq file and tracks number of codons seen
    with open(filename) as in_handle:
        # For every sequence in fastq file
        for title, seq, qual in tqdm(FastqGeneralIterator(in_handle), desc=filename):
            # Analyze all codons in current sequence, update number of codons seen
            codoncounter += track_codons(seq, site_map, 3, 9)

        # Creating a name for this read by removing the .fastq portion of the input file
        readname = filename[0:len(filename) - 6:]
        # SQL has a problem with file names containing '-'. Ensuring these get replaced with '_' instead
        readname = readname.replace("-", "_")

        # Creating and returning dataframe from current file
        df = create_df(readname, codoncounter, site_map)
        return df


def populate_site_dict(site_dict, site, input_df):
    for codon_key in site_dict.keys():
        temp_df = input_df[(input_df.codon == codon_key)]
        temp_df = temp_df[["codon", "site", "times_seen"]]

        #print(temp_df[(temp_df.site == site)])
        sum_seen = (temp_df[(temp_df.site == site)])["times_seen"].sum()
        site_dict[codon_key] = sum_seen
        #print(site_dict)

# Acquiring codons per site in sample type data
def get_site_occurrences(sample_df):
    sample_E_dict = {}
    sample_P_dict = {}
    sample_A_dict = {}

    for codon in sample_df["codon"]:
        sample_E_dict[codon] = 0
        sample_P_dict[codon] = 0
        sample_A_dict[codon] = 0
        populate_site_dict(sample_E_dict, "E", sample_df)
        populate_site_dict(sample_P_dict, "P", sample_df)
        populate_site_dict(sample_A_dict, "A", sample_df)

    print(sample_E_dict)
    print(sample_P_dict)
    print(sample_A_dict)

    return [sample_E_dict, sample_P_dict, sample_A_dict]

def compute_site_frequenices(sample_1, sample_2):
    # Index 0 = E, Index 1 = P, Index 2 = A
    codon_site_frequencies = [{}, {}, {}]
    # For sites E, P, and A (in that order)
    for site in range(0, len(sample_1)):
        # For every codon at current site between wild and mutant types
        for codon_key in sample_1[site].keys():
            # Save codon at current site as key, with (wild site occurrences)/(mutant site occurrences) as value
            codon_site_frequencies[site][codon_key] = sample_1[site][codon_key]/sample_2[site][codon_key]

    return codon_site_frequencies


def display_graph(df):
    # Index 0 = E, Index 1 = P, Index 2 = A
    sites = ["E", "P", "A"]

    for site in sites:
        temp_df = df[(df.site == site)]
        temp_df = temp_df.sort_values(by=['times_seen'], ascending=False)
        #temp_df = temp_df.sort_values(by=['codon'])
        #temp_df = temp_df.head(18)
        #print(temp_df)


        #sns.boxplot(x="codon", y="times_seen", hue="type", palette=["m", "g"], data=temp_df).set(title=site+' Site')


        sns.catplot(x="codon", y="times_seen", hue="type", palette=["m", "g"], data=temp_df, kind="box", height=8, aspect=2);
        matplotlib.pyplot.xticks(rotation=90)
        matplotlib.pyplot.show()




def frequencies_to_dataframe(codon_site_frequencies):
    # Index 0 = E, Index 1 = P, Index 2 = A
    sites = ["E", "P", "A"]
    frequency_df = pd.DataFrame(columns=['codon', 'site', 'frequency'])

    for site in range(0, len(sites)):
        for codon_key in codon_site_frequencies[site].keys():
            frequency_df.loc[len(frequency_df)] = [codon_key, sites[site], codon_site_frequencies[site][codon_key]]

    frequency_df = frequency_df.set_index('codon')
    return frequency_df


def main():
    # Creating dataframe to store all file-read results
    total_results = pd.DataFrame(columns=['read_name', 'codon', 'amino_acid', 'times_seen', 'site', 'contains_ambiguous', 'type'])
    total_results = total_results.set_index('read_name')

    # Reading all fastq files in current directory
    for entry in os.scandir():
        # If fastq file is found
        if entry.is_file() and entry.name.endswith('.fastq'):
            # Add latest file results to overall results
            total_results = pd.concat((total_results, read_file(entry.name)), axis=0)

    # Saving dataframe as a csv using the read name + ".csv"
    print(total_results)
    total_results.to_csv("total_results.csv")


# This is where the main code starts
if __name__ == "__main__":
    # main()
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    df = pd.read_csv("total_results.csv")
    display_graph(df)

    # Removing ambiguous nucleotide reads
    df = df[(df.contains_ambiguous == False)]

    #codon_freq = pd.read_csv("wild_vs_mutant_codon_site_frequencies.csv")

    #codon_freq = codon_freq.sort_values(by=['frequency'], ascending=False)
    #codon_freq = codon_freq.set_index("codon")
    #print(codon_freq)
    #quit(0)


    # Dataframe of only Wild Types
    wild_df = df[(df.type == "Wild")]
    mutant_df = df[(df.type == "Mutant")]

    # Acquiring codons per site in each sample type
    print("Computing Wild-Site-Occurrences")
    wild_site_list = get_site_occurrences(wild_df)
    print("Computing Mutant-Site-Occurrences")
    mutant_site_list = get_site_occurrences(mutant_df)
    print("Computing Codon-Site-Frequencies")
    codon_site_frequencies = compute_site_frequenices(wild_site_list, mutant_site_list)

    frequency_df = frequencies_to_dataframe(codon_site_frequencies)
    frequency_df.to_csv("wild_vs_mutant_codon_site_frequencies.csv")
