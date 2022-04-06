import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
from tqdm import tqdm
import seaborn as sns


# Function used to identify and keep track of incoming codons
# Input: sequence, codon frequency dictionary, upstream position -1 of last nucleotide in site A
# NOTE: ASSUMES DNA DATA IS 5-TO-3 PRIME CODING STRANDS
def track_codons(inputseq, site_map, codon_limit, offset):
    length = len(inputseq)  # Saving length of current sequence

    # If offset is larger than sequence size, skip
    if (offset >= length-1) or (offset+1 >= length-1) or (offset+2 >= length-1):
        return

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


# Function used to convert codon frequency dictionary to a dataframe
# Input: Name of read, length of read, codon frequency dictionary
def create_df(readname, site_map):
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

    # For each site (Index 0: A, Index 1: P, Index 2: E)
    for site in range(0, len(site_map)):
        # Iterate through all codons seen at current site
        for codon in site_map[site]:
            # If the codon read is not ambiguous, map to its respective amino acid
            try:
                amino_acid = amino_acid_table[codon]
            # Else, ambiguous codon does not have an amino acid
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

# Function used to read an individual fastq file
# Input: String of file name
def read_file(filename):
    # Dictionary used to keep track of seen codons (frequency dictionary)
    site_map = [{}, {}, {}]  # Index 0: A, Index 1: P, Index 2: E

    # Iterates through all entries in a fastq file and tracks number of codons seen
    with open(filename) as in_handle:
        # For every sequence in fastq file
        for title, seq, qual in tqdm(FastqGeneralIterator(in_handle), desc=filename):
            # Analyze all codons in current sequence, update number of codons seen
            track_codons(seq, site_map, 3, 9)

        # Creating a name for this read by removing the .fastq portion of the input file
        readname = filename[0:len(filename) - 6:]
        # SQL has a problem with file names containing '-'. Ensuring these get replaced with '_' instead
        readname = readname.replace("-", "_")

        # Creating and returning dataframe from current file
        df = create_df(readname, site_map)
        return df


# Acquiring codons per site in sample type data
# Input: wild or sample type data frame
def get_site_occurrences(sample_df):
    # Index 0 = E, Index 1 = P, Index 2 = A
    sites = ["E", "P", "A"]  # List of sites
    codon_occurrences = [{}, {}, {}]  # List of site dictionaries to contain codon:occurrence pairs

    # For each site
    for site in range(0, len(sites)):
        # Filter data frame to current site only
        site_df = sample_df[(sample_df.site == sites[site])]

        # For each codon in filtered data
        for codon in site_df["codon"]:
            # Filter data frame to current codon
            codon_df = site_df[(sample_df.codon == codon)]
            occurrences = 0  # Variable to keep track of occurrences of current codon

            # For each read at current site and codon, sum the occurrences
            for occurrence in codon_df["times_seen"]:
                occurrences += occurrence

            # Save the codon-occurrence data at this site
            codon_occurrences[site][codon] = occurrences

    # Return codon site occurences
    return codon_occurrences


# Computing frequencies of codon sites between two sample types
def compute_site_frequenices(sample_1, sample_2):
    # Index 0 = E, Index 1 = P, Index 2 = A
    codon_site_frequencies = [{}, {}, {}]

    # For sites E, P, and A (in that order)
    for site in range(0, len(sample_1)):
        # For every codon at current site between wild and mutant types
        for codon_key in sample_1[site].keys():
            # Save codon at current site as key, with (wild site occurrences)/(mutant site occurrences) as value
            codon_site_frequencies[site][codon_key] = sample_1[site][codon_key]/sample_2[site][codon_key]

    # Return the codon-site frequencies
    return codon_site_frequencies


# Displays trend of codon site occurences between wild and mutant type data (excludes ambiguous reads)
# Input: Data frame containing results from extracted fastq files
def display_graph(df):
    # Index 0 = E, Index 1 = P, Index 2 = A
    sites = ["E", "P", "A"]  # List of sites

    # For each site
    for site in tqdm(sites, desc="Saving Graphs"):
        # Filter data frame to current site
        site_df = df[(df.site == site)]
        # Sort filtered data in descending order to visualize trend
        site_df = site_df.sort_values(by=['times_seen'], ascending=False)

        # Displaying results as a box plot
        sns.catplot(x="codon", y="times_seen", hue="type", palette=["m", "g"], data=site_df, kind="box", height=8, aspect=2);
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(site + "_Site_boxplot.png")
        plt.close()

        # Displaying results as a strip plot
        sns.stripplot(x="codon", y="times_seen", hue="type", palette=["m", "g"], data=site_df, jitter=False).set(title=site + ' Site')
        plt.xticks(rotation=90)
        figure = plt.gcf()
        figure.set_size_inches(16, 6)
        plt.savefig(site + "_Site_stripplot.png", dpi=250)
        plt.close()


# Converts codon-site frequency data to a data frame
# Input: List of dictionaries as codon-site data
def frequencies_to_dataframe(codon_site_frequencies):
    # Index 0 = E, Index 1 = P, Index 2 = A
    sites = ["E", "P", "A"]

    # Creating an empty frequency data frame
    frequency_df = pd.DataFrame(columns=['codon', 'site', 'frequency'])

    # For each site
    for site in range(0, len(sites)):
        # For each codon in site
        for codon_key in codon_site_frequencies[site].keys():
            # Transfer codon frequency to data frame
            frequency_df.loc[len(frequency_df)] = [codon_key, sites[site], codon_site_frequencies[site][codon_key]]

    # Aligning data frame index to codon column
    frequency_df = frequency_df.set_index('codon')

    # Return the frequency data frame
    return frequency_df


# Extracts all fastq files in directory for codon data
def extract_fastqs():
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
    total_results.to_csv("total_results.csv")

    # Returns data frame of all codon data from files
    return total_results


# Main logic
def main():
    # If directory fastqs were already read, load results (mainly for testing purposes)
    total_results = extract_fastqs()
    total_results = pd.read_csv("total_results.csv")

    # Removing ambiguous nucleotide reads
    total_results = total_results[(total_results.contains_ambiguous == False)]

    # Saves graphs of codon-site occurrence trends
    display_graph(total_results)

    # Dataframe of only Wild Types
    wild_df = total_results[(total_results.type == "Wild")]
    # Dataframe of only Mutant Types
    mutant_df = total_results[(total_results.type == "Mutant")]

    # Acquiring codons per site in each sample type
    wild_site_list = get_site_occurrences(wild_df)
    mutant_site_list = get_site_occurrences(mutant_df)

    # Compute and save codon-site frequencies across wild and mutant types
    codon_site_frequencies = compute_site_frequenices(wild_site_list, mutant_site_list)

    # Converting codon-site frequency data to a data frame
    frequency_df = frequencies_to_dataframe(codon_site_frequencies)
    # Sorting the data frame by frequency
    frequency_df = frequency_df.sort_values(by=['frequency'], ascending=False)
    # Save frequencies
    frequency_df.to_csv("wild_vs_mutant_codon_site_frequencies.csv")

    # Printing codon-site frequencies
    # pd.set_option('display.max_rows', None)
    # print(frequency_df)


# This is where the main code starts
if __name__ == "__main__":
    main()
