import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os

# Function used to identify and keep track of incoming codons
# Input: sequence, codon frequency dictionary, leftmost nucleotide in reading frame = (sequence length - 1) - offset
# NOTE: ASSUMES DNA DATA IS 5-TO-3 PRIME CODING STRANDS
def track_codons(inputseq, codonmap, codon_limit, offset):
    length = len(inputseq)  # Saving length of current sequence
    end = length - offset  # Starting index for codon frame read
    start = end - 2  # Ending index for codon frame read
    codonsread = 0  # Keeps track of codons read

    # While the frame has three nucleotides to read
    for codon_read in range(0, codon_limit):
        # Ensures reading frame doesn't run-off DNA fragment (ie accessing 3rd codon when only 6 nucleotides were given)
        if (start < 0) or (end-1 < 0) or (end < 0):
            break

        # Saving the current codon in the reading frame
        codon = inputseq[start] + inputseq[end-1] + inputseq[end]

        # Try incrementing number of times seen of current codon in frequency dictionary
        try:
            codonmap[codon] = codonmap[codon] + 1
        # If an exception is thrown, current codon was not seen yet, add it to dictionary
        except:
            codonmap[codon] = 1

        # Shift reading frame up three nucleotides
        start -= 3
        end -= 3

        # Update number of codons read
        codonsread += 1

    # Return number of codons read
    return codonsread


# Function used to convert codon frequency dictionary to a dataframe
# Input: Name of read, length of read, codon frequency dictionary
def create_df(readname, total_codons, codonmap):
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

    # Iterating through all codons saved in codon map (aka frequency dictionary)
    for codon in codonmap:
        localfrequency = 0  # Variable for storing frequency of current codon's appearance in the read it came from

        # Case for handling a 0 length read: if sequence length is not zero, calculate frequency
        if total_codons != 0:
            localfrequency = round(float(codonmap[codon])/total_codons, 5)

        # Creating a list to store current entry in dataframe
        amino_acid = ""
        try:
            amino_acid = amino_acid_table[codon]
        except:
            amino_acid = None

        # Cols: Name of read codon came from, codon sequence, times seen, frequency in read, contains an ambiguous read
        codondata = [readname, codon, amino_acid, codonmap[codon], localfrequency, "N" in codon]

        # Adding entry to dataframe list
        data.append(codondata.copy())

    # Convert dataframe list into a dataframe, removing unnecessary index column from dataframe that is set on default
    df = pd.DataFrame(data, columns=['read_name', 'codon', 'amino_acid', 'times_seen', 'local_codon_frequency', 'contains_ambiguous'])
    # Print and return dataframe
    df = df.set_index('read_name')
    return df


# Main function
def read_file(filename):
    # Dictionary used to keep track of seen codons (frequency dictionary)
    codonmap = {}
    codoncounter = 0

    # Iterates through all entries in a fastq file and tracks number of codons seen
    with open(filename) as in_handle:
        # For every sequence in fastq file
        for title, seq, qual in FastqGeneralIterator(in_handle):
            # Analyze all codons in current sequence, update number of codons seen
            codoncounter += track_codons(seq, codonmap, 3, 10)

        # Creating a name for this read by removing the .fastq portion of the input file
        readname = filename[0:len(filename) - 6:]
        # SQL has a problem with file names containing '-'. Ensuring these get replaced with '_' instead
        readname = readname.replace("-", "_")

        # Creating and returning dataframe from current file
        df = create_df(readname, codoncounter, codonmap)
        return df


def main():
    # Creating dataframe to store all file-read results
    total_results = pd.DataFrame(columns=['read_name', 'codon', 'amino_acid', 'times_seen', 'local_codon_frequency', 'contains_ambiguous'])
    total_results = total_results.set_index('read_name')

    # Reading all fastq files in current directory
    for entry in os.scandir():
        # If fastq file is found
        if entry.is_file() and entry.name.endswith('.fastq'):
            # Notify wich file is being read
            print("Reading", entry.name)
            # Add latest file results to overall results
            total_results = pd.concat((total_results, read_file(entry.name)), axis=0)

    # Saving dataframe as a csv using the read name + ".csv"
    print(total_results)
    total_results.to_csv("total_results.csv")


# This is where the main code starts
if __name__ == "__main__":
    main()


# TODO: Update global codon frequency?
