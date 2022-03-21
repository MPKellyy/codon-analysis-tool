import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# Function used to identify and keep track of incoming codons
# Input: sequence, codon frequency dictionary, starting index
def track_codons(inputseq, codonmap, offset):
    length = len(inputseq)  # Saving length of current sequence
    start = offset  # Starting index for codon frame read
    end = start + 2  # Ending index for codon frame read
    codonsread = 0  # Keeps track of codons read

    # While the frame has three nucleotides to read
    while end < length:
        # Saving the current codon in the reading frame
        codon = inputseq[start] + inputseq[start+1] + inputseq[end]

        # Try incrementing number of times seen of current codon in frequency dictionary
        try:
            codonmap[codon] = codonmap[codon] + 1
        # If an exception is thrown, current codon was not seen yet, add it to dictionary
        except:
            codonmap[codon] = 1

        # Shift reading frame up three nucleotides
        start += 3
        end += 3

        # Update number of codons read
        codonsread += 1

    # Return number of codons read
    return codonsread


# Function used to convert codon frequency dictionary to a dataframe
# Input: Name of read, length of read, codon frequency dictionary
def create_df(readname, readlength, codonmap):
    # Creating an empty list to hold all entries in the dataframe
    data = []

    # Iterating through all codons saved in codon map (aka frequency dictionary)
    for codon in codonmap:
        localfrequency = 0  # Variable for storing frequency of current codon's appearance in the read it came from

        # Case for handling a 0 length read: if sequence length is not zero, calculate frequency
        if readlength != 0:
            localfrequency = round(float(codonmap[codon])/readlength, 5)

        # Creating a list to store current entry in dataframe
        # Cols: Name of read codon came from, codon sequence, times seen, frequency in read, contains an ambiguous read
        codondata = [readname, codon, codonmap[codon], localfrequency, "N" in codon]

        # Adding entry to dataframe list
        data.append(codondata.copy())

    # Convert dataframe list into a dataframe, removing unnecessary index column from dataframe that is set on default
    df = pd.DataFrame(data, columns=['read_name', 'codon', 'times_seen', 'locol_codon_frequency', 'contains_ambiguous'])
    # Print and return dataframe
    df = df.set_index('read_name')
    print(df)
    return df


# Main function
def main():
    # Dictionary used to keep track of seen codons (frequency dictionary)
    codonmap = {}
    codoncounter = 0

    # Hardcoded filename, will replace in future updates
    filename = "RIBO-K90D1_R1_trimmed.fastq"

    # Iterates through all entries in a fastq file and tracks number of codons seen
    with open(filename) as in_handle:
        # For every sequence in fastq file
        for title, seq, qual in FastqGeneralIterator(in_handle):
            # Analyze all codons in current sequence, update number of codons seen
            codoncounter += track_codons(seq, codonmap, 0)

        # Creating a name for this read by removing the .fastq portion of the input file
        readname = filename[0:len(filename) - 6:]
        # SQL has a problem with file names containing '-'. Ensuring these get replaced with '_' instead
        readname = readname.replace("-", "_")

        # Creating a dataframe of current sequence
        df = create_df(readname, codoncounter, codonmap)
        # Saving dataframe as a csv using the read name + ".csv"
        df.to_csv(readname + ".csv")


# This is where the main code starts
if __name__ == "__main__":
    main()


# Might be a rule where certain number of nucleotides must be skipped before read
# TODO: Update global codon frequency?
# TODO: Remove ambiguous reads?
# TODO: Replace hardcoded filename with an automated file read algorithm
