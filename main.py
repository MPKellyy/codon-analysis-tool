import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# Function used to identify and keep track of incoming codons
def TrackCodons(inputseq, codonmap):
    length = len(inputseq)
    start = 0
    end = 2
    codonsread = 0

    while end < length:
        codon = inputseq[start] + inputseq[start+1] + inputseq[end]

        try:
            codonmap[codon] = codonmap[codon] + 1
        except:
            codonmap[codon] = 1

        start += 3
        end += 3
        codonsread += 1

    return codonsread


def CreateDF(readname, readlength, codonmap):
    data = []
    for codon in codonmap:
        codondata = []
        codondata.append(readname)
        codondata.append(codon)
        codondata.append(codonmap[codon])
        localfrequency = 0;

        if readlength != 0:
            localfrequency = round(float(codonmap[codon])/readlength, 5)

        codondata.append(localfrequency)
        codondata.append("N" in codon)
        data.append(codondata.copy())

    df = pd.DataFrame(data, columns=['read_name', 'codon', 'times_seen', 'locol_codon_frequency', 'contains_ambiguous'])

    print(df)
    return df




# Dictionary used to keep track of seen codons
codonmap = {}
codoncounter = 0
filename = "RIBO-K90D1_R1_trimmed.fastq"

# Iterates through all entries in a fastq file and tracks number of codons seen
with open(filename) as in_handle:
     for title, seq, qual in FastqGeneralIterator(in_handle):
         codoncounter += TrackCodons(seq, codonmap)

     readname = filename[0:len(filename)-6:]
     readname = readname.replace("-", "_")
     df = CreateDF(readname, codoncounter, codonmap)
     df.to_csv(readname + ".csv")


# Might be a rule where certain number of nucleotides must be skipped before read
# TODO: Update codon frequency
# TODO: Remove ambiguous reads?







