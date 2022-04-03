# codon-analysis-tool: code logic flow

Phase 1: Read Input Fastq Files
* main()
   * Main function has loop that identifies all fastq files in directory
* read_file()
   * Inside the main() file loop, each file gets passed into read_file()
   * read_file() opens the current file and parses through each sequence
   * read_file() uses a list of three dictionaries to keep track of codons at sites A, P, and E

Phase 2: Data Curation
* After each file is read in main(), all sequence data is appended into the total_data data frame
* Sequence data includes the type of codon seen, number of times codon appeared at each target site, wild-mutant type status, and ambiguous read flags
* After all files have been read, the total_data data frame has all ambiguous reads removed
* Two copies of this dataframe are made, and each restricts the results to wild and mutant type reads respectively

Phase 3: Frequency Comparisons
* Both the wild and mutant type data frames have their codon site data extracted via three site dictionaries that are unique to each sample category (3 dictionaries for wild type, 3 dictionaries for mutant type)
* After data is extracted, three sequential loops are conducted to acquire codon frequencies of A, P, and E sites between wild and mutant types (frequency = wild_site_num / mutant_site_number)

Phase 4: Frequency Plots
* All calculated frequencies are then plotted via seaborn
