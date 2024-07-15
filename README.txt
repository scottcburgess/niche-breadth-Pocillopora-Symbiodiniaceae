This README.txt file was generated on 2024-Jul-15 by Scott Burgess

GENERAL INFORMATION

1. Title of Dataset: Niche breadth and divergence in sympatric cryptic coral species (Pocillopora spp.) across habitats within reefs and among algal symbionts 

2. Author Information
	A. Principal Investigator Contact Information
		Name: Scott Burgess
		Institution: Florida State University
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306
		Email: sburgess@bio.fsu.edu


3. Date of data collection (single date, range, approximate date): 2021-2022 

4. Geographic location of data collection: Moorea, French Polynesia

5. Information about funding sources that supported the collection of the data: National Science Foundation (NSF; OCE-1829867).



DATA & FILE OVERVIEW

1. File List: 
ALGAE.R (This code produces Table 1, Table 3, Table 4, Figure 4, and Figure 6)
DISTRIBUTION.R (This code produces Table 3, Figure 2, Figure 3, Figure S1)
Figure 5.R (This code produces Figure 5)
NICHE.R (This code produces Figure 7, Figure 8)

257_20230414T204804_DBV_20230415T005314.profiles.absolute.abund_and_meta.csv
257_20230414T204804_DBV_20230415T005314.profiles.relative.abund_and_meta.csv
257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv
257_20230414T204804_DBV_20230415T005314.seqs.relative.abund_and_meta.csv
Distribution data.csv
mt_ORF_Organism file.csv
mt_ORF_source-table.csv
psbA_Organism file.csv
psbA_source-table.csv
mtORF Alignment 2021_22 samples.fasta
Cladocopium_psbA_2021_alignment_rooted samples only.fasta
257_20230414T204804_DBV_20230415T005314.seqs.fasta
20230414T020653_SymPortal.zip

2. Relationship between files: 
ALGAE.R uses	257_20230414T204804_DBV_20230415T005314.profiles.relative.abund_and_meta.csv,
			257_20230414T204804_DBV_20230415T005314.seqs.relative.abund_and_meta.csv
			Distribution data.csv

DISTRIBUTION.R uses Distribution data.csv

Figure 5.R uses 	257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv
				Distribution data.csv
				257_20230414T204804_DBV_20230415T005314.seqs.fasta

NICHE.R uses 	Distribution data.csv
			257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv


3. Metadata
257_20230414T204804_DBV_20230415T005314.profiles.absolute.abund_and_meta.csv
Output from SymPortal. Absolute abundance of ITS2 Type Profiles

257_20230414T204804_DBV_20230415T005314.profiles.relative.abund_and_meta.csv
Output from SymPortal. Relative abundance of ITS2 Type Profiles

257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv
Output from SymPortal. Absolute abundance of ITS2 DIV sequences

257_20230414T204804_DBV_20230415T005314.seqs.relative.abund_and_meta.csv
Output from SymPortal. Relative abundance of ITS2 DIV sequences

Distribution data.csv
For each sample (rows), the depth, site, host ID, and symbiont ID.
Trip: Unique identifier for each visit to Moorea
Date collected.YYYY-MM-DD: Date the sample was collected in YYYY-MM-DD format.
Preservative: DNA preservative
Depth.m: Depth in meters where the sample was collected
Site: Site ID where the sample was collected
Coral.ID:	Unique identifier for each sample
Box ID: Box ID where the tissue for each sample is stored at Florida State University
mtORF.RFLP: The mtORF haplotype, and results from RFLP in the case of haplotype 1, for each sample	Cladocopium_species: The species of Cladocopium symbiont hosted in each sample, identified using psbA_ncr
Cladocopium_clade: The clade of Cladocopium symbiont hosted in each sample, identified in Suppl. Fig. 2 and 3 of the paper. 
Plate Number: Plate ID where the DNA for each sample is stored at Florida State University

mt_ORF_Organism file.csv
Organism file for mtORF Alignment 2021_22 samples.fasta showing Pocillopora species identification for each mtORF Sequence ID (sample).

mt_ORF_source-table.csv
Additional data associated with each Sequence ID (sample) in mtORF Alignment 2021_22 samples.fasta

psbA_Organism file.csv
Organism file for Cladocopium_psbA_2021_alignment_rooted samples only.fasta showing Cladocopium species identification for each psbA_ncr Sequence ID (sample).

psbA_source-table.csv
Additional data associated with each Sequence ID (sample) in Cladocopium_psbA_2021_alignment_rooted samples only.fasta 

mtORF Alignment 2021_22 samples.fasta
FASTA file of aligned mtORF sequences

Cladocopium_psbA_2021_alignment_rooted samples only.fasta
FASTA file of aligned psbA_ncr sequences

257_20230414T204804_DBV_20230415T005314.seqs.fasta
FASTA file of ITS2 sequences

20230414T020653_SymPortal.zip
Compressed file of complete SymPortal output


