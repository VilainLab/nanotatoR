#' Database of Genomic Variants
#'
#' A database containing Copy Number variant information
#' for humans. This is a subset of the database containing the 
#' first 300 lines.
#'
#' @format A data frame with 300 rows and 16 variables:
#' \describe{
#'   \item{variantaccession}{Accession number for the variant}
#'   \item{chr}{Chromosome Number}
#'   \item{start}{Structural Variant Breakpoint Start}
#'   \item{end}{Structural Variant Breakpoint End}
#'   \item{varianttype}{Type of Variant}
#'   \item{variantsubtype}{Subtype of variant}
#'   \item{reference}{Reference of the study which discovered the variant}
#'   \item{pubmedid}{Pubmed id of the study which discovered the variant}
#'   \item{method}{method used to discover variants}
#'   \item{supportingvariants}{supporting variants for the variants}
#'   \item{mergedorsample}{Are variants discovered from one sample?}
#'   \item{samplesize}{Size of the sample}
#'   \item{observedgains}{Number of observed gains/insertions}
#'   \item{observedlosses}{Number of observed losses/deletion}
#'   \item{genes}{Genes overlapping variants}
#'   \item{samples}{Name of samples where this is observed}
#' }
#' @source \url{http://dgv.tcag.ca/dgv/app/downloads?ref=GRCh37/hg19/GRCh37_hg19_variants_2016-05-15.txt}
'GRCh37_hg19_variants_2016-05-15.txt'
#' Decipher Database
#'
#' A database containing Copy Number variant information
#' for humans. This is a subset of the database containing the 
#' first 300 lines.
#'
#' @format A data frame with 300 rows and 16 variables:
#' \describe{
#'   \item{population_cnv_id}{Unique CNV ID for the variants}
#'   \item{chr}{Chromosome Number}
#'   \item{start}{Structural Variant Breakpoint Start}
#'   \item{end}{Structural Variant Breakpoint End}
#'   \item{deletion_observations}{Number of deletions observed}
#'   \item{deletion_frequency}{Frequency of the deletions}
#'   \item{deletion_standard_error}{Standard error in deletion frequency calculations}
#'   \item{duplication_observations}{Number of duplications observed}
#'   \item{duplication_frequency}{Frequency of the duplications}
#'   \item{duplication_standard_error}{Standard error in duplication frequency calculations}
#'   \item{observations}{Total number of observations (deletions + duplications) for the cnv id}
#'   \item{frequency}{Frequency of occurence for the cnv id}
#'   \item{standard_error}{Standard error of the same}
#'   \item{type}{Type of the variants}
#'   \item{sample_size}{Size of the samples}
#'   \item{study}{Name of studies where this is observed}
#' }
#' @source \url{https://decipher.sanger.ac.uk/about#downloads/data/population_cnv.txt}
"population_cnv.txt"
#' F1.1_TestSample1_solo_hg19
#'
#' A sample structural variant output file from bionano
#' for humans. This is a subset of the genome in a bottle 
#' entry for Ashkenazim family proband(HG002_NA24385_son).
#' Assembly (all.bnx) files were downloaded from the genome 
#' in a bottle ftp site 
#' (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/BioNano/all.bnx)
#' and Bionano SV caller Solve was used to call SVs. File was renamed to F1.1_TestSample1_solo_hg19.smap
#' to match naming convention for nanotatoR. This is the input file
#'
#' @format A data frame with 148 rows and 39 variables:
#' \describe{
#'   \item{SmapEntryID}{Unique smap entry ID for the variants}
#'   \item{QryContigID}{Unique contig id of the assembly}
#'   \item{RefcontigID1}{Chromosome Number corresponding to reference}
#'   \item{RefcontigID2}{Chromosome Number corresponding to reference. Same as RefcontigID1 for most cases 
#'   excepting translocation}
#'   \item{QryStartPos}{Start position of the Structural variant on the contig}
#'   \item{QryEndPos}{End position of the Structural variant on the contig}
#'   \item{RefStartPos}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{RefEndPos}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Confidence}{Confidence score generated for variants}
#'   \item{Type}{Type of structural variants}
#'   \item{Type1}{Type of variants observed from enzyme 1}
#'   \item{Type2}{Type of variants observed from enzyme 2}
#'   \item{XmapID1}{Bionano alignment ID for enzyme 1}
#'   \item{XmapID2}{Bionano alignment ID for enzyme 2}
#'   \item{LinkID}{If two Smap entries are linked, this ID links them}
#'   \item{QryStartIdx}{ID copied from smap entry for Query start during merging multiple smaps for SV calling}
#'   \item{QryEndIdx}{ID copied from smap entry for Query end during merging multiple smaps for SV calling}
#'   \item{RefStartIdx}{ID copied from smap entry for Reference start during merging multiple smaps for SV calling}
#'   \item{RefEndIdx}{ID copied from smap entry for Reference end during merging multiple smaps for SV calling}
#'   \item{Zygosity}{Zygosity of the Structural variants}
#'   \item{Genotype}{Genotype of the sample. '1' for homozygous SVs, '1' or '2' for heterozygous SVs, and '-1' for #'   unknown zygosity and SVs which are not indels or translocations}
#'   \item{GenotypeGroup}{Indels which overlap one another are assigned the same group}
#'   \item{RawConfidence}{Minimum of next three columns for indels. '-1' for other SV types}
#'   \item{RawConfidenceLeft}{Confidence of alignment to the left of indel or translocation}
#'   \item{RawConfidenceRight}{Confidence of alignment to the right of indel or translocation}
#'   \item{RawConfidenceCentre}{Confidence of alignment to the centre of indel or translocation}
#'   \item{Sample}{Name of Sample}
#'   \item{Algorithm}{Algorithm used.}
#'   \item{Size}{Size of the variant.}
#'   \item{Present_in_%_of_BNG_control_samples}{percentage of overlap between variants in the sample and in the 
#'   database}
#'   \item{Present_in_%_of_BNG_control_samples_with_the_same_enzyme}{percentage of overlap between variants in the #'   sample and in the database, for the same enzyme}
#'   \item{Fail_BSPQI_assembly_chimeric_score}{enzyme1 BSPQI chimeric score pass or fail}
#'   \item{Fail_BSSSI_assembly_chimeric_score}{enzyme2 BSSSI chimeric score pass or fail}
#'   \item{OverlapGenes}{Genes overlapping the SV region}
#'   \item{NearestNonOverlapGene}{Nearest genes to the the SV region}
#'   \item{NearestNonOverlapGeneDistance}{Nearest non overlaping gene distance from the SV breakpoints}
#'   \item{Found_in_parents_assemblies}{Whether they are found in parents assemblies. It can be "mother", 
#'   "father", "both" or "none"}
#'   \item{Found_in_parents_molecules}{Whether they are found in parents molecules. It can be "mother", "father", #'   "both" or "none"}
#'   \item{Found_in_self_BSPQI_molecules}{Whether they are found in parents enzyme 1 (BSPQI) molecules. It can be #'   "mother", "father", "both" or "none"}
#'   \item{Found_in_self_BSSSI_molecules}{Whether they are found in parents enzyme 1 (BSSSI) molecules. It can be #'   "mother", "father", "both" or "none"}
#' }
"F1.1_TestSample1_solo_hg19.smap"
#' F1.2_TestSample2_solo_hg19
#'
#' A sample structural variant output file from bionano
#' for humans. This is a subset of the genome in a bottle 
#' entry for Ashkenazim family proband(HG004_NA24143_mother).
#' Assembly (all.bnx) files were downloaded from the genome 
#' in a bottle ftp site 
#' (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/BioNano/all.bnx)
#' and Bionano SV caller Solve was used to call SVs. File was renamed to F1.1_TestSample1_solo_hg19.smap
#' to match naming convention for nanotatoR. This is the file in Solo file folder for internal frequency
#' calculations.
#'
#' @format A data frame with 1000 rows and 40 variables:
#' \describe{
#'   \item{SmapEntryID}{Unique smap entry ID for the variants}
#'   \item{QryContigID}{Unique contig id of the assembly}
#'   \item{RefcontigID1}{Chromosome Number corresponding to reference}
#'   \item{RefcontigID2}{Chromosome Number corresponding to reference. Same as RefcontigID1 for most cases 
#'   excepting translocation}
#'   \item{QryStartPos}{Start position of the Structural variant on the contig}
#'   \item{QryEndPos}{End position of the Structural variant on the contig}
#'   \item{RefStartPos}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{RefEndPos}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Confidence}{Confidence score generated for variants}
#'   \item{Type}{Type of structural variants}
#'   \item{Type1}{Type of variants observed from enzyme 1}
#'   \item{Type2}{Type of variants observed from enzyme 2}
#'   \item{XmapID1}{Bionano alignment ID for enzyme 1}
#'   \item{XmapID2}{Bionano alignment ID for enzyme 2}
#'   \item{LinkID}{If two Smap entries are linked, this ID links them}
#'   \item{QryStartIdx}{ID copied from smap entry for Query start during merging multiple smaps for SV calling}
#'   \item{QryEndIdx}{ID copied from smap entry for Query end during merging multiple smaps for SV calling}
#'   \item{RefStartIdx}{ID copied from smap entry for Reference start during merging multiple smaps for SV calling}
#'   \item{RefEndIdx}{ID copied from smap entry for Reference end during merging multiple smaps for SV calling}
#'   \item{Zygosity}{Zygosity of the Structural variants}
#'   \item{Genotype}{Genotype of the sample. '1' for homozygous SVs, '1' or '2' for heterozygous SVs, and '-1' for #'   unknown zygosity and SVs which are not indels or translocations}
#'   \item{GenotypeGroup}{Indels which overlap one another are assigned the same group}
#'   \item{RawConfidence}{Minimum of next three columns for indels. '-1' for other SV types}
#'   \item{RawConfidenceLeft}{Confidence of alignment to the left of indel or translocation}
#'   \item{RawConfidenceRight}{Confidence of alignment to the right of indel or translocation}
#'   \item{RawConfidenceCentre}{Confidence of alignment to the centre of indel or translocation}
#'   \item{Sample}{Name of Sample}
#'   \item{Size}{Size of the variant.}
#'   \item{Algorithm}{Algorithm used.}
#'   \item{Present_in_%_of_BNG_control_samples}{percentage of overlap between variants in the sample and in the 
#'   database}
#'   \item{Present_in_%_of_BNG_control_samples_with_the_same_enzyme}{percentage of overlap between variants in the #'   sample and in the database, for the same enzyme}
#'   \item{Fail_BSPQI_assembly_chimeric_score}{enzyme1 BSPQI chimeric score pass or fail}
#'   \item{Fail_BSSSI_assembly_chimeric_score}{enzyme2 BSSSI chimeric score pass or fail}
#'   \item{OverlapGenes}{Genes overlapping the SV region}
#'   \item{NearestNonOverlapGene}{Nearest genes to the the SV region}
#'   \item{NearestNonOverlapGeneDistance}{Nearest non overlaping gene distance from the SV breakpoints}
#'   \item{Found_in_parents_assemblies}{Whether they are found in parents assemblies. It can be "mother", 
#'   "father",    "both" or "none"}
#'   \item{Found_in_parents_molecules}{Whether they are found in parents molecules. It can be "mother", "father", #'   "both" or "none"}
#'   \item{Found_in_self_BSPQI_molecules}{Whether they are found in parents enzyme 1 (BSPQI) molecules. It can be #'   "mother", "father", "both" or "none"}
#'   \item{Found_in_self_BSSSI_molecules}{Whether they are found in parents enzyme 1 (BSSSI) molecules. It can be #'   "mother", "father", "both" or "none"}
#' }
"F1.2_TestSample1_solo_hg19.smap"
#' F1.3_TestSample2_solo_hg19
#'
#' A sample structural variant output file from bionano
#' for humans. This is a subset of the genome in a bottle 
#' entry for Ashkenazim family proband(HG004_NA24143_mother).
#' Assembly (all.bnx) files were downloaded from the genome 
#' in a bottle ftp site 
#' (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/BioNano/all.bnx)
#' and Bionano SV caller Solve was used to call SVs. File was renamed to F1.3_TestSample1_solo_hg19.smap
#' to match naming convention for nanotatoR. This is the file in Solo file folder for internal frequency
#' calculations.
#'
#' @format A data frame with 1000 rows and 40 variables:
#' \describe{
#'   \item{SmapEntryID}{Unique smap entry ID for the variants}
#'   \item{QryContigID}{Unique contig id of the assembly}
#'   \item{RefcontigID1}{Chromosome Number corresponding to reference}
#'   \item{RefcontigID2}{Chromosome Number corresponding to reference. Same as RefcontigID1 for most cases 
#'   excepting translocation}
#'   \item{QryStartPos}{Start position of the Structural variant on the contig}
#'   \item{QryEndPos}{End position of the Structural variant on the contig}
#'   \item{RefStartPos}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{RefEndPos}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Confidence}{Confidence score generated for variants}
#'   \item{Type}{Type of structural variants}
#'   \item{Type1}{Type of variants observed from enzyme 1}
#'   \item{Type2}{Type of variants observed from enzyme 2}
#'   \item{XmapID1}{Bionano alignment ID for enzyme 1}
#'   \item{XmapID2}{Bionano alignment ID for enzyme 2}
#'   \item{LinkID}{If two Smap entries are linked, this ID links them}
#'   \item{QryStartIdx}{ID copied from smap entry for Query start during merging multiple smaps for SV calling}
#'   \item{QryEndIdx}{ID copied from smap entry for Query end during merging multiple smaps for SV calling}
#'   \item{RefStartIdx}{ID copied from smap entry for Reference start during merging multiple smaps for SV calling}
#'   \item{RefEndIdx}{ID copied from smap entry for Reference end during merging multiple smaps for SV calling}
#'   \item{Zygosity}{Zygosity of the Structural variants}
#'   \item{Genotype}{Genotype of the sample. '1' for homozygous SVs, '1' or '2' for heterozygous SVs, and '-1' for #'   unknown zygosity and SVs which are not indels or translocations}
#'   \item{GenotypeGroup}{Indels which overlap one another are assigned the same group}
#'   \item{RawConfidence}{Minimum of next three columns for indels. '-1' for other SV types}
#'   \item{RawConfidenceLeft}{Confidence of alignment to the left of indel or translocation}
#'   \item{RawConfidenceRight}{Confidence of alignment to the right of indel or translocation}
#'   \item{RawConfidenceCentre}{Confidence of alignment to the centre of indel or translocation}
#'   \item{Sample}{Name of Sample}
#'   \item{Algorithm}{Algorithm used.}
#'   \item{Size}{Size of the variant.}
#'   \item{Present_in_%_of_BNG_control_samples}{percentage of overlap between variants in the sample and in the 
#'   database}
#'   \item{Present_in_%_of_BNG_control_samples_with_the_same_enzyme}{percentage of overlap between variants in the #'   sample and in the database, for the same enzyme}
#'   \item{Fail_BSPQI_assembly_chimeric_score}{enzyme1 BSPQI chimeric score pass or fail}
#'   \item{Fail_BSSSI_assembly_chimeric_score}{enzyme2 BSSSI chimeric score pass or fail}
#'   \item{OverlapGenes}{Genes overlapping the SV region}
#'   \item{NearestNonOverlapGene}{Nearest genes to the the SV region}
#'   \item{NearestNonOverlapGeneDistance}{Nearest non overlaping gene distance from the SV breakpoints}
#'   \item{Found_in_parents_assemblies}{Whether they are found in parents assemblies. It can be "mother", 
#'   "father", "both" or "none"}
#'   \item{Found_in_parents_molecules}{Whether they are found in parents molecules. It can be "mother", "father", #'   "both" or "none"}
#'   \item{Found_in_self_BSPQI_molecules}{Whether they are found in parents enzyme 1 (BSPQI) molecules. It can be #'   "mother", "father", "both" or "none"}
#'   \item{Found_in_self_BSSSI_molecules}{Whether they are found in parents enzyme 1 (BSSSI) molecules. It can be #'   "mother", "father", "both" or "none"}
#' }
"F1.3_TestSample2_solo_hg19.smap"
#' F2.2_TestSample2_solo_hg19
#'
#' A sample structural variant output file from bionano
#' for humans. This is a subset of the genome in a bottle 
#' entry for Ashkenazim family proband(HG004_NA24143_mother).
#' Assembly (all.bnx) files were downloaded from the genome 
#' in a bottle ftp site 
#' (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/MtSinai_BioNano/all.bnx)
#' and Bionano SV caller Solve was used to call SVs. File was renamed to F2.2_TestSample1_solo_hg19.smap
#' to match naming convention for nanotatoR. This is the file in Solo file folder for internal frequency
#' calculations.
#'
#' @format A data frame with 1000 rows and 40 variables:
#' \describe{
#'   \item{SmapEntryID}{Unique smap entry ID for the variants}
#'   \item{QryContigID}{Unique contig id of the assembly}
#'   \item{RefcontigID1}{Chromosome Number corresponding to reference}
#'   \item{RefcontigID2}{Chromosome Number corresponding to reference. Same as RefcontigID1 for most cases 
#'   excepting translocation}
#'   \item{QryStartPos}{Start position of the Structural variant on the contig}
#'   \item{QryEndPos}{End position of the Structural variant on the contig}
#'   \item{RefStartPos}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{RefEndPos}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Confidence}{Confidence score generated for variants}
#'   \item{Type}{Type of structural variants}
#'   \item{Type1}{Type of variants observed from enzyme 1}
#'   \item{Type2}{Type of variants observed from enzyme 2}
#'   \item{XmapID1}{Bionano alignment ID for enzyme 1}
#'   \item{XmapID2}{Bionano alignment ID for enzyme 2}
#'   \item{LinkID}{If two Smap entries are linked, this ID links them}
#'   \item{QryStartIdx}{ID copied from smap entry for Query start during merging multiple smaps for SV calling}
#'   \item{QryEndIdx}{ID copied from smap entry for Query end during merging multiple smaps for SV calling}
#'   \item{RefStartIdx}{ID copied from smap entry for Reference start during merging multiple smaps for SV calling}
#'   \item{RefEndIdx}{ID copied from smap entry for Reference end during merging multiple smaps for SV calling}
#'   \item{Zygosity}{Zygosity of the Structural variants}
#'   \item{Genotype}{Genotype of the sample. '1' for homozygous SVs, '1' or '2' for heterozygous SVs, and '-1' for #'   unknown zygosity and SVs  which are not indels or translocations}
#'   \item{GenotypeGroup}{Indels which overlap one another are assigned the same group}
#'   \item{RawConfidence}{Minimum of next three columns for indels. '-1' for other SV types}
#'   \item{RawConfidenceLeft}{Confidence of alignment to the left of indel or translocation}
#'   \item{RawConfidenceRight}{Confidence of alignment to the right of indel or translocation}
#'   \item{RawConfidenceCentre}{Confidence of alignment to the centre of indel or translocation}
#'   \item{Sample}{Name of Sample}
#'   \item{Algorithm}{Algorithm used.}
#'   \item{Size}{Size of the variant.}
#'   \item{Present_in_%_of_BNG_control_samples}{percentage of overlap between variants in the sample and in the 
#'   database}
#'   \item{Present_in_%_of_BNG_control_samples_with_the_same_enzyme}{percentage of overlap between variants in the #'   sample and in the database, for the same enzyme}
#'   \item{Fail_BSPQI_assembly_chimeric_score}{enzyme1 BSPQI chimeric score pass or fail}
#'   \item{Fail_BSSSI_assembly_chimeric_score}{enzyme2 BSSSI chimeric score pass or fail}
#'   \item{OverlapGenes}{Genes overlapping the SV region}
#'   \item{NearestNonOverlapGene}{Nearest genes to the the SV region}
#'   \item{NearestNonOverlapGeneDistance}{Nearest non overlaping gene distance from the SV breakpoints}
#'   \item{Found_in_parents_assemblies}{Whether they are found in parents assemblies. It can be "mother",
#'   "father", "both" or "none"}
#'   \item{Found_in_parents_molecules}{Whether they are found in parents molecules. It can be "mother", "father", #'   "both" or "none"}
#'   \item{Found_in_self_BSPQI_molecules}{Whether they are found in parents enzyme 1 (BSPQI) molecules. It can be #'   "mother", "father", "both" or "none"}
#'   \item{Found_in_self_BSSSI_molecules}{Whether they are found in parents enzyme 1 (BSSSI) molecules. It can be #'   "mother", "father", "both" or "none"}
#' }
"F2.2_TestSample2_solo_hg19.smap"
#' nanotatoR_merged
#'
#' Merged variant map file created by combining all the smap files in the SoloFile folder
#' using the makeMergedSVData function from nanotatoR.
#'
#' @format A data frame with 3965 rows and 40 variables:
#' \describe{
#'   \item{SmapEntryID}{Unique smap entry ID for the variants}
#'   \item{QryContigID}{Unique contig id of the assembly}
#'   \item{RefcontigID1}{Chromosome Number corresponding to reference}
#'   \item{RefcontigID2}{Chromosome Number corresponding to reference. Same as RefcontigID1 for most cases 
#'   excepting translocation}
#'   \item{QryStartPos}{Start position of the Structural variant on the contig}
#'   \item{QryEndPos}{End position of the Structural variant on the contig}
#'   \item{RefStartPos}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{RefEndPos}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Confidence}{Confidence score generated for variants}
#'   \item{Type}{Type of structural variants}
#'   \item{Type1}{Type of variants observed from enzyme 1}
#'   \item{Type2}{Type of variants observed from enzyme 2}
#'   \item{XmapID1}{Bionano alignment ID for enzyme 1}
#'   \item{XmapID2}{Bionano alignment ID for enzyme 2}
#'   \item{LinkID}{If two Smap entries are linked, this ID links them}
#'   \item{QryStartIdx}{ID copied from smap entry for Query start during merging multiple smaps for SV calling}
#'   \item{QryEndIdx}{ID copied from smap entry for Query end during merging multiple smaps for SV calling}
#'   \item{RefStartIdx}{ID copied from smap entry for Reference start during merging multiple smaps for SV calling}
#'   \item{RefEndIdx}{ID copied from smap entry for Reference end during merging multiple smaps for SV calling}
#'   \item{Zygosity}{Zygosity of the Structural variants}
#'   \item{Genotype}{Genotype of the sample. '1' for homozygous SVs, '1' or '2' for heterozygous SVs, and '-1' for #'   unknown zygosity and SVs which are not indels or translocations}
#'   \item{GenotypeGroup}{Indels which overlap one another are assigned the same group}
#'   \item{RawConfidence}{Minimum of next three columns for indels. '-1' for other SV types}
#'   \item{RawConfidenceLeft}{Confidence of alignment to the left of indel or translocation}
#'   \item{RawConfidenceRight}{Confidence of alignment to the right of indel or translocation}
#'   \item{RawConfidenceCentre}{Confidence of alignment to the centre of indel or translocation}
#'   \item{Sample}{Name of Sample}
#'   \item{Algorithm}{Algorithm used.}
#'   \item{Size}{Size of the variant.}
#'   \item{Present_in_%_of_BNG_control_samples}{percentage of overlap between variants in the sample and in the 
#'   database}
#'   \item{Present_in_%_of_BNG_control_samples_with_the_same_enzyme}{percentage of overlap between variants in the #'   sample and in the database, for the same enzyme}
#'   \item{Fail_BSPQI_assembly_chimeric_score}{enzyme1 BSPQI chimeric score pass or fail}
#'   \item{Fail_BSSSI_assembly_chimeric_score}{enzyme2 BSSSI chimeric score pass or fail}
#'   \item{OverlapGenes}{Genes overlapping the SV region}
#'   \item{NearestNonOverlapGene}{Nearest genes to the the SV region}
#'   \item{NearestNonOverlapGeneDistance}{Nearest non overlaping gene distance from the SV breakpoints}
#'   \item{Found_in_parents_assemblies}{Whether they are found in parents assemblies. It can be "mother", 
#'   "father", "both" or "none"}
#'   \item{Found_in_parents_molecules}{Whether they are found in parents molecules. It can be "mother", "father", #'   "both" or "none"}
#'   \item{Found_in_self_BSPQI_molecules}{Whether they are found in parents enzyme 1 (BSPQI) molecules. It can be #'   "mother", "father", "both" or "none"}
#'   \item{Found_in_self_BSSSI_molecules}{Whether they are found in parents enzyme 1 (BSSSI) molecules. It can be #'   "mother", "father", "both" or "none"}
#' }
"nanotatoR_merged.txt"
#' current_ctrl_dup_hg19_anonymize
#'
#' A subset of the Bionano internal control data set of 204 patients, for duplication.
#' The Bionano genomic reference database can be downloaded using the following 
#' command wget http://bnxinstall.com/solve/Solve3.3_10252018.tar.gz , followed by 
#' decompressing it using the tar -xvzf Solve3.3_10252018.tar.gz. 
#' The folder containing the database is in the config file 
#' which can be found in the folowing directory 
#' $PWD/Solve3.3_10252018/VariantAnnotation/10252018/config/. 
#' The reference files are named  based on the variant type and reference genome. 
#' For example: current_ctrl_dup_hg19_anonymize.txt, 
#' would be the duplication reference variant file for hg19 reference genome.
#'
#' @format A data frame with 199 rows and 10 variables:
#' \describe{
#'   \item{sample}{Name of Sample}
#'   \item{type}{Type of structural variants}
#'   \item{algorithm}{Algorithm used.}
#'   \item{chr}{Chromosome Nummber.}
#'   \item{start}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{end}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{zygosity}{Zygosity of the variant.}
#'   \item{size}{Size of the variant.}
#'   \item{score}{Confidence score generated for variants}
#'   \item{smapId}{Unique smap ID}
#' }
"current_ctrl_dup_hg19_anonymize.txt"
#' current_ctrl_ins_del_hg19_anonymize
#'
#' A subset of the Bionano internal control data set of 204 patients, for insertions
#' and duplications.
#' This was downloaded from bionano support database.
#' See Vignettes or duplication data header for download instruction.
#'
#' @format A data frame with 200 rows and 10 variables:
#' \describe{
#'   \item{sample}{Name of Sample}
#'   \item{type}{Type of structural variants}
#'   \item{algorithm}{Algorithm used.}
#'   \item{chr}{Chromosome Nummber.}
#'   \item{start}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{end}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{zygosity}{Zygosity of the variant.}
#'   \item{size}{Size of the variant.}
#'   \item{score}{Confidence score generated for variants}
#'   \item{smapId}{Unique smap ID}
#' }
"current_ctrl_ins_del_hg19_anonymize.txt"
#' current_ctrl_inv_hg19_anonymize
#'
#' This was downloaded from bionano support database.
#' See Vignettes or duplication data header for download instruction.
#'
#' @format A data frame with 200 rows and 11 variables:
#' \describe{
#'   \item{sample}{Name of Sample}
#'   \item{type}{Type of structural variants}
#'   \item{algorithm}{Algorithm used.}
#'   \item{chr}{Chromosome Nummber.}
#'   \item{start}{Start position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{end}{End position of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{size}{Size of the variant.}
#'   \item{zygosity}{Zygosity of the variant.}
#'   \item{score}{Confidence score generated for variants}
#'   \item{smapId}{Unique smap ID}
#'   \item{linkSmapId}{smap ID of linking variants}
#' }
"current_ctrl_inv_hg19_anonymize.txt"
#' current_ctrl_trans_hg19_anonymize
#'
#' A subset of the Bionano internal control data set of 204 patients, for translocations
#' This was downloaded from bionano support database.
#' See Vignettes or duplication data header for download instruction.
#'
#' @format A data frame with 200 rows and 10 variables:
#' \describe{
#'   \item{sample}{Name of Sample}
#'   \item{type}{Type of structural variants}
#'   \item{algorithm}{Algorithm used.}
#'   \item{chrA}{Chromosome Nummber for the first breakpoint.}
#'   \item{bkptA}{start location of the Structural variant on the human Reference genome (eg. hg19/hg38) for 1st #'   point}
#'   \item{chrB}{Chromosome Nummber for the second breakpoint.}
#'   \item{bkptB}{End location of the Structural variant on the human Reference genome (eg. hg19/hg38) for 2nd 
#'   point}
#'   \item{zygosity}{Zygosity of the variant.}
#'   \item{score}{Confidence score generated for variants}
#'   \item{smapId}{Unique smap ID}
#' }
"current_ctrl_trans_hg19_anonymize.txt"
#' BNSOLO2_merged
#'
#' Combining individual bionano internals datasets downloaded from the bionano support 
#' site into one combined dataset using makeMergedSmapData function.
#'
#' @format A data frame with 250 rows and 10 variables:
#' \describe{
#'   \item{SmapId}{Unique smap ID}
#'   \item{Sample}{Name of Sample}
#'   \item{Type}{Type of structural variants}
#'   \item{Chromosome1}{Chromosome Nummber for the first breakpoint.}
#'   \item{Chromosome2}{Chromosome Nummber for the second breakpoint.}
#'   \item{BreakPntStrt}{start location of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{BreakPntEnd}{end location of the Structural variant on the human Reference genome (eg. hg19/hg38)}
#'   \item{Size}{Size of the variant.}
#'   \item{Zygosity}{Zygosity of the variant.}
#'   \item{Score}{Confidence score generated for variants}
#' }
"BNSOLO2_merged.txt"
#' Homo_sapiens.Hg19
#'
#' Bed file generated from gencode gtf file using the linux command line  
#' awk '{print $1,$4,$5,$18,$7}' gencode.v19.annotation.gtf>Homo_sapiens.Hg19.bed
#'
#' @format A data frame with 331 rows and 5 variables:
#' \describe{
#'   \item{chromosome}{chromosome number}
#'   \item{ChromsomeStart}{Gene Start}
#'   \item{ChromsomeEnd}{Gene Start}
#'   \item{Gene}{Gene Symbol}
#'   \item{strand}{Strand information}
#' }
"Homo_sapiens.Hg19.bed"
#' Homo_sapiens.Hg19_BN_bed.txt
#'
#' Bed file generated from by buildrunBNBedFiles from nanotatoR  
#'
#' @format A data frame with 331 rows and 5 variables:
#' \describe{
#'   \item{chromosome}{chromosome number}
#'   \item{ChromsomeStart}{Gene Start}
#'   \item{ChromsomeEnd}{Gene Start}
#'   \item{Gene}{Gene Symbol}
#'   \item{strand}{Strand information}
#'   \item{ChromsomeStart}{Gene Start}
#'   \item{ChromsomeEnd}{Gene Start}
#'   \item{RGB}{RGB codes}
#' }
"Homo_sapiens.Hg19_BN_bed.txt"
