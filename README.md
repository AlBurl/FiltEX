# FiltEX

**FiltEX filters and formats whole exome sequencing data, with optional application of custom gene panels and Exomiser for variant prioritisation.**

**Please note that FiltEX is currently under development.**



## Introduction

FiltEX filters, sorts, and formats ANNOVAR CSV files for facile analysis and interpretation of whole exome sequencing results. In addition to standard filtering, FiltEX can perform custom filtering of variants to a custom gene list. The Exomiser tool can also be run within FiltEX, allowing Exomiser variant prioritisation results to be formatted and combined into the results file. The final output is an Excel workbook, allowing filtered and formatted lists to be easily browsed and interpreted.

## Requirements

FiltEX uses the following software:
1.	[RStudio](https://www.rstudio.com/products/rstudio/)
2.	[Exomiser](https://github.com/exomiser/Exomiser)

**NB: Exomiser is licensed under the GNU Affero General Public License v3.0**

FiltEX uses the following R packages:
1.	Dplyr
2.	Readr
3.	Stringr
4.	Tidyr
5.	Writexl

## Functionality

FiltEX performs the following functions:
1.	wANNOVAR CSV file imported  
2.	Columns removed to tidy data 
3.	Synonymous variants filtered out
4.	Variants with population frequency greater than the defined frequency (default 1%) in 1000G, ExAC and gnomAD databases filtered out (generates ***All rare*** variant list)
5.	Gene reference sequence imported from online HGNC database, and used to select correct transcript and extract corresponding AA variant and DNA variants into separate columns 
6.	***All rare*** variants split into heterozygous and homozygous lists
7.	Heterozygous variant list split into those that appear once in any given gene (single) and those that appear more than once in any given gene (multiple). 
8.	OPTIONAL: ***All rare*** variants filtered using custom panel gene list, and split into separate list
9.	OPTIONAL: Exomiser is run on raw vcf file using yml file generated by the user separately. Exomiser results automatically imported, and variants in prioritised genes selected from the raw csv file and separated into new ‘Exomiser’ list.
10.	All lists exported into Excel file in output folder, one list per sheet. FiltEX stats table is included, summarising the number of variants in each list.

## Setup 

1.	Create a master folder for FiltEX use. Within this folder, create ‘input’ and ‘output’ folders. 
2.	Pull the latest version of the FiltEX script from GitHub.
3.	FiltEX takes input of a wANNOVAR csv file. Before use, go to [wANNOVAR](http://wannovar.wglab.org), enter the VCF and download the resulting CSV file to the ‘input’ folder.
4.	Install packages if necessary using install.packages()
5.	Within FiltEX, enter the master folder location as indicated
6.	Edit options as desired 

## Optional: Custom gene panel setup

1.	Create the list of genes in csv format. FiltEX uses exact match of gene names to filter, so if the gene has alternative names/notation be sure to include these too. 
2.	Save the gene list into the master folder. 
3.	Edit the options in FiltEX to ‘Y’ for custom gene panel use and enter the file name of the gene panel. 

**NB: FiltEX is set up to filter the custom gene panel from the *All rare* list, not the raw data. Hence any variants filtered out in earlier stages will not be included in the custom filtering.** 

## Optional: Exomiser setup

1.	Install [Exomiser](https://github.com/exomiser/Exomiser)
2.	Test Exomiser in the command line and ensure it works
3.	Edit the options in FilEX to ‘Y’ for Exomiser, and enter working directory path (the same as used in the command line)
4.	Ensure that the Exomiser command is correct in the Exomiser section of script. This may need adapting depending on the set up of your Exomiser.

**NB: The use of Exomiser through FiltEX has only been tested through Windows currently. The code may need slight adaptation to run through Mac/Linux.**  

## Usage

1.	If starting with a VCF file, use [wANNOVAR](http://wannovar.wglab.org) to annotate the file, which generates a csv file. 
2.	Ensure the wANNOVAR file is in the ‘input’ folder
3.	Set up options as desired (see above) 
4.	Run whole script 
5.	Enter wANNOVAR file name as prompted in console (**without ‘.csv’**)
6.	If using Exomiser, enter yml file name as prompted in console (**without ‘.yml’**)
7.	When the script completes the results Excel file can be found in the output folder. 

## Test data 

A test wANNOVAR whole exome sequencing csv file has been included. This data was obtained in vcf format from Exomiser download, and was run through wANNOVAR to generate the csv.  

## Contact 

Alice Burleigh

alice.burleigh.19@ucl.ac.uk

University College London Great Ormond Street Institute of Child Health

## References 

If you use FiltEX, please reference Exomiser and ANNOVAR accordingly:

Smedley, Damian, Julius O. B. Jacobsen, Marten Jäger, Sebastian Köhler, Manuel Holtgrewe, Max Schubach, Enrico Siragusa, Tomasz Zemojtel, Orion J. Buske, Nicole L. Washington, William P. Bone, Melissa A. Haendel, and Peter N. Robinson. 2015. 'Next-generation diagnostics and disease-gene discovery with the Exomiser', Nature Protocols, 10: 2004-15.

Wang, K., M. Li, and H. Hakonarson. 2010. 'ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data', Nucleic Acids Research, 38: e164-e64.








