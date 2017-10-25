# Index_Swp_DR
Repository for Collaborating on Data and Results for Index swapping presentation

## File list:

1. Index_Swp_sort.py - Python script that was used to generate Index hopped read counts, and used to sort reads before FastQC analysis
    
    a. Output sorted read files can be found at: hpc:/home/abubie/qual_ind_swp/pres/

2. Index_Sort_sc.srun - BASH script used to execute python script; includes passed arguements
3. FastQC_anal_sc.srun - BASH script used to run FastQC analysis on the sorted reads
4. FastQC_Output - Folder contains the output files from FastQC analysis 
5. Formatted_Sorted_Indexes - Table of all Index pairs, with number of reads per pair.
6. demultiplexed_tile_counts.txt - Text file of read counts per tile (UNIX output)
7. 20170802_fragment_analyzer_project_summary.pdf - Fragment analyzer traces from GC3F

8. Data_Analysis_Ideas - R-markdown file establishing some ideas for results figures (Meant as scratch file)
