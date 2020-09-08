# Utility: blob of function frequently used in my research

- Author: Bansho Masutani<banmasutani@gmail.com>
- Language: Rust

## Current implimentation

+ get_sam: parse sam file and return a HashMap with read id as its key and (flag,mapped location) as its value.
+ setup_template_complement: parse a given fasta file and return the sequence and its revcomp.
+ setup_template_complement_autofit: The same as above exept that its lengths are restricted by the second argument.
+ extend_strand: extend a given string untils its length reaches a specified length.
+ convert_to_squiggle: convert a string to its corresponds event stream by using Squiggler model(for more detail, see Squiggler crate).
+ get_mode: parse argument and returns its corresponds DTW mode.
+ get_queries: Extract event streams from files in specified folder.
+ merge_queries_and_sam: Marge given sam file and queries by using read id.
+ get_dataset: compute DTW metric for each query and generate dataset.
+ get_score: wrapper function of DTW.
+ cross_validation: wrapper function of get_dataset and get_score.


## Requirements

- rust > 1.2 or later