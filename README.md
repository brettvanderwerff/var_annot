# VAR_ANNOT
##### A simple VCF annotation tool.

## What does it do?

    
Parses a VCF file and outputs an annotated VCF file with the following 
information for each variant:

- 'TYPE' (The type of allele, either snp, mnp, ins, del, or complex)
- 'DP' (Depth of sequence coverage at the site of variation)
- 'AO' (Alternate allele observation count)
- 'AF' (Estimated allele frequency in the range (0,1))
- '#CHROM' (Chromosome that the variant is on)
- 'POS' (Position of the variant on the chromosome)
- 'REF' (Allele of the reference genome)
- 'ALT' (Alternative allele detected by the variant caller)
- 'ALLELE_FREQUENCY(ExAC)' (Allele frequency as reported by the ExAC Browser 
API (http://exac.hms.harvard.edu/))

## Setup

```
$git clone https://github.com/brettvanderwerff/var_annot 
$cd var_annot
$pip install -r requirements.txt
```

## Usage Case

var_annot has a simple command line interface and accepts two command 
line arguments

- path to the VCF file to be annotated
- the row number of the VCF header (usually the first row following the
 meta-information rows in the VCF file)

From the top level folder:

```
$python3 var_annot.py ./example_data/Challenge_data.vcf 142
```

Processing takes a few minutes.
When complete an annotated VCF should appear in the `output` folder


## Requirements

##### Python Version:

- Tested on 3.6.8

##### Python Packages:

- numpy==1.16.2
- pandas==0.24.2
- requests==2.21.0

##### OS:

- Tested on Ubuntu 18.04

## Future Directions:

I did not get the opportunity to implement all of the requested features.
I do have some thoughts on how to complete this project and improve 
the existing functionality.

1. It has a fragile CLI right now, but could be improved with one of 
the nice CLI packages like 'argparse'.

2. The `expand_df` function, which turns DataFrame columns that contain 
delimited key value pairs into their own DataFrames is a bottleneck. It could be
parallelized easily, but there are also a lot of loops in this function, which may need to be rethought.

3. Part of the challenge was to annotate:
 
     "Type of variation (Substitution, 
    Insertion, Silent, Intergenic, etc.) If there are multiple
    possibilities, annotate with the most deleterious possibility."
    
    I just did not get this far, but my logic was to have a DataFrame with
    multiple rows dedicated to a different one of these possibilities. I
    could then drop duplicates based upon a function that ranks the
    different deleterious possibilities. I would keep only the most deleterious after dropping duplicates.
    I was headed in this direction, that is why in some cases the final annotated VCF has
    multiple rows representing the same chromosome position.

4. The functions that call to the ExAC server are a little fragile.
It would be better if they handled exceptions from a
non-responsive server and gave informative readout to the user.

5. Of course having a "gold standard" VCF file and building tests around it would be good.

    




