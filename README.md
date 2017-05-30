
# MSFractionMerge

## Philosophy

This software merges fractionated mass spectrometry abundance profiles from
replicate experiments where fraction data is not immediately comparable between
replicates and must first be aligned to combine equivalent fractions.

Merging is achieved by progressive, pairwise aggregation of abundance profiles,
using a scheme that maximises the sum of spectral count Pearson correlation
coefficients of components (typically proteins) shared between aligned
fractions. Fraction alignment involves an exhaustive search of relative offset
and linear width scaling.

## Installation

MSFractionMerge is written entirely in Python and requires no special
installation. The executable "ms_fraction_merge" shell script may be run
directly at the location where it is saved or placed on the system's executable
path. Note that this should remain alongside the Python file
"ms_fraction_merge.py"

### Python Module Requirements

This software uses Python version 2 or 3 and requires that the NumPy is
installed and available to the Python version that runs ms_fraction_merge.

The NumPy module is available in bundled Python packages like Anaconda or
Canopy, in most Linux distributions' package managers or can be installed on
most UNIX-like systems using pip:

        pip install numpy
  
## Operation

MSFractionMerge may be run on the command line, specifying the input file names,
output file name and other options after the script name:

        ./ms_fraction_merge input_spec_counts_*.csv merged_replicates.csv -n 5

See the command-line summary below for a detailed listing of all available options.
Test data is provided in the "example_data" directory.


## File Formats

### Input 

The input spectrum count data should be in tab or comma separated text files and
there should be one file for each replicate experiment. The input files must
have an initial header line, containing the fraction/column numbers, and
subsequent data lines, containing the spectral counts.

The first field in the header line is ignored and the subsequent fields are
expected to be integer fraction numbers corresponding the column data below.
Fraction numbers need not be ordered or sequential; fractions may be missing.

Data lines must begin with the component (e.g. protein) ID and then be followed
by the spectrum count for the fraction columns, as specified in the header.
Missing values must be represented by a zero. 

For example the header and first two data lines of a file may look like:
  
        NULL,13,14,15,18,19,20,22,23,24,25,26,27,28,29,31,33,34,35,37,40,41,43,44,46,49,52,56
        AT5G28540.1,29,22,25,13,18,28,34,148,140,163,176,141,125,168,135,163,258,260,209,483,379,29,0,3,0,0,0
        AT5G42020.1,30,22,25,13,14,25,35,138,134,161,168,122,111,152,115,142,233,231,191,456,353,30,0,0,0,0,0

The optional component category CSV files should be of tab or comma separated
format, without a header line, and contain two columns corresponding to the
component (e.g. protein) ID and the category name for that component. The
component ID should match the input spectrum count data. For example:

        AT5G52040.3,cytoplasm
        AT5G64770.1,ER
        AT5G67480.1,vacuole

### Output

The output file will be a comma-separated text file containing an initial  header
line and then subsequent data lines. Each data line consists of:

 * A component ID, as given in the input data
 * A component tag/label, as specified in the optional category file. This value defaults to "unknown"
 * The number of input (replicate) data sets that had data for this protein
 * Further lines containing merged, relative protein abundances for each pseudo-fraction
 
For example the first three lines may look something like:
 
        component,category,num_reps,0,1,2,3,4
        AT1G04710.1,PRX,3,0.0000000,0.9230769,1.8963752,1.1573064,1.1357401
        AT1G04910.1,GLG,5,0.4381594,0.6021773,0.9947047,1.0191765,1.0079090


---

## Command-line Summary

usage:

    ms_fraction_merge [-h] [-o OUT_CSV_FILE]
                           [-c CATEGORY_CSV_FILE [CATEGORY_CSV_FILE ...]]
                           [-n NUM_FRACTIONS] [-m MIN_REPS] [-omax MAX_OFFSET]
                           [-ostep STEP] [-smax MAX_SCALE] [-sstep STEP]
                           [-clip FRACTION_NUMS [FRACTION_NUMS ...]] [-q] [-log]
                           CSV_FILES [CSV_FILES ...]

positional arguments:

    CSV_FILES             Input file paths of CSV files containing spectral
                          count data, one for each replicate (may contain
                          wildcards)

optional arguments:

    -h, --help            show this help message and exit
    -o OUT_CSV_FILE       Optional output file path for CSV file containing
                          merged spectral count data; Default
                          merged_replicates.csv
    -c CATEGORY_CSV_FILE [CATEGORY_CSV_FILE ...]
                          One or more optional CSV files describing
                          categories/classes, e.g. for marker proteins, relating
                          to the input data (may contain wildcards)
    -n NUM_FRACTIONS      Number of output (pseudo-)fraction columns to output
                          for the merged data. Default: 25
    -m MIN_REPS           Minimum number of replicate occurrences required for a
                          protein to be accepted. Default: 3
    -omax MAX_OFFSET      The maximum number of whole columns that could offset
                          the first values between different replicate data;
                          defines the width of the offset parameter search
                          space. Default: 4.00
    -ostep STEP           The fractional increment in column widths to use when
                          aligning replicate data; defines the granularity of the
                          offset parameter search space. Default: 0.50
    -smax MAX_SCALE       The maximum stretch/compression scale factor between
                          different replicate data; defines the width of the
                          scale parameter search space. Default: 2.00
    -sstep STEP           The increment in stretch scale factor use when
                          aligning replicate data; defines the granularity of the
                          scale parameter search space. Default: 0.20
    -clip FRACTION_NUMS [FRACTION_NUMS ...]
                          Optional fraction numbers to clip input data at (one
                          for each replicate in input order), i.e. beyond which
                          there is no useful data
    -q                    Sets quiet mode to suppress on-screen reporting.
    -log                  Log all reported output to a file.

For further help on running this program please email tjs23@cam.ac.uk.

Example use:

        ms_fraction_merge example_data/test_spec_count_data_*.csv -o merged_replicates.csv -c example_data/test_categories.tsv -m 2 -clip 50 47 37

