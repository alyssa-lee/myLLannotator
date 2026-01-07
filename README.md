# myLLannotator

Alyssa Lee and Rohan Maddamsetti

![A muscular cyborg rainbow llama with a face stripe like ziggy stardust and a long rainbow mane working hard on a laptop](img/llama-1.jpg)


## Overview

User-friendly tool for automated annotation of metadata with open-source LLM

## Requirements

- python>=3.12
- numpy>=2.4.0
- ollama>=0.6.1
- pandas>=2.3.3
- tqdm>=4.67.1

To reproduce paper figures:
- R==4.2+ for generating figures and re-running analyses in this paper

## Installation

There are two options for using this code. The first way is to install the prebuilt package, which should automatically install dependencies. The second way is to download the script `src/myllannotator/main.py` and run it directly (but you will have to install dependencies yourself). (If you are doing this, skip this section.)

The package can be installed from PyPI or from source (tarball).

### How to install as a package from PyPI:

This is the test version on TestPyPI.
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple myllannotator
```
(Once it is published to PyPI, the command `pip install myllannotator` should work)

### How to install as a package from source:

First download the compressed binary from the latest [release](https://github.com/alyssa-lee/myLLannotator/tags). Then run:
```
pip install myllannotator-0.1.1.tar.gz
```


## Downloading the llama3.2 model
Download the llama3.2 model from ollama, after installing the ollama package (This is a required step; working on getting it to download automatically):
```
ollama run llama3.2:latest
```

## How to run

The following sample commands use the example data under `input/` and write the output to a new file `annotated_data.csv`.

How to run with package installed:
```
myllannotator input/valid_categories.txt input/system_prompt.txt input/per_sample_prompt.txt input/input_data.csv annotated_data.csv
```

How to run without package installed (assuming all dependencies are installed):
```
python main.py input/valid_categories.txt input/system_prompt.txt input/per_sample_prompt.txt input/input_data.csv annotated_data.csv
```

## Usage
Brief overview of command line usage:
```
usage: myllannotator [-h]
                     valid_categories system_prompt per_sample_prompt
                     input_csv output_csv

positional arguments:
  valid_categories   .txt file of valid categories, separated by line breaks.
  system_prompt      .txt file containing system prompt
  per_sample_prompt  .txt file containing per-sample prompt
  input_csv          .csv file of input data
  output_csv         .csv file for output data

options:
  -h, --help         show this help message and exit
```

Also see `input/` for examples of each input format.

### valid_categories (.txt)
List of categories separated by line breaks. Make sure to include an NA category if you want the model to have the option to assign no annotation.

Example:
```
Human
Animal
NA
```

### system_prompt (.txt)
The system prompt guides the LLM's overall behavior. Here you should give specific instructions for the annotation task.

Optionally, you can include `{categories}` somewhere in the text, which will be replaced by a comma-separated list of the values in `valid_categories` (for example, `"Human", "Animal", "NA"`). The tool will print the properly formatted version upon running so you can check if it is what you expected.

Example:
> ```You are an annotation tool for labeling the environment category that a microbial sample came from, given the host and isolation source metadata reported for this genome. Label the sample as one of the following categories: {categories} by following the following criteria. Samples from a human body should be labeled 'Humans'. Samples from domesticated or farm animals [...] Give a strictly one-word response that exactly matches of these categories, omitting punctuation marks.```

### per_sample_prompt (.txt)
The per-sample prompt tells the LLM the relevant metadata for each sample.

The way you write this prompt **will depend on the columns in your input data**. Where you write `{0}` in the prompt, it will be replaced by the value in column 0 of the input data, `{1}` will be replaced by the value in column 1, etc. See the example below. The tool will print the properly formatted version upon running so you can check if it is what you expected.

Optionally, you can include `{categories}` somewhere in the text, which will be replaced by a comma-separated list of the values in `valid_categories` (for example, `"Human", "Animal", "NA"`). The tool will print the properly formatted version upon running so you can check if it is what you expected.

Example (for an input dataset with three columns `Annotation_Accession,host,isolation_source`):
> ```Consider a microbial sample from the host "{1}" and the isolation source "{2}". Label the sample as one of the following categories: {categories}. Give a strictly one-word response without punctuation marks.```

For the first sample, the prompt received by the LLM will be:
> ```Consider a microbial sample from the host "chicken" and the isolation source "Epidemic materials". Label the sample as one of the following categories: "Humans", "Livestock", "Food", "Freshwater", "Anthropogenic", "Marine", "Sediment", "Agriculture", "Soil", "Terrestrial", "Plants", "Animals", "NA". Give a strictly one-word response without punctuation marks.```


### input_csv (.csv)

This is your input data. It can have any number of columns. You will need to write your per-sample prompt according to the column order (see above).

```
Annotation_Accession,host,isolation_source
GCF_019552145.1_ASM1955214v1,chicken,Epidemic materials
GCF_001635975.1_ASM163597v1,Homo sapiens,NA
GCF_900636445.1_41965_G01,NA,Oral Cavity
```

### output_csv (.csv)
Path to the output file, which will be created as the program runs.

The format of the output file will be the same as the input file, with one additional column for the annotation.

Example:
```
Annotation_Accession,host,isolation_source,Annotation
GCF_019552145.1_ASM1955214v1,chicken,Epidemic materials,Livestock
GCF_001635975.1_ASM163597v1,Homo sapiens,NA,Humans
GCF_900636445.1_41965_G01,NA,Oral Cavity,Humans
```


## Important notes on the behavior of myLLannotator
- The tool is not deterministic. Different answers may be produced on the same input data.
- The tool will give up on labeling a particular sample after 5 failed attempts. In that case `NoAnnotation` will show up as the annotation.


## Replicating results in the paper

1. **Create a project directory**:
   ```
   mkdir myLLannotator/
   cd myLLannotator/
   mkdir results
   ```
   
   2. **Download data into the project directory**
   Go to https://zenodo.org/records/18110824 and save the download as a directory `data/` in the project directory.
   3. **Download this github repository**
   Save the download as a directory `src/` in the project directory.
   4. Run the R script `src/simple-ARG-duplication-analysis.R` to generate figures and run the analysis.
   
