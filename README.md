# myLLannotator

Alyssa Lee and Rohan Maddamsetti

![A muscular cyborg rainbow llama with a face stripe like ziggy stardust and a long rainbow mane working hard on a laptop](img/llama-1.jpg)
## Overview

details go here

## Requirements

- Python==3.12.9+
- DETAILS OF PYTHON LIBRARIES ETC GO HERE.
- R==4.2+ for generating figures and re-running analyses in this paper

## Usage

How to run (without package built):

uv run main.py input/valid_categories.txt input/system_prompt.txt input/per_sample_prompt.txt input/input_data.csv output/annotated_data.csv


### Replicating results in the paper

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
   
