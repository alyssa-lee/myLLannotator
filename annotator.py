#!/usr/bin/env python

"""
annotator.py by Alyssa Lee and Rohan Maddamsetti.

Usage: Run this script in parallel in separate terminal windows as follows:
python annotator.py 1
python annotator.py 2
""" 

import ollama
import numpy as np
import pandas as pd
import re
from tqdm import tqdm ## for progress bar
import argparse


def annotate_genomes_with_ollama(model_name, output_path, TEST_MODE=True):

    ## Input files for annotation.
    df = pd.read_csv("../data/Maddamsetti2024/gbk-annotation-table.csv", na_filter=False)
    ground_truth = pd.read_csv("../data/Maddamsetti2024/computationally-annotated-gbk-annotation-table.csv", na_filter=False)
    
    ## configuration variables
    N_EXAMPLE = 10
    DEBUGGING = False

    if TEST_MODE: ## then just annotate a few genomes.
        df = df.head(n=N_EXAMPLE)
        ground_truth = ground_truth.head(n=N_EXAMPLE)

    nrow = df.shape[0]

    annotations = np.full(nrow, "NoAnnotation")
    categories = [
        "Agriculture",
        "Animals",
        "Food",
        "Freshwater",
        "Human-impacted",
        "Humans",
        "Livestock",
        "Marine",
        "Plants",
        "Sediment",
        "Soil",
        "Terrestrial",
        "NA"]

    categories_str = ', '.join(['"' + category + '"' for category in categories])

    for i in tqdm(range(nrow)):
        input_host = df.iloc[i]['host']
        input_isolation_source = df.iloc[i]['isolation_source']

        cur_row_annotated = False
        ## retry until the model does it right.
        while not cur_row_annotated: 

            system_prompt = "You are an annotation tool for labeling the environment category that a microbial sample came from, given the host and isolation source metadata reported for this genome. Label the sample as one of the following categories: " + categories_str + " by following the following criteria. Samples from a human body should be labeled 'Humans'. Samples from domesticated or farm animals should be labeled 'Livestock'. Samples from food should be labeled 'Food'. Samples from freshwater should be labeled 'Freshwater'.Samples from a human-impacted environment or an anthropogenic environmental source should be labeled 'Human-impacted'. Samples from the ocean, including the deep ocean should be labeled 'Marine'. Samples from anoxic sediments, including aquatic sediments should be labeled 'Sediment'. Samples from domesticated plants and crops should be labeled 'Agriculture'. Samples from soil, including farm soil should be labeled 'Soil'. Samples from extreme terrestrial environments (extreme temperature, pH, or salinity) should be labeled 'Terrestrial'. Samples from non-domesticated or wild plants should be labeled 'Plants'. Samples from non-domesticated or wild animals, including invertebrates, and also including protists and fungi even though these are not strictly animals should be labeled 'Animals'. Samples that lack enough information in the host metadata and isolation source metadata provided or have missing or incomplete information in these fields should be labeled 'NA'. Give a strictly one-word response that exactly matches of these categories, omitting punctuation marks."
            
            LLMresponse  = ollama.chat(
            model=model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": "Consider a microbial sample from the host " + input_host + " and the isolation source " + input_isolation_source + ". Label the sample as one of the following categories: " + categories_str + ". Give a strictly one-word response without punctuation marks."},
            ]
        )

            ## check if the LLM annotation matches one of the given categories for annotation.
            if LLMresponse.message.content in categories:
                annotations[i] = LLMresponse.message.content
                cur_row_annotated = True
            
            if DEBUGGING:
                ## check the output:
                print("INPUT FIELDS")
                print(input_host, input_isolation_source)
                print("PRINTING OUTPUT:")
                print(LLMresponse.message.content)
            
    ## save the results
    df['Annotation'] = annotations
    df.to_csv(output_path, index=False)

    n_correct = sum(df['Annotation'] == ground_truth['Annotation'])
    print(model_name + " accuracy: " + str(n_correct/nrow) + " (" + str(n_correct) + " out of " + str(nrow) + ")")
    return


def main():
    parser = argparse.ArgumentParser(description="My LiL Annotator :-).")
    parser.add_argument(
        "model",
        type=int,
        help="A required integer argument."
    )
    args = parser.parse_args()

    ## I chose these two models because they run rather quickly.
    if args.model == 1:
        annotate_genomes_with_ollama("llama3.2:latest", "../results/llama3.2_latest_gbk-annotation-table.csv", TEST_MODE=False)
    elif args.model == 2:
        annotate_genomes_with_ollama("gpt-oss:120b-cloud", "../results/gpt-oss_120b-cloud_gbk-annotation-table.csv", TEST_MODE=False)
    quit()
    
    return


if __name__ == "__main__":
    main()
