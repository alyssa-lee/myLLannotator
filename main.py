import argparse
import ollama
from tqdm import tqdm
import numpy as np
import pandas as pd

def annotate_samples(args):
    df = pd.read_csv(args.input_csv, na_filter=False)

    with open(args.valid_categories) as f:
        categories = [line.rstrip() for line in f]

    categories_str = ', '.join(['"' + category + '"' for category in categories])

    with open(args.system_prompt) as f:
        system_prompt = f.read().strip().format(categories=categories_str)
    with open(args.per_sample_prompt) as f:
        per_sample_prompt = f.read().strip()


    print("System prompt:")
    print(system_prompt)

    print("Example prompt using the first line:")
    print(per_sample_prompt.format(*df.loc[0, :].values.flatten().tolist(), categories=categories_str))


    model_name = "llama3.2:latest"

    nrow = df.shape[0]
    annotations = np.full(nrow, "NoAnnotation")

    with open(args.output_csv, "w") as f:
        f.write(','.join([*list(df), "Annotation"]) + '\n')

    for i in tqdm(range(nrow)):
        cur_row_annotated = False

        row_values = df.loc[i, :].values.flatten().tolist()
        tries = 0
        max_tries = 5

        ## retry until the model does it right
        while not cur_row_annotated and tries < max_tries:
            LLMresponse = ollama.chat(
                model=model_name,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": per_sample_prompt.format(*row_values, categories=categories_str)},
                ]
            )
                
            ## check if the LLM annotation matches one of the given categories for annotation.
            if LLMresponse.message.content in categories:
                annotations[i] = LLMresponse.message.content
                cur_row_annotated = True

            tries = tries + 1

        
        with open(args.output_csv, "a") as f:
            f.write(','.join([*row_values, annotations[i]]) + '\n')

    return

def main():
    print("Hello from myllannotator!")
    parser = argparse.ArgumentParser(description="myLLannotator")
    parser.add_argument(
        "valid_categories",
        type=str,
        help=".txt file of valid categories, separated by line breaks."
    )
    parser.add_argument(
        "system_prompt",
        type=str,
        help=".txt file containing system prompt"
    )
    parser.add_argument(
        "per_sample_prompt",
        type=str,
        help=".txt file containing per-sample prompt"
    )
    parser.add_argument(
        "input_csv",
        type=str,
        help=".csv file of input data"
    )
    parser.add_argument(
        "output_csv",
        type=str,
        help=".csv file for output data"
    )
    
    #TODO: optional output column name
    #TODO: model name
    #TODO: debug mode
    #TODO: max tries
    #TODO: record metrics (time, max tries)

    args = parser.parse_args()

    annotate_samples(args)

    return


if __name__ == "__main__":
    main()
