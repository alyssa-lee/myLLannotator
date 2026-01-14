import argparse
import ollama
from tqdm import tqdm
import pandas as pd

def annotate_samples(args):
    if args.debug:
        print("Settings: ")
        print(args)

    df = pd.read_csv(args.input_csv, na_filter=False)

    with open(args.valid_categories) as f:
        categories = [line.rstrip() for line in f]

    categories_str = ', '.join(['"' + category + '"' for category in categories])

    with open(args.system_prompt) as f:
        system_prompt = f.read().strip().format(categories=categories_str)
    with open(args.per_sample_prompt) as f:
        per_sample_prompt = f.read().strip()

    if not args.silent:
        print("System prompt:")
        print(system_prompt)
        print("\n")

        print("Example prompt using the first line:")
        print(per_sample_prompt.format(*df.loc[0, :].values.flatten().tolist(), categories=categories_str))
        print("\n")

    nrow = df.shape[0]
    if args.debug and nrow > 5:
        nrow = 5

    annotations = ["NoAnnotation" for i in range(nrow)]

    with open(args.output_csv, "w") as f:
        f.write(','.join([*list(df), "Annotation"]) + '\n')

    for i in tqdm(range(nrow)):
        cur_row_annotated = False

        row_values = df.loc[i, :].values.flatten().tolist()
        tries = 0

        ## retry until the model does it right
        while not cur_row_annotated and tries < args.max_tries:
            LLMresponse = ollama.chat(
                model=args.model_name,
                messages=[
                    {"role": "user" if args.disable_system_role else "system", "content": system_prompt},
                    {"role": "user", "content": per_sample_prompt.format(*row_values, categories=categories_str)},
                ]
            )
                
            ## check if the LLM annotation matches one of the given categories for annotation.
            if LLMresponse.message.content.strip() in categories:
                annotations[i] = LLMresponse.message.content.strip()
                cur_row_annotated = True

            tries = tries + 1

            if args.debug:
                print("Prompt: ")
                print(per_sample_prompt.format(*row_values, categories=categories_str))
                print("Response: ")
                print(LLMresponse.message.content)
                print("Tries: ")
                print(tries)

        with open(args.output_csv, "a") as f:
            f.write(','.join([*row_values, annotations[i]]) + '\n')

    return

def main():

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

    
    parser.add_argument(
        "--model-name",
        type=str,
        help="ollama model name, default is llama3.2:latest",
        default="llama3.2:latest"
    )
    parser.add_argument(
        "--max-tries",
        type=int,
        help="maximum number of attempts per sample if the LLM response is invalid, default is 5",
        default=5
    )
    parser.add_argument(
        "--silent",
        action="store_true",
        help="if enabled, do not print usual prompt output"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="if enabled, print debug output, and only annotate the first 5 samples"
    )
    parser.add_argument(
        "--disable-system-role",
        action="store_true",
        help="Disables the system role, instead having the system prompt come from the user. Set this option when using LLMs that do not have a system role."
    )

    args = parser.parse_args()

    annotate_samples(args)

    return


if __name__ == "__main__":
    main()
