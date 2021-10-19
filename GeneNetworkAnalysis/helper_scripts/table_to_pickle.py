import pandas as pd
import sys


def main():
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    df = pd.read_csv(input_file_path, sep="\t", index_col=0)
    print("table size:", df.shape)
    print("header:", df.columns[:5])
    print("rows:", df.index[:5])

    df.to_pickle(output_file_path)


if __name__ == '__main__':
    main()
