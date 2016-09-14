"""A script to divide a table based on column names. """

import pandas as pd
import argparse
import sys

def main(args):
    df = pd.read_table(args.input_table, index_col=0)
    
    column_names = []
    with open(args.column_names) as cn_fh:
        for line in cn_fh.readlines():
            column_names.append(line.strip())

    assert len(column_names) != 0

    for output_col in column_names:
        try:
            assert output_col in df.columns
        except AssertionError:
            print(output_col)
            raise
    
    df[column_names].to_csv(sys.stdout, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("input_table", help="The input table")
    parser.add_argument("column_names", help="A file containint the column names that should be present in the output table")
    args = parser.parse_args()
    main(args)

