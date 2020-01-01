import argparse
import os
import sys


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', help='Path to files', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=True)
    parser.add_argument('--blacklist', help='Path to the blacklist file', required=False)
    parser.add_argument('--cpgi', help='Path to cpg island file', required=False)
    args = parser.parse_args()

    input_files = args.files
    if not os.path.exists(input_files):
        print("Couldn't find the input path")
        sys.exit(-1)

    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        print("Couldn't find the output path")
        sys.exit(-1)

    blacklist_path = args.blacklist
    cpgi_path = args.cpgi

    return input_files, output_folder, blacklist_path, cpgi_path


def remove_cpg_islands():
    pass


def remove_blacklist():
    pass


def main():
    remove_cpg_islands()
    remove_blacklist()


if __name__ == '__main__':
    main()
