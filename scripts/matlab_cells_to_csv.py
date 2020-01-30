from scipy.io import loadmat
import pandas as pd
import argparse
import os
from pathlib import Path

def convert_file(source, destination=""):
    matlab_file = loadmat(source)
    key = ""
    if source.endswith("cells_physicell.mat"):
        key = "cells"
    elif source.endswith("cells.mat"):
        key = "basic_agents"
    elif source.endswith("microenvironment0.mat"):
        key = "multiscale_microenvironment"
    elif source.endswith("mesh0.mat"):
        key = "mesh"
    data = pd.DataFrame(matlab_file[key])
    data = data.transpose()
    data.to_csv(os.path.join(os.path.dirname(source), destination, "%s.csv" %os.path.splitext(os.path.basename(source))[0]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to change from matlab to csv")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        required=True,
                        help="Input file/folder containing matlab files")

    options = parser.parse_args()
    source = os.path.abspath(options.input)
    if os.path.isfile(source):
        convert_file(source)
    else:
        Path(source+"/csv/").mkdir(parents=True, exist_ok=True)
        for file in os.listdir(source):
            if file.endswith(".mat"):
                file = os.path.join(source, file)
                convert_file(file, "csv")