from scipy.io import loadmat
import pandas as pd
import argparse
import os
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to change from matlab to csv")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        required=True,
                        help="Input folder containing matlab files")
    parser.add_argument('-t', '--time-step',
                        dest="time",
                        action="store",
                        type=int,
                        required=True,
                        help="Time step over files.")

    options = parser.parse_args()
    source = os.path.abspath(options.input)
    if os.path.isdir(source):
        result = {}
        time = 0
        for file in sorted(os.listdir(source)):
            if file.startswith("output") and file.endswith("cells_physicell.mat"):
                cells = {
                    "alive": 0,
                    "dead": 0,
                    "necrosis": 0,
                    "total": 0,    
                }
                matlab_file = loadmat(os.path.join(source, file))
                column_data = pd.Series(matlab_file["cells"][7])
                counts = column_data.value_counts()
                for phase, count in counts.items():
                    if phase < 100:
                        cells["alive"] += count
                    elif phase == 100:
                        cells["dead"] += count
                    elif phase > 100:
                        cells["necrosis"] += count
                    cells["total"] += count
                result[time] = cells
                time += options.time
        with open('count.json', 'w') as fp:
            json.dump(result, fp)
    else:
        pass # rise error