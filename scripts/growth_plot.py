import argparse
import os
import json
import matplotlib
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to change from matlab to csv")
    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        required=True,
                        help="Input folder containing matlab files")

    options = parser.parse_args()
    source = os.path.abspath(options.input)
    if os.path.isfile(source):
        with open(source, 'r') as f:
            data = json.load(f)
        time = list(map(int,data.keys()))
        num_cells = list(map(lambda x: x["alive"],data.values()))

        fig, ax = plt.subplots()
        ax.plot(time, num_cells)

        ax.set(xlabel='time (min)', ylabel='number of cells alive',
            title='Growth curve')
        ax.grid()
        plt.show()