import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import Counter
import numpy as np

def create_dict(number_of_files, folder):
    "create a dictionary with the states file in the folder 'output', half of the dict is used to calculate the percentage of the node, the other half is for the states"
    file_dict = {}
    for i in range (0, number_of_files):
        nodes_dict = {}
        states_dict = {}
        with open('%sstates_%08u.csv' %(folder,i), newline='') as csvfile:
            has_header = csv.Sniffer().has_header(csvfile.read(1024))
            csvfile.seek(0)
            states_reader = csv.reader(csvfile, delimiter=',')
            if has_header:
                next(states_reader)
            for row in states_reader:
                states_dict[row[0]] = row[1]
                nodes_dict[row[0]] = row[1].replace("--", "").split()
        file_dict["node_step{0}".format(i)] = nodes_dict
        file_dict["state_step{0}".format(i)] = states_dict

    print(file_dict["state_step3"]["5"])
    return file_dict

def node_counter(number_of_files, file_dict):
    "create a dict with the count of the nodes in the network, it can be used to print percentage pie chart"
    count_dict = {}
    for i in range (0, number_of_files):
        node_list = []
        for key in file_dict["node_step{0}".format(i)]:
            for value in file_dict["node_step{0}".format(i)][key]:
                node_list.append(value)
        node_counts = Counter(node_list)
        count_dict["node_count{0}".format(i)] = node_counts
    return count_dict

def print_all_nodes_pie(node_counter_dict):
    "print a pie chart for each file in the dict, with the percentage of active nodes"
    for k in node_counter_dict:
        labels = node_counter_dict[k].keys()
        sizes = node_counter_dict[k].values()
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                shadow=True)
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.tight_layout()
        plt.show()

def state_counter(number_of_files, file_dict):
    "create a dict with the states of the network, it can be used to print states pie chart"
    count_dict = {}
    for i in range (0, number_of_files):
        state_list = []
        for key in file_dict["state_step{0}".format(i)]:
            for value in file_dict["state_step{0}".format(i)][key]:
                state_list.append(file_dict["state_step{0}".format(i)][key])
        state_counts = Counter(state_list)
        count_dict["state_count{0}".format(i)] = state_counts
    return count_dict

def print_all_states_pie(count_dict):
    "print a pie chart for each file in the dict, with the percentage of the network's states"
    for k in count_dict:
        labels = count_dict[k].keys()
        sizes = count_dict[k].values()
        fig1, ax1 = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

        wedges, texts, fig1 = ax1.pie(sizes, textprops=dict(color="w"), autopct='%1.1f%%')
        ax1.legend(wedges, labels=labels,
            title="Node States",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1))


        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.tight_layout()
        plt.show()

def print_states_pie(count_dict, step):

    "print a pie chart with the percentage of the network's states for a certain simulation step"

    k = "state_count{0}".format(step)
    labels = count_dict[k].keys()
    sizes = count_dict[k].values()
    fig1, ax1 = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

    wedges, texts, fig1 = ax1.pie(sizes, textprops=dict(color="w"), autopct='%1.1f%%')
    ax1.legend(wedges, labels=labels,
        title="Node States",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1))


    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.tight_layout()
    plt.show()


def print_nodes_pie(node_counter_dict, step):

    "print a pie chart with the percentage of active nodes for a certain simulation step"
    
    k = "node_count{0}".format(step)
    labels = node_counter_dict[k].keys()
    sizes = node_counter_dict[k].values()
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.tight_layout()
    plt.show()










