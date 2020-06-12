import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import Counter
import numpy as np
import itertools
import copy

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

    return file_dict

def node_counter(number_of_files, file_dict, percentage):
    "create a dict with the count of the nodes in the network, it can be used to print percentage pie chart"
    count_dict = {}
    fixed_dict = {}
    max_cell = 0
    for i in range (0, number_of_files):
        node_list = []
        for key in file_dict["node_step{0}".format(i)]:
            for value in file_dict["node_step{0}".format(i)][key]:
                node_list.append(value)
        node_counts = Counter(node_list)
        max_cell = max_cell + sum(node_counts.values())
        #fix_count_dict = {}
        #for key, group in itertools.groupby(node_counts, lambda k: 'others' if (node_counts[k]<(0.01* (len(file_dict)))) else k):
            #fix_count_dict[key] = sum([node_counts[k] for k in list(group)])
        count_dict["node_count{0}".format(i)] = node_counts
    fixed_dict = filter_states(max_cell, count_dict, percentage)
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

def state_counter(number_of_files, file_dict, percentage):
    "create a dict with the states of the network, it can be used to print states pie chart"
    count_dict = {}
    fixed_dict = {}
    max_cell = 0
    for i in range (0, number_of_files):
        state_list = []
        for key in file_dict["state_step{0}".format(i)]:
            state_list.append(file_dict["state_step{0}".format(i)][key])
        state_counts = Counter(state_list)
        max_cell = max_cell + sum(state_counts.values())
        #fix_count_dict = {}
        #for key, group in itertools.groupby(state_counts, lambda k: 'others' if (state_counts[k]<(0.01* (len(file_dict)))) else k):
            #fix_count_dict[key] = sum([state_counts[k] for k in list(group)])
        count_dict["state_count{0}".format(i)] = state_counts
    fixed_dict = filter_states(max_cell, count_dict, percentage)    
    return fixed_dict

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
    fig1, ax1 = plt.subplots(figsize=(9, 8), subplot_kw=dict(aspect="equal"))
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




def create_area_chart(count_dict):
    "plot an area chart with the evolution of the network states during the simulation"

    state_list = []
    all_state = []
    a = []
    for k in count_dict:
        state_list.append(list(count_dict[k].keys()))
        for l in state_list:
            for state in l:
                all_state.append(state)
    all_state = list(dict.fromkeys(all_state))
    
    for state_count in count_dict:
        b = []
        for states in all_state:
            try:
                b.append(count_dict[state_count][states])
            except:
                b.append(0)
        a.append(b)
    a = np.array(a)
    #print(a)
    a = np.transpose(a)
    #percent = a /  a.sum(axis=0).astype(float) * 100
    percent = a
    x = np.arange(len(count_dict))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.stackplot(x, percent, labels=all_state)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax.set_title('100 % stacked area chart')
    ax.set_ylabel('Percent (%)')
    ax.margins(0, 0) # Set margins to avoid "whitespace"

    plt.show()

def filter_states(max_cell, all_counts, percentage):
    """max_cell = 0
    all_counts = {}
    for i in range (0, number_of_files):
        state_list = []
        for key in file_dict["state_step{0}".format(i)]:
            state_list.append(file_dict["state_step{0}".format(i)][key])
        state_counts = Counter(state_list)
        max_cell = max_cell + sum(state_counts.values())
        all_counts[i] = state_counts"""

    copy_all_counts = copy.deepcopy(all_counts)

    state_list = []
    all_state = []
    for k in all_counts:
        state_list.append(list(all_counts[k].keys()))
        for l in state_list:
            for state in l:
                all_state.append(state)
    all_state = list(dict.fromkeys(all_state))
    
    banned_list = []
    for state in all_state:
        a = 0
        for i in all_counts.keys():
            try: 
                a = a + all_counts[i][state]
            except:
                a = a + 0
        if (a < (percentage/100) * max_cell):
            banned_list.append(state)
            for i in all_counts.keys():
                del all_counts[i][state]
    
    for i in all_counts.keys():
        b = 0
        for state in banned_list:
            try:
                b = b + copy_all_counts[i][state]
            except:
                b = b + 0
        all_counts[i]["others"] = b

    return all_counts






