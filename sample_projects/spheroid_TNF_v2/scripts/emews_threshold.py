import argparse
import os
import copy
import xml.etree.ElementTree as ET
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to generate physicell settings")
    parser.add_argument('-t', '--template',
                        dest="template",
                        action="store",
                        required=True,
                        help="Template to use on creating xml files.")
    parser.add_argument('-o', '--output',
                        dest="output",
                        action="store",
                        required=True,
                        help="Folder to use to save generated xml files.")
    parser.add_argument('--min-range',
                        dest="min",
                        action="store",
                        required=True,
                        type=float,
                        help="Minimum range of TNF threshold")
    parser.add_argument('--max-range',
                        dest="max",
                        action="store",
                        required=True,
                        type=float,
                        help="Maximum range of TNF threshold")
    parser.add_argument('--step-range',
                        dest="step",
                        action="store",
                        required=True,
                        type=float,
                        help="Range step of threshold")

    options = parser.parse_args()
    template_source = os.path.abspath(options.template)
    destination_folder = os.path.abspath(options.output)

    tnf_thresholds = np.arange(options.min, options.max, options.step)

    xml_tree = ET.parse(template_source)

    xml_filenames = []

    for tnf_threshold in tnf_thresholds:
        tnf_threshold = round(tnf_threshold, 3)
        new_xml = copy.deepcopy(xml_tree)
        root = new_xml.getroot()
        root.find("user_parameters").find("tnf_threshold").text = str(tnf_threshold)
        
        new_xml.write(destination_folder+"/PhysiCell_settings_"+str(tnf_threshold)+".xml")

        xml_filenames.append("./config/PhysiCell_settings_"+str(tnf_threshold)+".xml")
    
    # write xml_filenames into a txt file
    with open('input.txt', 'w') as f:
        for item in xml_filenames:
            f.write("%s\n" % item)