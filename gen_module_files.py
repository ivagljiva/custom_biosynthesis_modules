## program to convert custom modules from spreadsheet to modules file
## USAGE: ./gen_module_files.py Amino_Acids.txt|Vitamins.txt /PATH/TO/KEGG/ [OUTPUT_DIR]

import sys
import os
import re

from anvio import utils as utils

# expected params
input_file = sys.argv[1] if len(sys.argv) > 1 else None
kegg_dir_path = sys.argv[2] if len(sys.argv) > 2 else None
# optional params
output_dir = sys.argv[3] if len(sys.argv) > 3 else "./"
output_dir = os.path.join(output_dir, "modules")

if not input_file or not kegg_dir_path:
    print("Improper arguments provided. USAGE: ./gen_module_files.py [Amino_Acids.txt|Vitamins.txt] [/PATH/TO/KEGG] [OUTPUT_DIR]")
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ko_file_path = os.path.join(kegg_dir_path, "ko_list.txt")
if not os.path.exists(ko_file_path):
    print(f"No KO list file at {ko_file_path}")
    sys.exit(1)
ko_dict = utils.get_TAB_delimited_file_as_dictionary(ko_file_path, expected_fields = ['definition'])

lines = None
compound_type = None
with open(input_file, 'r', errors='ignore') as f:
    header = f.readline().split("\t")
    compound_type = header[1]
    lines = f.readlines()

current_module = None
target_compound = None
KOs_set = set([])
kegg_module_set = set([])
pmap_set = set([])
def_string = ""

def make_modules_file():
    """Prints info about the current custom module into a module file."""

    try:
        enzyme_dict = {k:{'source': "KOfam", 'orthology': ko_dict[k]['definition']} for k in KOs_set}
    except KeyError:
        print(f"The ko_list.txt file in the provided KEGG data directory doesn't contain the following KO: {k}")
    orthology_lines = []
    source_lines = []
    for acc,d in enzyme_dict.items():
        # first line gets the data name
        if len(orthology_lines) == 0:
            o_line_prefix = "ORTHOLOGY   "
            a_line_prefix = "ANNOTATION_SOURCE  "
        else:
            o_line_prefix = a_line_prefix = "            "
            a_line_prefix = "                   "

        ortho = d['orthology']
        src = d['source']

        orthology_lines.append(o_line_prefix + acc + "  " + ortho)
        source_lines.append(a_line_prefix + acc + "  " + src)

    all_lines = [f"ENTRY       {current_module}",
                 f"NAME        {target_compound} biosynthesis",
                 f"DEFINITION  {def_string}"] + \
              orthology_lines + \
                [f"CLASS       User modules; Biosynthesis; {compound_type} Biosynthesis",
                 f"PATHWAY     {';'.join(list(pmap_set))}"] + \
              source_lines + \
             ['///' + '\n']

    mod_file_path = os.path.join(output_dir, current_module)
    with open(mod_file_path, 'w') as f:
        f.write("\n".join(all_lines))
    
    return mod_file_path

for l in lines:
    if re.match("^\s+$", l): # empty line between different custom modules
        # print info to modules file
        of = make_modules_file()
        print(f"Printed {current_module} to output file path: {of}")
        # reset for next module
        current_module = None
        target_compound = None
        KOs_set = set([])
        kegg_module_set = set([])
        pmap_set = set([])
        def_string = ""
        continue

    fields = l.split("\t")

    if current_module == None and fields[0] != "":
        current_module = fields[0]
        target_compound = fields[1]
    
    if 'K' not in fields[3]:
        # we handle special connector strings like OR, AND
        if "AND" in fields[3]:
            connector = re.sub("( )?AND( )?", " ", fields[3]) # AND => space
        elif "OR" in fields[3]:
            connector = re.sub("( )?OR( )?", ",", fields[3]) # OR => comma
        elif fields[3] == "--":
            connector = " " + fields[3]
        else:
            connector = fields[3]
        def_string += connector
    else:
        ko_list = re.split('\(|\)|,|\+|-|"| ', fields[3])
        KOs_set.update(set([k for k in ko_list if k != ""]))
        if def_string == "" or def_string[-1] == '(' or def_string[-1] == ',' or def_string[-1] == ' ':
            def_string += fields[3].strip('"')
        else:
            def_string += " " + fields[3].strip('"')
    
    module_list = re.split('/', fields[5])
    kegg_module_set.update(set([m for m in module_list if m != ""]))
    map_list = re.split('/', fields[6])
    pmap_set.update(set([m for m in map_list if m != ""]))

# print info of final module to modules file
of = make_modules_file()
print(f"Printed {current_module} to output file path: {of}")
