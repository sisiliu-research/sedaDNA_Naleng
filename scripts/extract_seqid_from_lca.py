#!/usr/bin/env python3

# contact: sisi.liu@awi.de; sisi.liu.research@gmail.com

import os

# Define the list of input strings to search for
input_strings = ["1152:Pseudanabaena:genus", "1158:Oscillatoria:genus", "1177:Nostoc:genus", "3798:Saxifraga:genus", "4479:Poaceae:family", "5748:Nannochloropsis:genus", "7953:Cyprinidae:family", "8015:Salmonidae:family", "9789:Equus:genus", "9859:Cervus:genus", "9903:Bos:genus", "9977:Ochotona:genus", "13228:Potamogeton:genus", "13398:Carex:genus", "21880:Salvia:genus", "23204:Potentilla:genus", "24952:Myriophyllum:genus", "40685:Salix:genus", "43174:Pedicularis:genus", "47251:Leptolyngbya:genus", "54304:Planktothrix:genus", "102804:Asteroideae:subfamily", "167375:Cyanobium:genus", "202994:Rhodiola:genus", "217161:Chamaesiphon:genus", "337677:Cricetidae:family"]  # Add other strings as needed

# Define the directory path where your *.lca files are located
directory_path = "/path/to/lca"

# Iterate through all *.lca files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".lca"):
        # Parse the sample name from the input filename
        sample_name = filename.split(".")[0]

        # Create dictionaries to store lines for each input string
        lines_by_input_string = {string: [] for string in input_strings}

        # Open the *.lca file for reading
        with open(os.path.join(directory_path, filename), 'r') as lca_file:
            for line in lca_file:
                # Check if any of the input strings is present in the line
                for input_string in input_strings:
                    if input_string in line:
                        lines_by_input_string[input_string].append(line)

        # Process each input string
        for input_string, lines in lines_by_input_string.items():
            # Create a file to store lines for this input string
            output_filename = f"{input_string.replace(':', '_')}_{sample_name}.txt"
            with open(output_filename, 'w') as output_file:
                output_file.writelines(lines)

            # Extract the part of the 'HISEQ' ID before the seventh ':'
            seqid_filename = f"{input_string.replace(':', '_')}_{sample_name}_seqid.txt"
            with open(seqid_filename, 'w') as seqid_file:
                for line in lines:
                    parts = line.split(':')
                    if len(parts) >= 7:
                        seqid = ':'.join(parts[:7])
                        seqid_file.write(seqid + '\n')

            print(f"Extracted lines containing '{input_string}' in {sample_name} to {output_filename} and {seqid_filename}")

