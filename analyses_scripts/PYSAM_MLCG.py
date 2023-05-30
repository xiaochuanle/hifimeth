# Import pysam and multiprocessing modules
import pysam
import sys
import re

# Define a function to find all CG motifs in a sequence and return their coordinates in a list
def find_CG(seq):
  # Initialize an empty list to store the coordinates
    coords = []
  # Use regular expression to find all matches of CG in the sequence
    matches = re.finditer("CG", seq)
  # Loop through the matches and append their start positions to the list
    for match in matches:
        coords.append(match.start())
  # Return the list of coordinates
    return coords

def process_record(record):
    # Get the ML tag value and the read sequence
    ml_tag = record.get_tag("ML")
    cg_coords = find_CG(record.query_sequence)
    
    if record.is_reverse:
        cg_coords = [record.query_length - i -2 for i in cg_coords][::-1]
        

    
    # Initialize a variable to store the current position in the read
    results=[]
    # Loop through the ML tag value and the read sequence simultaneously
    for ml, coord in zip(ml_tag, cg_coords):
        # If the ML value is 1, it means a CG motif is present
        # Increment the position by 1
        results.append([record.query_name, coord, ml])
    # Return the record name and the CG coordinates list as a tuple
    return results

# Define a function to write the results to a tsv file
def write_results(results, output_file):
    # Open the output file in write mode
    with open(output_file, "w") as f:
        # Loop through the results list
        for result in results:
            # Unpack the record name and the CG coordinates list
            for line in result:
                record_name, coord, ml = line
            # Write the record name to the file
                f.write(f'{record_name}\t{coord}\t{ml}')
            # Loop through the CG coordinates list
            # Write a newline character to the file
                f.write("\n")

# Create a pysam.AlignmentFile object to read the bam file
bamfile = pysam.AlignmentFile(sys.argv[1], "rb", check_sq=False)
RESULTS = []
# Create a multiprocessing.Pool object to use all available cores
for record in bamfile:
	    RESULTS.append(process_record(record))
# Apply the process_record function to each bam record using the pool.map method
write_results(RESULTS, sys.argv[2])


