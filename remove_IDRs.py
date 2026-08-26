import argparse

def filter_pdb_by_confidence(input_pdb, output_pdb, threshold):
    atom_count = 1  # Initialize atom counter
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    confidence = float(line[61:66].strip())  # Extract B-factor/confidence score (fixed index)
                    if confidence >= threshold:
                        new_line = line[:6] + f"{atom_count:5d}" + line[11:]  # Renumber atom serial
                        outfile.write(new_line)
                        atom_count += 1
                except ValueError:
                    outfile.write(line)  # Keep if parsing fails
            else:
                outfile.write(line)  # Keep non-atom lines like headers and remarks

# Manually define input values for interactive testing
input_pdb = "spe18AF.pdb"  # Replace with file name
output_pdb = "spe18noIDR.pdb"
threshold = 50.0  # Set confidence threshold

filter_pdb_by_confidence(input_pdb, output_pdb, threshold)
print(f"Filtered PDB saved to {output_pdb}")

