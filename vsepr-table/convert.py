import csv
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert.py <input.csv> <output.h>")
        return

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    with open(in_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    with open(out_file, 'w') as f:
        f.write("#ifndef VSEPR_DATA_H\n")
        f.write("#define VSEPR_DATA_H\n\n")

        f.write("typedef struct {\n")
        f.write("    const char *valence_pairs;\n")
        f.write("    const char *ep_geometry;\n")
        f.write("    const char *bond_pairs;\n")
        f.write("    const char *lone_pairs;\n")
        f.write("    const char *shape;\n")
        f.write("    const char *hybridization;\n")
        f.write("} vsepr_row_t;\n\n")

        f.write(f"#define VSEPR_NUM_ROWS {len(rows)}\n\n")

        f.write("const vsepr_row_t vsepr_data[VSEPR_NUM_ROWS] = {\n")
        for r in rows:
            f.write(f"    {{\"{r[0]}\", \"{r[1]}\", \"{r[2]}\", \"{r[3]}\", \"{r[4]}\", \"{r[5]}\"}},\n")
        f.write("};\n\n")
        
        # also export headers to display them
        f.write("const char *vsepr_headers[6] = {\n")
        for h in header:
            f.write(f"    \"{h}\",\n")
        f.write("};\n\n")

        f.write("#endif\n")

if __name__ == "__main__":
    main()
