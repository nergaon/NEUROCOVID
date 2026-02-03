from pathlib import Path

def main():
    dirs = ["Cortical_neurons_microglia", "Cortical_neurons", "Microglia", "organoids"]
    base_dir = Path("/gpfs0/tals/projects/Analysis/NEUROCOVID/hadas_harschnitz/bulkRNA/leafcutter_0.2.9")
    input_file = base_dir / "NEUROCOVID_perind_numers.counts"
    with open(input_file) as f:
        first_line = f.readline().strip()
    samples = [s.strip() for s in first_line.split(" ")]
    for d in dirs:
        out_dir = base_dir / d
        out_dir.mkdir(exist_ok=True)
        out_file = out_dir / "group_file.txt"
        with open(out_file, "w") as f:
            for s in samples:
                # exact match: directory + '_' + condition
                #print(d, s)
                if s.startswith(d + "_"):
                    if "_mock_" in s:
                        f.write(f"{s} mock\n")
                    elif "_SARS-CoV_" in s:
                        f.write(f"{s} SARS-CoV\n")
    return

if __name__ == "__main__":
    main() 