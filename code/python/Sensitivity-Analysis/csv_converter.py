import re
import pandas as pd

INPUT_TXT = "./SA_results/temp9-no-r.txt"
OUTPUT_CSV = "sobol_results.csv"

rows = []
current_n = None

with open(INPUT_TXT) as f:
    for line in f:
        line = line.strip()

        # block header like "2^4: ----" sets the sample size for rows that follow
        header = re.match(r"^2\^(\d+):", line)
        if header:
            current_n = int(header.group(1))
            continue

        # skip blank lines and repeated column-header row
        if not line or line.startswith("Params"):
            continue

        tokens = line.split()
        if len(tokens) != 6:
            continue

        param, s1, st, st_conf, high_acc, below_cutoff = tokens
        rows.append({
            "Params": param,
            "Sobol First Order (S1)": float(s1),
            "Sobol Total (ST)": float(st),
            "Sobol Total Conf": float(st_conf),
            "High Accuracy": high_acc,
            "Below Relative Cutoff": below_cutoff,
            "N": current_n,
        })

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)
print(f"saved {len(df)} rows to {OUTPUT_CSV}")
