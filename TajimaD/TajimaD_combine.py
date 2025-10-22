import os
import pandas as pd

input_dir = r"D:\software\ParaAT2.0\D\TajimaD_results_region"
output_file = r"D:\software\ParaAT2.0\D\TajimaD_summary_region.csv"

all_rows = []

for f in os.listdir(input_dir):
    if f.endswith(".txt"):
        filepath = os.path.join(input_dir, f)
        gene_name = f.replace("_TajimaD.txt","")
        vals = {}
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if "=" in line:
                    key, value = line.split("=")
                    key = key.strip()
                    value = value.strip()
                    try:
                        value = float(value)
                    except:
                        continue
                    vals[key] = value
        all_rows.append({
            "gene": gene_name,
            "pS": vals.get("pS", None),
            "pi": vals.get("pi", None),
            "S": vals.get("S", None),
            "Theta": vals.get("Theta", None),
            "Diff": vals.get("Diff", None),
            "SE": vals.get("SE", None),
            "D": vals.get("D", None)
        })

df = pd.DataFrame(all_rows)
df.to_csv(output_file, index=False)
print(f"汇总完成，保存到：{output_file}")
