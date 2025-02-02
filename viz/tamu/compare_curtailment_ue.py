import matplotlib.pyplot as plt
import pandas as pd
import sys
import json
import numpy as np
import toml

# Parse multiple simdir arguments
simdirs = sys.argv[1:]

# Check if the length of simdirs is 3 and assign hardcoded labels
if len(simdirs) == 3:
    labels = ["TEP+Storage", "Transmission Only", "Storage Only"]

if len(simdirs) == 2:
    labels = ["TEP+Storage", "Transmission Only"]

with open(f"../../{simdirs[0]}/data.json", 'r') as file:
    data = json.load(file)

with open(f'../../{simdirs[0]}/config.toml', 'r') as toml_file:
    config = toml.load(toml_file)

NUM_HOURS = config["num_hours"]
nonrenewable_types = set(config["nonrenewable_types"])
renewable_types = set(config["renewable_types"])
hours = np.arange(1, NUM_HOURS + 1)

for rep_index in range(1, len(data["param"]["dates"]) + 1):
    date = data["param"]["dates"][rep_index - 1]

    # Calculate renewable production
    renewable_production = np.zeros(NUM_HOURS)
    for key in data["gen"]:
        if data["gen"][key]["gen_type"] in renewable_types:
            renewable_production += data["gen"][key]["profile"][f"{rep_index}"]

    ues = []
    curtailments = []
    curtailment_percentages = []

    for idx, simdir in enumerate(simdirs):
        df_energy = pd.read_csv(f"../../{simdir}/output/{date}/energy.csv")

        # Calculate unserved energy
        ue = np.zeros(NUM_HOURS)
        df_filtered = df_energy[df_energy['Energy_Imbalance'] <= 0]
        df_grouped = df_filtered.groupby('Hour')['Energy_Imbalance'].sum().reset_index()
        for _, row in df_grouped.iterrows():
            hour = int(row['Hour'])
            ue[hour - 1] = -1 * row['Energy_Imbalance']
        ues.append(ue)

        # Compute curtailments
        summed_renewable_df = df_energy.groupby('Hour')[list(renewable_types & set(df_energy.columns))].sum()
        renewable_outputs = summed_renewable_df.sum(axis=1).values
        curtailment = renewable_production - renewable_outputs
        curtailments.append(curtailment)

        # Compute curtailment percentages
        curtailment_percentage = np.zeros(NUM_HOURS)
        for h in range(NUM_HOURS):
            if renewable_production[h] > 0:
                curtailment_percentage[h] = (curtailment[h] / renewable_production[h]) * 100
            else:
                curtailment_percentage[h] = 0
        curtailment_percentages.append(curtailment_percentage)

    # Plot curtailed renewable energy (PUH)
    plt.figure(figsize=(10, 5))
    plt.xlabel("Hour")
    plt.xlim(1, 24)
    plt.ylabel("Curtailed Renewable Energy (puh)")
    plt.title(f"Curtailed Renewable Energy for {date}")
    for curtailment, label in zip(curtailments, labels):
        plt.plot(hours, curtailment, label=label)

    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdirs[0]}/output/{date}/curtailment_combined.png')
    plt.close()

    # Plot percentage of curtailed renewable energy
    plt.figure(figsize=(10, 5))
    plt.xlabel("Hour")
    plt.xlim(1, 24)
    plt.ylabel("Curtailment Percentage (%)")
    plt.title(f"Percentage of Renewable Curtailment for {date}")
    for curtailment_percentage, label in zip(curtailment_percentages, labels):
        plt.plot(hours, curtailment_percentage, label=label)

    plt.legend()
    plt.grid(True)
    plt.savefig(f'../../{simdirs[0]}/output/{date}/curtailment_percentage_combined.png')
    plt.close()
