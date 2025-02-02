import matplotlib.pyplot as plt
import pandas as pd
import sys
import json
import numpy as np
import toml

# Parse multiple simdir arguments
simdirs = sys.argv[1:]

# Check if the length of simdirs is 4 and assign hardcoded labels
if len(simdirs) == 4:
    labels = ["no_upgrades", "only_transmission", "only_storage", "all_upgrades"]
else:
    raise ValueError("Expected exactly 4 simdirs for hardcoded labels.")

# Dictionary to store padded unserved energy values for each simdir
all_padded_unserved_energy = {}

# Process each simdir
for idx, simdir in enumerate(simdirs):
    # Read data json
    with open(f"../../{simdir}/data.json", 'r') as file:
        data = json.load(file)

    with open(f'../../{simdir}/config.toml', 'r') as toml_file:
        config = toml.load(toml_file)

    NUM_HOURS = config["num_hours"]

    # Initialize a list to store the unserved energy values for histogram input
    unserved_energy_values = []

    for (i, each_date) in enumerate(data["param"]["dates"]):
        prob = data["param"]["representative_prob"][f'{i+1}']
        df_energy = pd.read_csv(f'../../{simdir}/output/{each_date}/energy.csv')

        # Filter and group data
        df_filtered = df_energy[df_energy['Energy_Imbalance'] <= 0]
        df_grouped = df_filtered.groupby('Hour')['Energy_Imbalance'].sum().reset_index()
        
        # Collect unserved energy values weighted by probability
        for _, row in df_grouped.iterrows():
            unserved_energy = -1 * row['Energy_Imbalance']
            unserved_energy_values.extend([unserved_energy] * int(prob * 365))

    # Sort and pad the unserved energy array
    sorted_unserved_energy = np.sort(unserved_energy_values)[::-1]
    total_hours = 8760
    if len(sorted_unserved_energy) < total_hours:
        padded_unserved_energy = np.pad(sorted_unserved_energy, (0, total_hours - len(sorted_unserved_energy)), 'constant')
    else:
        padded_unserved_energy = sorted_unserved_energy[:total_hours]
    
    # Store padded unserved energy for this simdir
    all_padded_unserved_energy[labels[idx]] = padded_unserved_energy

# Determine the maximum y-axis value across all simulations
max_y_value = max(max(values) for values in all_padded_unserved_energy.values())

# Create a single plot with all simdirs
plt.figure(figsize=(10, 5))
for label, padded_unserved_energy in all_padded_unserved_energy.items():
    plt.plot(padded_unserved_energy, marker='o', markersize=1, linestyle='-', label=label)

# Set the title, labels, and y-axis limit
plt.title('Unserved Energy Sorted by Hour (Descending)')
plt.xlabel('Hour (sorted by descending unserved energy)')
plt.ylabel('Unserved Energy (puh)')
plt.ylim(0, max_y_value)  # Set uniform y-axis limit

# Set x-axis limit to focus on the region with significant values
plt.xlim(0, 2000)  # Adjust this value based on where the data becomes mostly zero

plt.grid(True)
plt.legend(title="Simulations")  # Add legend with hardcoded labels

# Save the plot with all simdirs on the same plot
plt.savefig(f'../../{simdirs[0]}/visual/unserved_energy_sorted_combined.png')
plt.close()
