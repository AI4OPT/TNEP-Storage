import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import json
from collections import Counter

TAMU_NETWORK = "/storage/home/hcoda1/1/kwu381/TNEP-Storage/data/topology/tamu/power_system_data.json"

with open(TAMU_NETWORK, 'r') as file:
    data = json.load(file)


iso_list = [bus["iso"] for bus in data["bus"].values()]
iso_counts = Counter(iso_list)
sorted_iso_counts = iso_counts.most_common()
iso_names, bus_counts = zip(*sorted_iso_counts)

# Create a bar plot
plt.figure(figsize=(10, 6))
plt.bar(iso_names, bus_counts, color='skyblue')
plt.xlabel('ISO')
plt.ylabel('Number of Buses')
plt.title('Number of Buses in Each ISO')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
plt.savefig("/storage/home/hcoda1/1/kwu381/TNEP-Storage/data/topology/tamu/iso_size_counts.png")