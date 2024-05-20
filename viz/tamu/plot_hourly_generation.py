import matplotlib.pyplot as plt
import pandas as pd
import sys
import json
import numpy as np
import toml

def average_profile(profile_dict, num_hours):
    num_scenarios = len(profile_dict)
    scenario_sum = np.zeros(num_hours)
    # Sum up all the arrays across scenarios
    for scenario in profile_dict.values():
        scenario_sum += np.array(scenario)
    # Calculate the average of scenarios
    return scenario_sum / num_scenarios

simdir = sys.argv[1]
# Read CSV file
# df_storage = pd.read_csv("../../{}/output/storage_investments.csv".format(simdir))
# df_lines = pd.read_csv("../../{}/output/line_investments.csv".format(simdir))
df_energy = pd.read_csv("../../{}/output/energy.csv".format(simdir))

with open(f'../../{simdir}/config.toml', 'r') as toml_file:
    config = toml.load(toml_file)

with open(f'../../{simdir}/data.json', 'r') as file:
    data = json.load(file)

NUM_HOURS = config["num_hours"]

# Compute the total nonrenewable capacity
nonrenewable_types = set(config["nonrenewable_types"])
renewable_types = set(config["renewable_types"])
total_pmax_nonrenewable = 0
for key in data["gen"]:
    # Check if the 'gen_type' of this generator is in the set of nonrenewable types
    if data["gen"][key]["gen_type"] in nonrenewable_types:
        # Add the 'pmax' value of this generator to the total
        total_pmax_nonrenewable += data["gen"][key].get("pmax", 0)

# Initialize a NumPy array to store the sum of average loads for each hour (24 hours)
total_load = np.zeros(NUM_HOURS)

# Calculate the average load per hour for each bus and sum them
for key in data["bus"]:
    total_load += average_profile(data["bus"][key]["load"], NUM_HOURS)

renewable_production = np.zeros(NUM_HOURS)
for key in data["gen"]:
    if data["gen"][key]["gen_type"] in renewable_types:
        renewable_production += average_profile(data["gen"][key]["profile"], NUM_HOURS)

# Now compute the actual renewable output
summed_renewable_df = df_energy.groupby('Hour')[list(renewable_types & set(df_energy.columns))].sum()
renewable_outputs = summed_renewable_df.sum(axis=1).values
# Now compute the actual nonrenewable output
summed_nonrenewable_df = df_energy.groupby('Hour')[list(nonrenewable_types & set(df_energy.columns))].sum()
nonrenewable_outputs = summed_nonrenewable_df.sum(axis=1).values
# Now compute charge and discharge
summed_charge_df = df_energy.groupby('Hour')['Charge'].sum()
charge_amount = summed_charge_df.values
summed_discharge_df = df_energy.groupby('Hour')['Discharge'].sum()
discharge_amount = summed_discharge_df.values

hours = np.arange(1, NUM_HOURS+1)
plt.figure(figsize=(10, 5))  # Set the figure size
plt.plot(hours, total_load, label='Total Load', marker='o')  # Plot total load
plt.plot(hours, renewable_production, label='Renewable Production', marker='o')  # Plot renewable production
plt.plot(hours, renewable_outputs, label='Renewable Outputs', marker='o')  # Plot renewable outputs
plt.plot(hours, np.array([total_pmax_nonrenewable] * NUM_HOURS), label='Nonrenewable Capacity', marker='o')  # Plot nonrenewable capacity
plt.plot(hours, nonrenewable_outputs, label='Nonrenewable Outputs', marker='o')  # Plot renewable outputs
plt.plot(hours, charge_amount, label='Charging Amount', marker='o')
plt.plot(hours, discharge_amount, label='Discharge Amount', marker='o')

plt.title('Energy Production and Load Over Time')
plt.xlabel('Hour of the Day')
plt.ylabel('Energy (MW)')
plt.legend()
plt.grid(True)
plt.savefig(f'../../{simdir}/visual/hourly_generation.png')

# Second plot: stacked plot
stacked_production = renewable_production + nonrenewable_outputs + discharge_amount
stacked_no_discharge = renewable_production + nonrenewable_outputs
plt.figure(figsize=(10, 5))  # Set the figure size
plt.plot(hours, total_load, label='Total Load', marker='o', linestyle='-')  # Plot total load
plt.plot(hours, stacked_production, label='Stacked Production', marker='o', linestyle='-')  # Plot stacked production
plt.plot(hours, stacked_no_discharge, label='Stacked No Discharge', marker='o', linestyle='-')  # Plot stacked production
plt.title('Stacked Energy Production vs Total Load')
plt.xlabel('Hour of the Day')
plt.ylabel('Energy (MW)')
plt.legend()
plt.grid(True)
plt.savefig(f'../../{simdir}/visual/stacked_hourly_generation.png')



