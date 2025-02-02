import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
# Read data json
with open(f"../../{simdir}/data.json", 'r') as file:
    data = json.load(file)

renewable_types = data["param"]["renewable_types"]
nonrenewable_types = data["param"]["nonrenewable_types"]

for date in data["param"]["dates"]:
    # Read CSV file
    df_energy = pd.read_csv(f"../../{simdir}/output/{date}/energy.csv")

    # Check if required columns exist
    required_cols = {'Lon', 'Lat', 'Node_Name', 'Hour', 'Charge', 'Discharge'}
    if not required_cols.issubset(df_energy.columns):
        print(f"Missing columns in {date}/energy.csv. Skipping this file.")
        continue

    # Ensure 'Hour' is treated as an integer
    df_energy['Hour'] = pd.to_numeric(df_energy['Hour'], errors='coerce')

    # Reshape the DataFrame to include both Charge and Discharge in a single plot
    df_charge = df_energy[['Lon', 'Lat', 'Node_Name', 'Hour', 'Charge']].copy()
    df_charge['Type'] = 'Charge'
    df_charge['Value'] = df_charge['Charge']

    df_discharge = df_energy[['Lon', 'Lat', 'Node_Name', 'Hour', 'Discharge']].copy()
    df_discharge['Type'] = 'Discharge'
    df_discharge['Value'] = df_discharge['Discharge']

    # Combine the two DataFrames
    df_combined = pd.concat([df_charge, df_discharge], ignore_index=True)

    # Create the scatter plot with animation
    fig = px.scatter_geo(
        df_combined,
        lon='Lon',
        lat='Lat',
        hover_name='Node_Name',
        size='Value',
        color='Type',
        color_discrete_map={'Charge': 'red', 'Discharge': 'blue'},
        animation_frame='Hour',
        scope='usa',
        title='Charge and Discharge Animation'
    )

    # Save the plot as an HTML file
    fig.write_html(f"../../{simdir}/output/{date}/charge_discharge.html")


for date in data["param"]["dates"]:
    # Read CSV file
    df_energy = pd.read_csv(f"../../{simdir}/output/{date}/energy.csv")

    # Calculate total renewable and nonrenewable production
    df_energy['Renewable'] = df_energy[renewable_types].sum(axis=1)
    df_energy['Nonrenewable'] = df_energy[nonrenewable_types].sum(axis=1)

    # Ensure 'Hour' is treated as an integer
    df_energy['Hour'] = pd.to_numeric(df_energy['Hour'], errors='coerce')

    # Reshape the DataFrame to include both Renewable and Nonrenewable in a single plot
    df_renewable = df_energy[['Lon', 'Lat', 'Node_Name', 'Hour', 'Renewable']].copy()
    df_renewable['Type'] = 'Renewable'
    df_renewable['Value'] = df_renewable['Renewable']

    df_nonrenewable = df_energy[['Lon', 'Lat', 'Node_Name', 'Hour', 'Nonrenewable']].copy()
    df_nonrenewable['Type'] = 'Nonrenewable'
    df_nonrenewable['Value'] = df_nonrenewable['Nonrenewable']

    # Combine the two DataFrames
    df_combined = pd.concat([df_renewable, df_nonrenewable], ignore_index=True)

    # Sort by 'Hour' to ensure proper order
    df_combined = df_combined.sort_values(by='Hour')

    # Create the scatter plot with animation
    fig = px.scatter_geo(
        df_combined,
        lon='Lon',
        lat='Lat',
        hover_name='Node_Name',
        size='Value',
        color='Type',
        color_discrete_map={'Renewable': 'green', 'Nonrenewable': 'brown'},
        animation_frame='Hour',
        scope='usa',
        title='Renewable vs. Nonrenewable Energy Production Animation'
    )

    # Save the plot as an HTML file
    fig.write_html(f"../../{simdir}/output/{date}/renewable_nonrenewable.html")

    
"""
for (i, date) in enumerate(data["param"]["dates"]):
    records = []
    for node_id, node_data in data["bus"].items():
        lat = node_data["lat"]
        lon = node_data["lon"]
        # Loop through each hour (24 hours) of the load data
        for hour, load in enumerate(node_data["load"][f'{i+1}']):
            records.append({
                "node": node_id,
                "lat": lat,
                "lon": lon,
                "load": load,
                "hour": hour
            })
    
    df = pd.DataFrame(records)

    # Create the scatter plot with animation
    fig = px.scatter_geo(
        df,
        lon='lon',
        lat='lat',
        hover_name='node',
        size='load',
        animation_frame='hour',
        scope='usa',
        title='Load Animation'
    )

    # Save the plot as an HTML file
    fig.write_html(f"../../{simdir}/output/{date}/load.html")"""

