import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

simdir = sys.argv[1]
# Read CSV file
df_storage = pd.read_csv(f"../../../{simdir}/output/storage_investments.csv")
df_lines = pd.read_csv(f"../../../{simdir}/output/line_investments.csv")

with open(f'../../../{simdir}/data.json', 'r') as file:
    data = json.load(file)

generator_types = set()

for bus in data["bus"]:
    generator_types.update(data["bus"][bus]["gen"].keys())

for gen_type in generator_types:
    df_storage[gen_type + "_capacity"] = 0.0

# Iterate over each row in df_storage
for index, row in df_storage.iterrows():
    bus_index = str(row["Node_Index"])  # Match bus index as a string
    gen_data = data["bus"][bus_index]["gen"]

    # Sum capacity for each generator type
    for gen_type, gen_list in gen_data.items():
        total_capacity = sum(data["gen"][str(gen_id)]["pmax"] for gen_id in gen_list)
        df_storage.at[index, gen_type + "_capacity"] = total_capacity

df_storage["nonrenewable_capacity"] = df_storage[["coal_capacity", "nuclear_capacity", "ng_capacity"]].sum(axis=1)

# Filter data for each generator type
df_solar = df_storage[df_storage["solar_capacity"] > 0]
df_wind = df_storage[df_storage["wind_capacity"] > 0]
df_nonrenewable = df_storage[df_storage["nonrenewable_capacity"] > 0]

# Plot for Solar
fig = px.scatter_mapbox(
    df_solar,
    lat='Lat',
    lon='Lon',
    size='solar_capacity',               # Size based on solar capacity
    color='solar_capacity',              # Color based on solar capacity
    color_continuous_scale='YlOrRd_r',    # Use viridis color scale
    center=dict(lat=31.5, lon=-99.5),    # Center the map
    zoom=4,                              # Set zoom level
    labels={'solar_capacity': 'Solar Capacity (100 MW)'},
    mapbox_style="carto-positron",       # Map style
)

# Function to add lines to the map
def add_lines(df, color="black", width=0.2):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattermapbox(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width)
        ))
    return traces

# Add all lines (black, width 0.2)
fig.add_traces(add_lines(df_lines, color="black", width=0.2))

# Save the map
fig.write_image(f"../../../{simdir}/visual/solar_cap.png")