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

fig = go.Figure()
# Add nonrenewable plants (e.g., gray circles)
fig.add_trace(go.Scattergeo(
    lat=df_nonrenewable['Lat'],
    lon=df_nonrenewable['Lon'],
    mode='markers',
    marker=dict(
        size=18,           # Static size
        color='gray',      # Color
        opacity=0.8,       # Opacity
        line=dict(width=1, color='black')  # Optional: border
    ),
    name='Nonrenewable',
    showlegend=True
))

# Add wind plants (e.g., blue circles)
fig.add_trace(go.Scattergeo(
    lat=df_wind['Lat'],
    lon=df_wind['Lon'],
    mode='markers',
    marker=dict(
        size=16,           # Different size
        color='blue',      # Color
        opacity=0.7,       # Opacity
        line=dict(width=1, color='darkblue')
    ),
    name='Wind',
    showlegend=True
))

# Add solar plants (e.g., yellow circles)
fig.add_trace(go.Scattergeo(
    lat=df_solar['Lat'],
    lon=df_solar['Lon'],
    mode='markers',
    marker=dict(
        size=14,            # Different size
        color='red',    # Color
        opacity=0.9,       # Opacity
        line=dict(width=1, color='orange')
    ),
    name='Solar',
    showlegend=True
))

# Updated parameter names
fig.update_geos(
    lonaxis_range=[-105, -94],
    lataxis_range=[25.5, 36], 
    showland=True,
    showocean=True,
    oceancolor="lightblue",
    showlakes=True,
    lakecolor="lightblue",
    showcountries=True,
    countrycolor="lightgray",
    countrywidth=0,  # Width of country borders
    showsubunits=True,  # This shows states/provinces
    subunitcolor="lightgray",
    subunitwidth=3
)

# Update layout with colorbar customization
fig.update_layout(
    legend=dict(
        x=0.78,                    # Horizontal position (0=left, 1=right)
        font=dict(
            size=20,               # Font size for legend text
            color="black",         # Font color
        ),
        orientation="v",           # "v" for vertical, "h" for horizontal
        xanchor="left",           # Anchor point for x position
        yanchor="top"             # Anchor point for y position
    )
)

# Function to add lines - use Scattergeo for lines too
def add_lines(df, color, width_func):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattergeo(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width_func(row))
        ))
    return traces

# Add all lines (black, width 0.2)
fig.add_traces(add_lines(df_lines, "black", lambda row: 0.4))

# Save the map
fig.write_html(f"../../../{simdir}/cap.html")