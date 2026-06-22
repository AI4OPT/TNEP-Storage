import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import json

# EXAMPLE USAGE: 
# cd viz/tamu/topology
# python plot_capacities.py examples/example_simdir

# OUTPUT: an html file of the grid

simdir = sys.argv[1]

with open(f'../../../{simdir}/data.json', 'r') as file:
    data = json.load(file)

generator_types = set()

for bus in data["bus"]:
    generator_types.update(data["bus"][bus]["gen"].keys())

# Build a row per bus with lat/lon and generator capacity by type
bus_records = []
for bus_index, bus in data["bus"].items():
    capacities = {gen_type + "_capacity": 0.0 for gen_type in generator_types}
    for gen_type, gen_list in bus["gen"].items():
        capacities[gen_type + "_capacity"] = sum(data["gen"][str(gen_id)]["pmax"] for gen_id in gen_list)
    bus_records.append({"Lat": bus["lat"], "Lon": bus["lon"], **capacities})

df_storage = pd.DataFrame(bus_records)
df_storage["nonrenewable_capacity"] = df_storage[["coal_capacity", "nuclear_capacity", "ng_capacity"]].sum(axis=1)

# Build a row per branch with the lat/lon of its two endpoint buses
line_records = []
for branch in data["branch"].values():
    f_bus = data["bus"][str(branch["f_bus"])]
    t_bus = data["bus"][str(branch["t_bus"])]
    line_records.append({
        "Lat1": f_bus["lat"], "Lon1": f_bus["lon"],
        "Lat2": t_bus["lat"], "Lon2": t_bus["lon"],
    })

df_lines = pd.DataFrame(line_records)

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