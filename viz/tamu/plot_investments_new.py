import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
df_storage = pd.read_csv("../../{}/output/storage_investments.csv".format(simdir))
df_storage = df_storage[df_storage['Storage_Energy'] > 0]  # Filter rows with Storage_Energy > 0
df_lines = pd.read_csv("../../{}/output/line_investments.csv".format(simdir))

# Add an invisible dummy point with Storage_Energy = 30
df_storage = pd.concat([
    df_storage,
    pd.DataFrame({'Lat': [0], 'Lon': [0], 'Storage_Energy': [30]})  # Add dummy point
])

# Create the map plot
fig = px.scatter_mapbox(
    df_storage, 
    lat='Lat', 
    lon='Lon', 
    size='Storage_Energy',            # Size based on Storage_Energy
    color='Storage_Energy',           # Color based on Storage_Energy
    center=dict(lat=31.5, lon=-99.5), 
    zoom=4,
    color_continuous_scale="Viridis",  # Color scale for Storage_Energy
    labels={'Storage_Energy': 'Energy (100 MWh)'},  # Adjust labels for clarity
    mapbox_style="carto-positron",
    size_max=30                       # Maximum marker size
)

# Ensure consistent color scale
fig.update_layout(
    coloraxis=dict(cmin=0, cmax=30)  # Set consistent color scale range
)

upgraded_lines = df_lines[df_lines['Upgrade_Lvl'] > 0]
regular_lines = df_lines[df_lines['Upgrade_Lvl'] <= 0]


# Function to add lines to the map
def add_lines(df, color, width_func):
    traces = []
    for _, row in df.iterrows():
        traces.append(go.Scattermapbox(
            lat=[row['Lat1'], row['Lat2']],
            lon=[row['Lon1'], row['Lon2']],
            mode='lines',
            showlegend=False,
            line=dict(color=color, width=width_func(row))
        ))
    return traces

# Add regular lines (thin black lines)
fig.add_traces(add_lines(regular_lines, "black", lambda row: 0.2))

# Add upgraded lines (blue lines with varying width based on Upgrade_Lvl)
fig.add_traces(add_lines(upgraded_lines, "red", lambda row: row['Upgrade_Lvl']))

fig.write_image(f"../../{simdir}/visual/upgrades_new.png")