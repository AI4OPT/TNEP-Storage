import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
df_storage = pd.read_csv("../../{}/storage_investments.csv".format(simdir))
df_storage = df_storage[df_storage['Storage_Energy'] > 0]  # Filter rows with Storage_Energy > 0
df_lines = pd.read_csv("../../{}/line_investments.csv".format(simdir))


# Add an invisible dummy point with Storage_Energy = 30
df_storage = pd.concat([
    df_storage,
    pd.DataFrame({'Lat': [0], 'Lon': [0], 'Storage_Energy': [12]})  # Add dummy point
])

# Create the map plot
"""
fig = px.scatter_geo(
    df_storage, 
    lat='Lat', 
    lon='Lon', 
    size='Storage_Energy',            # Size based on Storage_Energy
    color='Storage_Energy',           # Color based on Storage_Energy
    color_continuous_scale="Viridis",  # Color scale for Storage_Energy
    labels={'Storage_Energy': 'Energy (250 MWh)'},  # Adjust labels for clarity
    size_max=30,                       # Maximum marker size
    scope="usa"
)

# Updated parameter names
fig.update_geos(
center=dict(lat=31.1255, lon=-98.8142),
projection_scale=2.6,
showland=True,
landcolor="blanchedalmond",
showocean=True,
oceancolor="lightblue",
showlakes=True,
lakecolor="lightblue",
showcountries=True,
countrycolor="lightgray",
countrywidth=0, # Width of country borders
showsubunits=True, # This shows states/provinces
subunitcolor="lightgray",
subunitwidth=3
)"""

# Create the map plot
fig = px.scatter_geo(
    df_storage, 
    lat='Lat', 
    lon='Lon', 
    size='Storage_Energy',            # Size based on Storage_Energy
    color='Storage_Energy',           # Color based on Storage_Energy
    color_continuous_scale="Viridis",  # Color scale for Storage_Energy
    labels={'Storage_Energy': 'Energy (250 MWh)'},  # Adjust labels for clarity
    size_max=30,                       # Maximum marker size
)

# Updated parameter names
fig.update_geos(
    lonaxis_range=[-105, -94],
    lataxis_range=[25.5, 36], 
    showland=True,
    # landcolor="blanchedalmond",
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

# Ensure consistent color scale
# Update layout with colorbar customization
fig.update_layout(
    coloraxis=dict(
        cmin=0, 
        cmax=12,  # Set consistent color scale range
        colorbar=dict(
            x=0.78,           # Move colorbar closer to plot (0-1, where 1 is far right)
            len=0.95,          # Length of colorbar (0-1)
            thickness=30,     # Thickness of colorbar in pixels
            title_font_size=24,  # Title font size
            tickfont_size=20     # Tick label font size
        )
    )
)

upgraded_lines = df_lines[df_lines['Upgrade_Lvl'] > 0]
regular_lines = df_lines[df_lines['Upgrade_Lvl'] <= 0]


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

# Add regular lines (thin black lines)
fig.add_traces(add_lines(regular_lines, "black", lambda row: 0.4))

# Add upgraded lines (blue lines with varying width based on Upgrade_Lvl)
# fig.add_traces(add_lines(upgraded_lines, "red", lambda row: row['Upgrade_Lvl']))
fig.add_traces(add_lines(upgraded_lines, "red", lambda row: row['Upgrade_Lvl']/1.25))

fig.write_html(f"../../{simdir}/upgrades_new.html")