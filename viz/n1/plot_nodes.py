import plotly.graph_objects as go
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_storage = pd.read_csv(f"../../{simdir}/output/storage_investments.csv")
df_lines = pd.read_csv(f"../../{simdir}/output/line_investments.csv")

fig = go.Figure()

# Add lines
for _, row in df_lines.iterrows():
    fig.add_trace(go.Scattergeo(
        lon = [row['Lon1'], row['Lon2']],
        lat = [row['Lat1'], row['Lat2']],
        mode = 'lines',
        line=dict(width=1, color='blue'),
    ))

fig.add_trace(go.Scattergeo(
    lon = df_storage['Lon'],
    lat = df_storage['Lat'],
    mode = 'markers',
    marker_color = 'red',
    ))

fig.update_layout(
    title='Texas Transmission and Storage Investments',
    geo=dict(
        scope='usa',  # Still uses USA map for better rendering
        resolution=50,  # High-resolution map
        showland=True,
        landcolor='lightgray',
        subunitcolor='white',
        countrycolor='white',
        center=dict(lat=31.9686, lon=-99.9018),  # Center on Texas
        projection=dict(type='mercator'),  # Mercator projection for clarity
        lataxis=dict(range=[25.8, 36.5]),  # Latitude range for Texas
        lonaxis=dict(range=[-106.6, -93.5])  # Longitude range for Texas
    )
)
fig.write_image(f"../../{simdir}/visual/nodes.png")