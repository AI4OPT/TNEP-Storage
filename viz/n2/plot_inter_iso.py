import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import numpy as np

INTER_ISO_PATH = "../../data/topology/tamu/n2/inter_iso_branches.csv"
df = pd.read_csv(INTER_ISO_PATH)

unique_inter_iso = df['inter_iso'].unique()
# Map each unique value to an integer
mapping = {iso: i for i, iso in enumerate(unique_inter_iso)}
# Normalize these integers to [0, 1] for continuous color scale
normalized_values = np.array([mapping[iso] for iso in df['inter_iso']]) / max(mapping.values())

# Use the viridis continuous color scale from plotly
color_scale = px.colors.sequential.Viridis
# Map normalized values to actual colors
colors = np.interp(normalized_values, [0, 1], [0, len(color_scale)-1]).astype(int)
color_values = [color_scale[i] for i in colors]

fig = go.Figure()

# Add lines with colors depending on 'inter_iso'
for (index, row), color in zip(df.iterrows(), color_values):
    fig.add_trace(go.Scattergeo(
        lon = [row['lon1'], row['lon2']],
        lat = [row['lat1'], row['lat2']],
        mode = 'lines',
        line=dict(width=1, color=color),
    ))

fig.update_layout(
    title='Nodes',
    geo_scope='usa',
    showlegend=False
)

fig.write_image("../../data/topology/tamu/n2/inter_iso.png")
