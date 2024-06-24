import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_energy = pd.read_csv(f"../../{simdir}/output/energy.csv")
production_cols = [col for col in df_energy.columns if col.endswith('_production')]
df_energy['Total_Curtailment'] = 0
# Calculate the row-wise sum of differences
for prod_col in production_cols:
    original_col = prod_col.replace('_production', '')
    df_energy['Total_Curtailment'] += df_energy[prod_col] - df_energy[original_col]

fig = px.scatter_geo(df_energy, 
                     lon='Lon',
                     lat='Lat',
                     hover_name='Node_Name', 
                     size='Total_Curtailment',
                     color_discrete_sequence=['blue'],  # Use a single color for clarity
                     animation_frame='Hour',
                     scope='usa')

# Save the plot as an HTML file
fig.write_html(f"../../{simdir}/visual/curtailed_energy.html")



