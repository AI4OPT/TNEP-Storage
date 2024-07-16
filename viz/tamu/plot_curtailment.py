import plotly.express as px
import pandas as pd
import sys

simdir = sys.argv[1]
# Read CSV file
df_energy = pd.read_csv(f"../../{simdir}/output/energy.csv")
production_cols = [col for col in df_energy.columns if col.endswith('_production')]

# Calculate curtailment for each production type
for prod_col in production_cols:
    original_col = prod_col.replace('_production', '')
    curtailment_col = f'Curtailment_{original_col}'
    df_energy[curtailment_col] = df_energy[prod_col] - df_energy[original_col]

# Melt the DataFrame to long format for plotting
df_melted = df_energy.melt(id_vars=['Lon', 'Lat', 'Node_Name', 'Hour'], 
                           value_vars=[f'Curtailment_{col.replace("_production", "")}' for col in production_cols],
                           var_name='Curtailment_Type', 
                           value_name='Curtailment')

# Create a mapping for colors
color_map = {
    'Curtailment_wind': 'green',
    'Curtailment_solar': 'yellow',
    'Curtailment_wind_offshore': 'green',
    'Curtailment_hydro': 'blue'
    # Add more colors if there are more types
}

# Plot with different colors for each curtailment type
fig = px.scatter_geo(df_melted, 
                     lon='Lon',
                     lat='Lat',
                     hover_name='Node_Name', 
                     size='Curtailment',
                     color='Curtailment_Type',
                     color_discrete_map=color_map,
                     animation_frame='Hour',
                     scope='usa')

# Save the plot as an HTML file
fig.write_html(f"../../{simdir}/visual/curtailed_energy.html")



