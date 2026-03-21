using CSV
using DataFrames

LINE_MAX = 3.0
STORAGE_MAX = 30.0

function count_max_upgrades(simdir)
    line_path = "$(simdir)/output/line_investments.csv"
    storage_path = "$(simdir)/output/storage_investments.csv"
    df_lines = CSV.read(line_path, DataFrame)
    df_storage = CSV.read(storage_path, DataFrame)

    lines_max = count(row -> row[:Upgrade_Lvl] == LINE_MAX, eachrow(df_lines))
    storage_max = count(row -> row[:Storage_Energy] == STORAGE_MAX, eachrow(df_storage))

    return lines_max, storage_max
end