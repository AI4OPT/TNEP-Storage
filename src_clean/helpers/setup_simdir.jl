function setup_simdir(simdir)
    if !isdir(joinpath(simdir, "output"))
        mkdir(joinpath(simdir, "output"))
    end

    if !isdir(joinpath(simdir, "visual"))
        mkdir(joinpath(simdir, "visual"))
    end
end