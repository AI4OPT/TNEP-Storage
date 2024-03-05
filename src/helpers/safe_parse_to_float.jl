function safe_parse_to_float(x)
    if ismissing(x)
        return 0.0
    else
        x_no_commas = replace(x, "," => "")
        return tryparse(Float64, x_no_commas) === nothing ? 0.0 : parse(Float64, x_no_commas)
    end
end