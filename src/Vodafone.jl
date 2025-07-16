module Vodafone

using GLMakie

const DATA_PATH = joinpath(@__DIR__, "..", "data")
export DATA_PATH

# Coordinate transformations

using CoordRefSystems
using GeoInterface
using GADM

const BNG = CoordRefSystems.get(EPSG{27700})

function latlon2tb(lat, lon; bin_size = 25)
    (x, y) = latlon2xy(lat, lon)
    x = x - (x % bin_size)
    y = y - (y % bin_size)
    return "$x$y"
end

function latlon2xy(lat, lon)
    p = convert(BNG, LatLon(lat, lon))
    x = floor(Int, p.x.val)
    y = floor(Int, p.y.val)
    return (x, y)
end

# Data processing

using CSV
using DataFrames
using Arrow

# Outline of GBR

function get_outline(gbr = GADM.get("GBR"))
    coords = GeoInterface.coordinates(gbr)
    # Convert coordinates to BNG
    return [[[[latlon2xy(pt[2], pt[1]) for pt in pts] for pts in ptss] for ptss in ptsss]
            for ptsss in coords]
end

function plot_outline!(ax, coords)
    for coord in coords
        for pts in coord
            for pt in pts
                lines!(ax, [xx[1] for xx in pt], [xx[2] for xx in pt],
                    color = :red, linewidth = 0.5)
            end
        end
    end
end

# Plotting OFCOM data from 2024

function read_ofcom_2024_raw()
    columns = ["latitude", "longitude", "pci_top1_vf", "rsrp_top1_vf", "earfcn_top1_vf",
        "pci_top2_vf", "rsrp_top2_vf", "earfcn_top2_vf", "pci_top3_vf", "rsrp_top3_vf",
        "earfcn_top3_vf", "pci_top4_vf", "rsrp_top4_vf", "earfcn_top4_vf"]
    data = CSV.read(
        joinpath(DATA_PATH, "Tony_Vernon_2024 Ofcom 4G Measurement Data",
            "4g-lte-2024-mobile-signal-measurement-data.csv"),
        DataFrame;
        select = columns)
    pts = latlon2xy.(data.latitude, data.longitude)
    data.x = [x[1] for x in pts]
    data.y = [x[2] for x in pts]
    return data
end

function read_ofcom_2024_arrow()
    return DataFrame(
        Arrow.Table(joinpath(
            DATA_PATH, "4g-lte-2024-mobile-signal-measurement-data.arrow"));
        copycols = false)
end

function get_heatmap(data; bin_size = 25)
    x_range = (minimum(data.x), maximum(data.x))
    y_range = (minimum(data.y), maximum(data.y))
    rsrp = Matrix{Union{Float64, Missing}}(
        missing, ((x_range[2] - x_range[1]) ÷ bin_size + 1),
        (y_range[2] - y_range[1]) ÷ bin_size + 1)
    fill!(rsrp, missing)
    for row in eachrow(data)
        rsrp_top1_vf = ((row.earfcn_top1_vf !== missing) &&
                        (row.rsrp_top1_vf !== missing) && (row.earfcn_top1_vf == 6300)) ?
                       row.rsrp_top1_vf : missing
        rsrp_top2_vf = ((row.earfcn_top2_vf !== missing) &&
                        (row.rsrp_top2_vf !== missing) && (row.earfcn_top2_vf == 6300)) ?
                       row.rsrp_top2_vf : missing
        rsrp_top3_vf = ((row.earfcn_top3_vf !== missing) &&
                        (row.rsrp_top3_vf !== missing) && (row.earfcn_top3_vf == 6300)) ?
                       row.rsrp_top3_vf : missing
        rsrp_top4_vf = ((row.earfcn_top4_vf !== missing) &&
                        (row.rsrp_top4_vf !== missing) && (row.earfcn_top4_vf == 6300)) ?
                       row.rsrp_top4_vf : missing
        rsrp_all = collect(skipmissing([
            rsrp_top1_vf, rsrp_top2_vf, rsrp_top3_vf, rsrp_top4_vf]))
        if isempty(rsrp_all)
            continue
        end
        rsrp_max = maximum(rsrp_all)
        x_bin = (row.x - x_range[1]) ÷ bin_size + 1
        y_bin = (row.y - y_range[1]) ÷ bin_size + 1
        rsrp[x_bin, y_bin] = rsrp_max
    end
    x = x_range[1]:bin_size:x_range[2]
    y = y_range[1]:bin_size:y_range[2]
    return (x, y, rsrp)
end

# Note: can use `construct_plot(read_ofcom_2024_raw(), GADM.get("GBR"))` to avoid having to create the Arrow file.
function construct_plot(data = read_ofcom_2024_arrow(), gbr = GADM.get("GBR"))
    coords = get_outline(gbr)
    (x, y, rsrp) = get_heatmap(data; bin_size = 1000)
    fig = Figure(size = (1200, 800))
    ax = Axis(fig[1, 1], title = "Vodafone 4G RSRP Heatmap",
        xlabel = "Easting (m)", ylabel = "Northing (m)", aspect = 1)
    plot_outline!(ax, coords)
    heatmap!(ax, x, y, rsrp)
    return fig
end

# Extract difference between OFCOM data and Vodafone predictions

function combine_data_2024()
    local pred_rsrp, pred_pci, x_origin, y_origin, bin_size
    open(joinpath(DATA_PATH, "Tony_Vernon_2024 LTE MIMO Terminal [No Loss] (DL) 800MHz",
        "LTE  MIMO Terminal [No Loss] (DL) 800MHz.asc")) do f
        ncols = parse(Int, split(readline(f))[2])
        nrows = parse(Int, split(readline(f))[2])
        x_origin = parse(Float64, split(readline(f))[2])
        y_origin = parse(Float64, split(readline(f))[2])
        bin_size = parse(Float64, split(readline(f))[2])
        nodata_value = parse(Float64, split(readline(f))[2])
        @assert ncols == 6860
        @assert nrows == 12380
        pred_rsrp = Matrix{Union{Float64, Missing}}(undef, nrows, ncols)
        for i in 1:nrows
            line = readline(f)
            values = split(line)
            for j in 1:ncols
                value = parse(Float64, values[j])
                pred_rsrp[i, j] = value == nodata_value ? missing : value
            end
        end
    end
    open(joinpath(DATA_PATH, "Tony_Vernon_2024 LTE MIMO Terminal [No Loss] (DL) 800MHz",
        "LTE  MIMO Terminal [No Loss] (DL) 800MHz.svr.asc")) do f
        ncols = parse(Int, split(readline(f))[2])
        nrows = parse(Int, split(readline(f))[2])
        x_origin = parse(Float64, split(readline(f))[2])
        y_origin = parse(Float64, split(readline(f))[2])
        bin_size = parse(Float64, split(readline(f))[2])
        nodata_value = parse(Int, split(readline(f))[2])
        @assert ncols == 6860
        @assert nrows == 12380
        pred_pci = Matrix{Union{Int, Missing}}(undef, nrows, ncols)
        for i in 1:nrows
            line = readline(f)
            values = split(line)
            for j in 1:ncols
                value = parse(Int, values[j])
                pred_pci[i, j] = value == nodata_value ? missing : value
            end
        end
    end
    index2cell_data = CSV.read(
        joinpath(DATA_PATH, "Tony_Vernon_2024 LTE MIMO Terminal [No Loss] (DL) 800MHz",
            "LTE  MIMO Terminal [No Loss] (DL) 800MHz.svr.csv"), DataFrame; header = [
            "index", "name"])
    index2cell = Dict{Int, String}()
    for row in eachrow(index2cell_data)
        index2cell[row.index] = row.name
    end
    cell2pci_data = CSV.read(
        joinpath(DATA_PATH, "L8_Sept24 Sites TXs and Cells", "L8_Sept24 Cells.txt"),
        DataFrame; select = ["Transmitter", "Physical Cell ID"])
    cell2pci = Dict{String, Int}()
    for row in eachrow(cell2pci_data)
        if row.var"Physical Cell ID" !== missing
            cell2pci[row.Transmitter] = row.var"Physical Cell ID"
        else
            cell2pci[row.Transmitter] = -1
        end
    end
    # Convert the predicted PCI values to the correct format
    for i in eachindex(pred_pci)
        if pred_pci[i] !== missing
            pred_pci[i] = cell2pci[index2cell[pred_pci[i]]]
        end
    end
    # Flip the matrices to match the expected coordinate system
    pred_rsrp = reverse(pred_rsrp, dims = 1)
    pred_pci = reverse(pred_pci, dims = 1)
    # Read the OFCOM data and combine it with the predictions
    columns = Dict(["latitude" => Float64, "longitude" => Float64, "pci_top1_vf" => Int,
        "rsrp_top1_vf" => Float64, "earfcn_top1_vf" => Int, "pci_top2_vf" => Int,
        "rsrp_top2_vf" => Float64, "earfcn_top2_vf" => Int, "pci_top3_vf" => Int,
        "rsrp_top3_vf" => Float64, "earfcn_top3_vf" => Int, "pci_top4_vf" => Int,
        "rsrp_top4_vf" => Float64, "earfcn_top4_vf" => Int, "month_year" => String,
        "hour_ref" => String])
    data = CSV.Rows(
        joinpath(DATA_PATH, "Tony_Vernon_2024 Ofcom 4G Measurement Data",
            "4g-lte-2024-mobile-signal-measurement-data.csv");
        select = collect(keys(columns)),
        types = columns)
    open(joinpath(DATA_PATH, "4g-lte-2024-mobile-signal-measurement-data-combined.csv"),
        "w") do f
        println(f,
            "x,y,month_year,hour_ref,voda_pci,voda_rsrp,ofcom_pci_top1,ofcom_rsrp_top1,ofcom_pci_top2,ofcom_rsrp_top2,ofcom_pci_top3,ofcom_rsrp_top3,ofcom_pci_top4,ofcom_rsrp_top4")
        for row in data
            pt = latlon2xy(row.latitude, row.longitude)
            i_x = round(Int, (pt[1] - x_origin) ÷ bin_size + 1)
            i_y = round(Int, (pt[2] - y_origin) ÷ bin_size + 1)
            if i_x < 1 || i_y < 1 || i_x > size(pred_pci, 2) || i_y > size(pred_pci, 1)
                throw(ArgumentError("Point $pt is out of bounds for the prediction grid (lat=$(row.latitude), lon=$(row.longitude))."))
            end
            x = pt[1]
            y = pt[2]
            month_year = row.month_year
            hour_ref = row.hour_ref
            voda_pci = pred_pci[i_y, i_x] !== missing ? pred_pci[i_y, i_x] : ""
            voda_rsrp = pred_rsrp[i_y, i_x] !== missing ? pred_rsrp[i_y, i_x] : ""
            ofcom_pci_top1 = (row.earfcn_top1_vf !== missing) &&
                             (row.earfcn_top1_vf == 6300) ?
                             row.pci_top1_vf : ""
            ofcom_pci_top2 = (row.earfcn_top2_vf !== missing) &&
                             (row.earfcn_top2_vf == 6300) ?
                             row.pci_top2_vf : ""
            ofcom_pci_top3 = (row.earfcn_top3_vf !== missing) &&
                             (row.earfcn_top3_vf == 6300) ?
                             row.pci_top3_vf : ""
            ofcom_pci_top4 = (row.earfcn_top4_vf !== missing) &&
                             (row.earfcn_top4_vf == 6300) ?
                             row.pci_top4_vf : ""
            ofcom_rsrp_top1 = (row.earfcn_top1_vf !== missing) &&
                              (row.earfcn_top1_vf == 6300) ?
                              row.rsrp_top1_vf : ""
            ofcom_rsrp_top2 = (row.earfcn_top2_vf !== missing) &&
                              (row.earfcn_top2_vf == 6300) ?
                              row.rsrp_top2_vf : ""
            ofcom_rsrp_top3 = (row.earfcn_top3_vf !== missing) &&
                              (row.earfcn_top3_vf == 6300) ?
                              row.rsrp_top3_vf : ""
            ofcom_rsrp_top4 = (row.earfcn_top4_vf !== missing) &&
                              (row.earfcn_top4_vf == 6300) ?
                              row.rsrp_top4_vf : ""
            println(f,
                "$x,$y,$month_year,$hour_ref,$voda_pci,$voda_rsrp,$ofcom_pci_top1,$ofcom_rsrp_top1,$ofcom_pci_top2,$ofcom_rsrp_top2,$ofcom_pci_top3,$ofcom_rsrp_top3,$ofcom_pci_top4,$ofcom_rsrp_top4")
        end
    end
end

end # module Vodafone
