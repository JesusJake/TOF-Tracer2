### CHANGE PARAMETERS HERE

#path = "D:/2605.01-2608.05_IP+NOx_-50/"
#start_filename = "2023-09-30-17h32m15.h5"
#end_filename = "2023-10-03-21h14m59.h5"
#parameter_name = "Inlet" # EoverN, RackTmp, 
#y_label = "Temperature [°C]"

### ----- HERE THE CODING STUFF BEGINS

module PlottingNepheleData
using HDF5
using Dates
using PyPlot
pygui(:tk)

export read_hdf5_data, plot_data, filter_files, find_parameter_location, plot_nephele_data

"""
    read_hdf5_data(path, filename, parameter_location, parameter_index)

Reads the data of interest from an HDF5 file.

Arguments:
- `path::String`: Path to the directory containing the HDF5 file.
- `filename::String`: Name of the HDF5 file.
- `parameter_location::String`: Location of the data of interest within the HDF5 file.
- `parameter_index::Int`: Index of the parameter to extract.

Returns a tuple containing:
- `file::HDF5.File`: HDF5 file object.
- `measurement_data::Array`: Array containing the measurement data.
- `timing_data_readed::Array`: Array containing the timing data.
- `acquisition_start::String`: Timestamp indicating the acquisition start time.

The function opens the specified HDF5 file, reads the measurement data corresponding to the specified parameter index,
as well as the timing data and the acquisition start time attribute, and returns them as separate entities.
"""
function read_hdf5_data(path, filename, parameter_location, parameter_index)

    # open the hdf5 file
    file = h5open(joinpath(path, filename), "r+")

    # read measurement data of hdf5 file
    measurement_data = h5read(joinpath(path, filename), "$(parameter_location)/TwData")[parameter_index, :, :]

    acquisition_start = read_attribute(file, "HDF5 File Creation Time")
    timing_data_readed = read(file["TimingData"])

    return file, measurement_data, timing_data_readed, acquisition_start

end

"""
    convert_measurement_data(a_measurement_dataset)

Converts a measurement dataset into an julia array.

Arguments:
- `a_measurement_dataset::AbstractArray`: The measurement dataset to be converted.

Returns:
- `Array`: julia array containing the measurement data.

The function iterates over each slice of the measurement dataset along the second dimension,
concatenates them into an array, and converts it to a regular Julia array.
"""
function convert_measurement_data(a_measurement_dataset)
    measurement_data = [subarray for subarray in eachslice(a_measurement_dataset, dims=2)]
    measurement_array = hcat(transpose.(measurement_data)...)
    return convert(Array, measurement_array)
end

"""
    convert_timing_array(a_timing_dataset)

Converts a timing dataset into an array.

Arguments:
- `a_timing_dataset::AbstractArray`: The timing dataset to be converted.

Returns:
- `Array`: julia array containing the timing data.

The function extracts the timing data from the 'BufTimes' field of the timing dataset,
reshapes it into an array, and converts it to a regular Julia array.
"""
function convert_timing_array(a_timing_dataset)
    buf_times = a_timing_dataset["BufTimes"]
    timing_array = reshape(buf_times, :, 1)'
    return convert(Array, timing_array)
end

"""
    create_datetime_array(acquisition_start, timing_array)

Creates a datetime array based on the acquisition start time and timing data.

Arguments:
- `acquisition_start::String`: Timestamp indicating the acquisition start time.
- `timing_array::Array`: array containing the timing data.

Returns:
- `Array`: 1D array containing datetime values.

The function parses the acquisition start time string into a DateTime object.
It then calculates datetime values based on the acquisition start time and timing data,
rounding the milliseconds to the nearest whole number and adding them to the acquisition start time.
"""
function create_datetime_array(acquisition_start, timing_array)
    datetime_acquisition_start = Dates.DateTime(acquisition_start, "dd/mm/yyyy HH:MM:SS")
    datetime_array = datetime_acquisition_start .+ Millisecond.(round.(1000 .* timing_array))
    return datetime_array
end

"""
    plot_data(a_measurement_array, a_datetime_array, y_label)

Plots measurement data against datetime values.

Arguments:
- `a_measurement_array::Array`: 1D array containing measurement values.
- `a_datetime_array::Array`: 1D array containing datetime values.
- `y_label::String`: Label for the y-axis.

The function creates a plot with datetime values on the x-axis and measurement values on the y-axis.
Measurement points are represented by markers without connecting lines.
"""
function plot_data(a_measurement_array, a_datetime_array, y_label)
    figure()
    plot(a_datetime_array, a_measurement_array, marker="o", linestyle="", label="Combined Measurements")
    xlabel("Date")
    ylabel(y_label)
    title("Nephele Measurement")
    legend([y_label])
    show()
end

"""
    filter_files(directory_path::String, start_filename::String, end_filename::String)

Filters files in a directory based on start and end filenames.

Arguments:
- `directory_path::String`: Path to the directory containing the files.
- `start_filename::String`: Start filename for filtering.
- `end_filename::String`: End filename for filtering.

Returns:
- `Array{String}`: Array containing filenames within the specified range.

The function retrieves all filenames in the specified directory, sorts them in lexicographical order,
and finds the indices of the start and end filenames. It then extracts filenames within the specified range
and returns them as an array. If the start or end filenames are not found, an empty array is returned.
"""
function filter_files(directory_path::String, start_filename::String, end_filename::String)
    all_files = readdir(directory_path)
    
    # Sort files based on lexicographical order
    sorted_files = sort(all_files)
    
    # Find indices of start and end filenames
    start_index = findfirst(file -> file == start_filename, sorted_files)
    end_index = findfirst(file -> file == end_filename, sorted_files)
    
    # Check if start and end filenames were found
    if start_index === nothing || end_index === nothing
        return String[]
    end
    
    # Extract files within the specified range
    filtered_files = sorted_files[start_index:end_index]
    
    return filtered_files
end

"""
    find_parameter_location(path::AbstractString, filename::AbstractString, parameter_name::AbstractString)

Finds the location and index of the parameter of interest within an HDF5 file generated by Nephele software.

Arguments:
- `path::AbstractString`: Path to the directory containing the HDF5 file.
- `filename::AbstractString`: Name of the HDF5 file.
- `parameter_name::AbstractString`: Name of the parameter to search for.

Returns a tuple containing:
- `parameter_location::Union{String, Nothing}`: Location of the parameter within the HDF5 file. If the parameter is not found, returns `nothing`.
- `parameter_index::Union{Int, Nothing}`: Index of the parameter within the parameter vector. If the parameter is not found, returns `nothing`.

The function opens the specified HDF5 file, searches for the parameter within the "NEPHELE" group,
and returns the location and index of the parameter if found. If the parameter is not found or if
the file does not contain "NEPHELE" or "TimingData" groups, appropriate messages are printed.
"""
function find_parameter_location(path::AbstractString, filename::AbstractString, parameter_name::AbstractString)
    # Open the HDF5 file
    h5file = h5open(joinpath(path, filename), "r")

    try
        # check if the NEPHELE data and TimingData is existent:
        if haskey(h5file, "NEPHELE") && haskey(h5file, "TimingData")
            # Initialize parameter location
            parameter_location = nothing
            # Navigate to the "NEPHELE" group
            nephele_group = h5file["/NEPHELE"]

            # Iterate through subgroups (e.g., "192_168_168_90")
            for subgroup in keys(nephele_group)
                if haskey(nephele_group[subgroup], "TwInfo")
                    twinfo_values = read(nephele_group[subgroup]["TwInfo"])

                    # Check if the parameter is in the vector (substring match)
                    parameter_index = findfirst(x -> contains(parameter_name, x), twinfo_values)

                    if parameter_index !== nothing
                        # Get the location of the parameter
                        parameter_location = string("/NEPHELE/", subgroup)
                        return parameter_location, parameter_index
                    end
                end
            end
            # check if a parameter was found
            if parameter_location === nothing
                println("Parameter not found in file $filename !")
            end
        else
            println("Warning: No Nephele or TimingData in file $filename")
        end
    finally
        # Close the HDF5 file
        close(h5file)
    end
end

"""
    plot_nephele_data(a_path, a_start_filename, a_end_filename, a_parameter_name, a_y_label)

Plots Nephele data from multiple files.

Arguments:
- `a_path::String`: Path to the directory containing the HDF5 files.
- `a_start_filename::String`: Start filename for plotting the data.
- `a_end_filename::String`: End filename for plotting the data.
- `a_parameter_name::String`: Name of the parameter to plot.
- `a_y_label::String`: Label for the y-axis.

The function filters files within the specified range, retrieves parameter locations for each file,
reads measurement and timing data from HDF5 files, converts them into usable arrays, and creates datetime arrays.
It then combines the measurement and datetime data from all files and plots them using the `plot_data` function.
"""
function plot_nephele_data(a_path, a_start_filename, a_end_filename, a_parameter_name, a_y_label)

    filtered_files = filter_files(a_path, a_start_filename, a_end_filename)

    # Initialize arrays to store the combined timing and measurement data
    combined_datetimes = DateTime[]
    combined_measurement = Float64[]

    for a_file in filtered_files

        # get the parameter location and check if data is there
        location_result = find_parameter_location(a_path, a_file, a_parameter_name)

        # check if the find_parameter_location function has reveived any inforamtion
        if location_result !== nothing && length(location_result) == 2
            parameter_location, parameter_index = location_result
            
            file, measurement_dataset, timing_data, acquisition_start = read_hdf5_data(a_path, a_file, parameter_location, parameter_index)
            # convert the data into usable arrays
            measurement_array = convert_measurement_data(measurement_dataset)
            timing_array = convert_timing_array(timing_data)

            datetime_array = create_datetime_array(acquisition_start, timing_array)

            # Append data to the combined arrays
            append!(combined_datetimes, datetime_array)
            append!(combined_measurement, measurement_array)
           
            close(file)
        end
    end

    plot_data(combined_measurement, combined_datetimes, a_y_label)

end

end

# call function here
#PlottingNepheleData.plot_nephele_data(path, start_filename, end_filename, parameter_name, y_label)