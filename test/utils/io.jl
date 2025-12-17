using HDF5

"""
    test_h5_supported_types()

Test that all HDF5-supported types are indeed supported by HDF5.

This function checks that `String` and all subtypes of `HDF5_SUPPORTED_NUMBER_TYPES`
    are supported by the underlying `HDF5` writer.
Tests are also conducted for 1-, 2- and 3-dimensional `Array` of those types.
"""
function test_h5_supported_types()
    supported_number_types = Base.uniontypes(PGLearn.HDF5_SUPPORTED_NUMBER_TYPES)
    supported_types = [String; supported_number_types]

    # Check that all HDF5-supported types are indeed supported
    # Saving to an HDF5 should not error nor throw any warning
    @testset "$T" for T in supported_types
        f5 = tempname()

        # First check scalars
        d = Dict("$T" => one(T))
        @test_logs PGLearn.save_h5(f5, d)
        d_ = HDF5.h5read(f5, "/")
        @test d == d_

        # Now do the same thing with arrays (up to 3-dimensional arrays)
        for N in 1:3
            v = Dict("$T" => ones(T, NTuple{N}(ones(Int, N))))
            @test_logs PGLearn.save_h5(f5, v)
            h = HDF5.h5read(f5, "/")
            @test h == v
        end
    end
    return nothing
end

"""
    test_h5_precision_warning()

Test that unsupported-but-convertible number types are handled as expected.
"""
function test_h5_precision_warning()
    f5 = tempname()

    x = 3.141592653589 * one(BigFloat)
    z = Complex(x) + im

    msg = """Unsupported data type for writing to an HDF5 group: \"x\"::BigFloat.
    This value was converted to Float64, which may incur a loss of precision."""
    @test_logs (:warn, msg) PGLearn.save_h5(f5, Dict("x" => x))
    h = HDF5.h5read(f5, "/")
    @test h["x"] == convert(Float64, x)

    msg = """Unsupported data type for writing to an HDF5 group: \"z\"::Complex{BigFloat}.
    This value was converted to ComplexF64, which may incur a loss of precision."""
    @test_logs (:warn, msg) PGLearn.save_h5(f5,  Dict("z" => z))
    h = HDF5.h5read(f5, "/")
    @test h["z"] == convert(Complex{Float64}, z)

    msg = """Unsupported data type for writing to an HDF5 group: \"v\"::Vector{BigFloat}.
    This value was converted to Vector{Float64}, which may incur a loss of precision."""
    @test_logs (:warn, msg) PGLearn.save_h5(f5, Dict("v" => [x]))
    h = HDF5.h5read(f5, "/")
    @test h["v"] == convert(Vector{Float64}, [x])

    msg = """Unsupported data type for writing to an HDF5 group: \"w\"::Vector{Complex{BigFloat}}.
    This value was converted to Vector{ComplexF64}, which may incur a loss of precision."""
    @test_logs (:warn, msg) PGLearn.save_h5(f5,  Dict("w" => [z]))
    h = HDF5.h5read(f5, "/")
    @test h["w"] == convert(Vector{Complex{Float64}}, [z])

    rm(f5)  # delete temp file

    return nothing
end

"""
    test_h5_invalid_types()

Test that unsupported HDF5 types throw an error.
"""
function test_h5_invalid_types()
    f5 = tempname()
 
    x = 1:4
    msg = "Unsupported data type for writing to an HDF5 group: \"key\"::UnitRange{Int64}."
    @test_throws msg PGLearn.save_h5(f5, Dict("key" => x))

    # Not converting JuMP status codes to String is a common mistake...
    st = MOI.OPTIMAL
    msg = "Unsupported data type for writing to an HDF5 group: \"termination_status\"::TerminationStatusCode."
    @test_throws msg PGLearn.save_h5(f5, Dict("termination_status" => st))

    return nothing
end

"""
    test_h5_io()

Test that a non-trivial dictionary gets written to HDF5, then parsed correctly.
"""
function test_h5_io()
    f5 = tempname()

    d = Dict(
        "a" => 1,
        "b" => "hello",
        "c" => [3.0, 4.1, 5.2, 6.3],
        "d" => Dict(
            "e" => exp(1),
            "f" => Float32(π),
            "g" => Dict(
                "h" => ones(Float32, 1, 2, 3),
                "i" => im,
                "j" => 1.0 * im,
                "k" => 1 + 1.0 * im,
            ),
        ),
    )
    PGLearn.save_h5(f5, d)
    h = HDF5.h5read(f5, "/")
    @test h == d
end

function test_json_io()
    d = Dict("a" => "a", "b" => 1, "c" => 2.0, "dd" => [3.0, 4.1])
    fjson = tempname()

    @testset "$ext" for ext in [".json", ".json.gz", ".json.bz2"]
        PGLearn.save_json(fjson * ext, d)
        @test isfile(fjson * ext)
        @test PGLearn.load_json(fjson * ext) == d
        rm(fjson * ext)
    end

    # Check that non-JSON extensions are not OK
    @testset "$ext" for ext in [".h5", ".json.xy"]
        msg = """
        Unsupported JSON extension: \"$(basename(fjson*ext))\"
        Supported extensions are: \".json\", \".json.bz2\" and \".json.gz\".
        """
        @test_throws msg PGLearn.save_json(fjson * ext, d)
        @test !isfile(fjson * ext)
        @test_throws msg PGLearn.load_json(fjson * ext)
    end

    return nothing
end

@testset "I/O" begin
    @testset "HDF5" begin
        @testset test_h5_supported_types()
        @testset test_h5_precision_warning()
        @testset test_h5_invalid_types()
        @testset test_h5_io()
    end
    @testset "JSON" begin
        @testset test_json_io()
    end
end
