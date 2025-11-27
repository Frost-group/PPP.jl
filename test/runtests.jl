using PPP
using Test
using StaticArrays
using Logging

# Silent mode (default - only warnings/errors)
global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "PPP.jl" begin
    #include("ValidateJorner2024.jl")
    #include("test_parameters.jl")
    #include("ValidateBedogni2024.jl")
    include("ValidateZhang2011.jl") 
end

