using PPP
using Test

@testset "PPP.jl" begin
    # Include all test files
    include("ValidateBedogni2024.jl")
    include("ValidateJorner2024.jl")
    #    include("test_parameters.jl")
end
