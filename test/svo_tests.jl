# General SVO usage is covered in doctests. Here we test additional special cases, other things not covered in doctests.

import PhotometricFilters
import Unitful
using Test

@testset "SVO special cases" begin
    # Test that Roman/WFI filters are correctly converted from effective area to throughput
    @testset "Roman/WFI filters" begin
        roman_wfi_filters = ["Roman/WFI.F062", "Roman/WFI.F087", "Roman/WFI.F106", "Roman/WFI.F129", "Roman/WFI.F146", "Roman/WFI.F158", "Roman/WFI.F184", "Roman/WFI.F213"]
        for filter_name in roman_wfi_filters
            filt = PhotometricFilters.get_filter(filter_name)
            @test filt isa PhotometricFilters.SVOFilter
            th = PhotometricFilters.throughput(filt)
            for i in th
                @test Unitful.unit(i) == Unitful.NoUnits
            end
            @test all(>=(0), th)
            @test all(<=(1), th)
        end
    end
end