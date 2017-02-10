#!/usr/bin/env julia

using OhMyJulia
using IntRangeSets

const gdna, cfdna, bed = ARGS

function output(somatic::Bool, ref::AbstractString, alt::AbstractString, info::AbstractString)
    info = Dict(split(x, '=') for x in split(info, ';'))
    prt(i32(somatic), split(info["ATCGNID"], ',')...,
        info["raw_ref_depth"], info["raw_alt_depth"], info["raw_frequency"],
        info["raw_forward_ref"], info["raw_reverse_ref"],
        info["raw_forward_alt"], info["raw_reverse_alt"],
        get(info, "unique_ref", 0), get(info, "unique_alt", 0), get(info, "unique_af", 0.),
        get(info, "P90", "0%")[1:end-1], get(info, "P98", "0%")[1:end-1])
end

const mask = groupby(car, (x,y)->push!(x, parse(i32, y[2]):parse(i32, y[3])), IntRangeSet{i32},
                     map(x->split(x, '\t'), eachline(bed)), dict=Dict{String, IntRangeSet{i32}}())

const gvar = Set(join(split(x, '\t')[[1,2,4,5]], '\t') for x in eachline(gdna) if !startswith(x, '#'))

for line in eachline(cfdna) @when !startswith(line, '#')
    chr, pos, ref, alt, info = split(line, '\t')[[1, 2, 4, 5, 8]]

    if join((chr, pos, ref, alt), '\t') in gvar
        output(false, ref, alt, info)
    elseif parse(i32, pos) in mask[chr]
        output(true, ref, alt, info)
    end
end
