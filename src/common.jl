const Maybe{T} = Union{Nothing,T}

gets(c::Dict, def, k, ks...) = haskey(c, k) ? gets(c[k], def, ks...) : def
gets(c::Tuple, def, k, ks...) = k in 1:length(c) ? gets(c[k], def, ks...) : def
gets(c, _) = c

test_dicts(match::Function, a::Dict, b::Dict) =
    for k in union(keys(a), keys(b))
        @testset "$k" begin
            @test haskey(a, k) && haskey(b, k)
            if haskey(a, k) && haskey(b, k)
                match(k, a[k], b[k])
            end
        end
    end

struct QuietTestSet <: Test.AbstractTestSet
    inner::Test.DefaultTestSet
    QuietTestSet(desc::AbstractString; kwargs...) =
        new(Test.DefaultTestSet(desc; kwargs...))
end

Quiet(args...; kwargs...) = QuietTestSet(Test.DefaultTestSet(args...; kwargs...))

function Test.record(ts::QuietTestSet, t::Union{Test.Result,Test.AbstractTestSet})
    Test.record(ts.inner, t)
end

function Test.record(
    ts::QuietTestSet,
    t::Test.Fail;
    print_result::Bool = Test.TESTSET_PRINT_ENABLE[],
)
    if print_result
        print(ts.inner.description, ": ")
        println()
        if t.test_type !== :test_interrupted
            println(t)
        end
    end
    push!(ts.inner.results, t)
    return t
end

Test.finish(ts::QuietTestSet; kwargs...) = Test.finish(ts.inner; kwargs...)
Test.print_test_errors(ts::QuietTestSet) = Test.print_test_errors(ts.inner)
Test.print_test_results(ts::QuietTestSet, depth_pad = 0) =
    Test.print_test_results(ts.inner, depth_pad)
Test.get_alignment(ts::QuietTestSet, depth = 0) = Test.get_alignment(ts.inner, depth)
Test.filter_errors(ts::QuietTestSet) = Test.filter_errors(ts.inner)
Test.get_test_counts(ts::QuietTestSet) = Test.get_test_counts(ts.inner)
Test.print_counts(ts::QuietTestSet, args...) = Test.print_counts(ts.inner, args...)

macro atest(ex, desc)
    result = quote
        try
            $(Test.Returned)($(esc(ex)), nothing, $(QuoteNode(__source__)))
        catch _e
            _e isa InterruptException && rethrow()
            $(Test.Threw)(_e, Base.current_exceptions(), $(QuoteNode(__source__)))
        end
    end
    Base.remove_linenums!(result)
    quote
        $(Test.do_test)($result, $desc)
    end
end
