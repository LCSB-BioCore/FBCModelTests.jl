
mutable struct CountTests <: Test.AbstractTestSet
    passes::Int
    fails::Int
    errs::Int
    broken::Int
    CountTests(_) = new(0, 0, 0, 0)
end

function Test.record(ts::CountTests, child::CountTests)
    ts.passes += child.passes
    ts.fails += child.fails
    ts.errs += child.errs
    ts.broken += child.broken
end

function Test.record(ts::CountTests, res::Test.Pass)
    ts.passes += 1
end

function Test.record(ts::CountTests, res::Test.Fail)
    ts.fails += 1
end

function Test.record(ts::CountTests, res::Test.Error)
    ts.errs += 1
end

function Test.record(ts::CountTests, res::Test.Broken)
    ts.broken += 1
end

function Test.finish(ts::CountTests)
    if Test.get_testset_depth() > 0
        Test.record(Test.get_testset(), ts)
    end
    ts
end
