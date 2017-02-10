using OhMyJulia
using LightGBM
using HDF5

const X, y = h5open("data.hdf5", "r") do f
    read(f, "X"), read(f, "y")
end

const model = LGBMBinary(num_iterations = 256,
                         learning_rate = .05,
                         early_stopping_round = 5,
                         bagging_fraction = .9,
                         bagging_freq = 1,
                         num_leaves = 127,
                         metric = ["auc", "binary_error"])

# const splits = let
#     l = length(y)
#     q1, q2, q3 = (round(Int, x*l) for x in (.25, .5, .75))
#     map(collect, (1:q1, q1:q2, q2:q3, q3:l))
# end

# cv(model, X, y, splits)

const split = round(Int, .8length(y))

fit(model, X[1:split, :], y[1:split])

const tX, ty = X[split+1:end, :], y[split+1:end]

const d = tX[:, end-10] .+ tX[:, end-11]
const maxd = maximum(d)

for (stop, start) in Δ(tuple, [-1, .1maxd, .2maxd, .3maxd, .4maxd, .5maxd, .6maxd, maxd])
    ind = start .< d .<= stop
    p = predict(model, tX[ind, :])
    prt(round(Int, start), round(Int, stop), sumabs(p .- ty[ind]) / sum(ind))
end
