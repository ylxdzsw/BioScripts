include("../../NNBuilder/server/model/wrapper/Keras.jl")
using OhMyJulia
using HDF5

const X, y = h5open("data.hdf5", "r") do f
    X = read(f, "X")
    y = read(f, "y")
    scale = maximum(abs(X), 1)
    for j in 1:ncol(X), i in 1:nrow(X)
        X[i, j] /= scale[j]
    end
    X, y
end

const model = let
    input = Keras.Input(shape=(40,))
    fc1   = Keras.Dense(90, activation="sigmoid")(input)
    fc2   = Keras.Dense(60, activation="sigmoid")(fc1)
    fc3   = Keras.Dense(40, activation="sigmoid")(fc2)
    fc4   = Keras.Dense(1, activation="sigmoid")(fc3)
    Keras.Model(input=input, output=fc4)
end

model[:compile](optimizer="SGD", loss="binary_crossentropy", metrics=["accuracy"])
model[:fit](X, y, nb_epoch=20, batch_size=256, validation_split=.3)
