include("Keras.jl")

using BinBlocks

const model = let
    input = Keras.Input(shape=(256, 64, 8))

    output = input |>
        Keras.Convolution2D(64, 5, 5, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Convolution2D(128, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((4, 1), strides=(4, 1)) |>

        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Flatten() |>
        Keras.Dense(1024, activation="sigmoid") |>
        Keras.Dense(1, activation="sigmoid")

    Keras.Model(input=input, output=output)
end

const X, y = let
    data = []
    for sample in readdir(".") |> filter(x->endswith(x, ".binblock")) |> map(BinBlock)
        while true
            try
                push!(data, read(sample, 1024))
            catch e
                isa(e, EOFError) ? break : rethrow(e)
            end
        end
    end
    [map(car, data)...;], [map(cadr, data)...;]
end

const callbacks = [Keras.ModelCheckpoint("weights.{epoch:02d}-{val_acc:.4f}.h5", monitor="val_acc", save_weights_only=true)]

const phase_one = readdir(".") |> filter(x->startswith(x, "weights.14"))

if !isempty(phase_one)
    model[:compile](Keras.SGD(lr=.0005, momentum=.95), "binary_crossentropy", metrics=["accuracy"])
    model[:load_weights](phase_one[])
    model[:fit](X, y, batch_size=256, nb_epoch=20, validation_split=0.01, callbacks=callbacks, initial_epoch=15)
else
    model[:compile](Keras.SGD(lr=.002, momentum=.9), "binary_crossentropy", metrics=["accuracy"])
    model[:fit](X, y, batch_size=256, nb_epoch=15, validation_split=0.01, callbacks=callbacks)
end


