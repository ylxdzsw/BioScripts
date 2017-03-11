#!/usr/bin/env julia
using OhMyJulia
using DataFrames
using StatsBase
using PyCall

@pyimport sklearn.ensemble as skensemble
@pyimport sklearn.linear_model as sklinear
@pyimport sklearn.cross_validation as skcv
@pyimport sklearn.metrics as skmetric

const models = (
    (sklinear   , :LogisticRegression          , ()),
    # (sklinear   , :SGDClassifier               , ()),
    # (sklinear   , :PassiveAggressiveClassifier , ()),
    (skensemble , :ExtraTreesClassifier        , ()),
    (skensemble , :GradientBoostingClassifier  , ()),
    (skensemble , :RandomForestClassifier      , ()))

const data = readtable(car(ARGS))

const cmr_features = filter(names(data)) do name
    endswith(string(name), "cmr")
end

const num_features = filter(names(data)) do name
    endswith(string(name), "mutpos")
end

function dummy!(data, var, featlist, f=identity)
    for i in data[var] |> dropna |> map(f) |> unique
        name = Symbol(var, '_', i)
        push!(featlist, name)
        data[name] = Int[!isna(x) && f(x) == i for x in data[var]]
    end
    delete!(data, var)
end

function bootstrap_sampling(len::Integer)
    @with fill(false, len) do p
        p[sample(1:len, round(Int, .9len), replace=false)] = true
    end
end

const other_features = [:tumor_size, :age]

dummy!(data, :gender, other_features, strip)
# dummy!(data, :smoke, other_features)

for feat in other_features
    nas = isna(data[feat])
    data[feat] = data[feat] / maximum(data[!nas, feat])
    data[nas, feat] = mean(data[!nas, feat])
end

show_args(x) = if isempty(x) "" else "($(join(map(x->"$k=$v", x), ", ")))" end

for (package, modelname, args) in models
    for feats in (cmr_features, num_features)
        feats = feats ++ other_features
        X = Array(data[feats])
        y = Int[x>2 for x in data[:stage]]

        model = getfield(package, modelname)(;args...)
        model = model[:fit](X, y)
        prediction = model[:predict](X)
        score_all = skmetric.accuracy_score(y, prediction)
        outlyzers = sign.((y .- .5) .* (prediction .- .5)) .< 0

        importance = try
            model[:feature_importances_]
        catch
            model[:coef_]
        end
        top15 = sortperm(abs(importance[:]), rev=true)[1:15]

        # priori = findfirst(x->startswith(string(x), "TP53_nonsense"), feats)
        # priori in top15 || push!(top15, priori)

        X = X[:, top15]

        model = getfield(package, modelname)(;args...)
        model = model[:fit](X, y)
        score_top15 = skmetric.accuracy_score(y, model[:predict](X))

        accuracy = skcv.cross_val_score(model, X, y, cv=10, scoring="accuracy")
        roc_auc  = skcv.cross_val_score(model, X, y, cv=10, scoring="roc_auc")

        cvoutlyzer, specificity = [], f64[]
        for (train, test) in skcv.KFold(size(X, 1), n_folds=7)
            train, test = train[:]+1, test[:]+1
            model = getfield(package, modelname)(;args...)
            model = model[:fit](X[train, :], y[train])
            prediction = model[:predict_proba](X[test, :])[:, 2]

            outlyzer = sign.((y[test] .- .5) .* (prediction .- .5)) .< 0
            append!(cvoutlyzer, data[test[outlyzer], :id])

            true_positives = sum(y[test])
            for i in sort(prediction)
                if sum(y[test][prediction .> i]) <= .8true_positives
                    s = if length(test) > true_positives
                        count(x->x==0, y[test][prediction .<= i]) / (length(test) - true_positives)
                    else
                        1
                    end
                    push!(specificity, s)
                    break
                end
            end
        end

        boot_accuracy = map(1:1000) do _
            train = bootstrap_sampling(length(y))
            test = !train
            model = getfield(package, modelname)(;args...)
            model = model[:fit](X[train, :], y[train])
            model[:score](X[test, :], y[test])
        end

        println("""
        Model: $modelname$(show_args(args))
        Features:
        $(join(map(i->"  $(feats[i]): $(importance[i])", top15), '\n'))
        Training Set Outlyzers ID: $(any(outlyzers) ? join(data[outlyzers, :id], ", ") : "(None)")
        Training Set Accuracy(All Features): $(score_all)
        Training Set Accuracy(Top15 Features): $(score_top15)
        Cross-Validation Outlyzers ID: $(isempty(cvoutlyzer) ? "(None)" : join(cvoutlyzer, ", "))
        Cross-Validation Accuracy: $(mean(accuracy)) (+/- $(2std(accuracy)))
        Cross-Validation ROC-AUC: $(mean(roc_auc)) (+/- $(2std(roc_auc)))
        Bootstrapping Accuracy: $(mean(boot_accuracy)) (+/- $(2std(boot_accuracy)))
        Specificity when 80% Sensitivity: $(mean(specificity)) (+/- $(2std(specificity)))
        """)
    end
end
