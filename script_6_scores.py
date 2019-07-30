#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
# vim: ts=8 et sw=4 sts=4 guioptions= modeline

# {{{ Libraries
# Module level imports
import numpy as np
import pandas as pd
import itertools as it
import sklearn.neighbors as nei
import sklearn.ensemble as ens
import multiprocessing as mp

# Base libraries
import time
import os
import sys
import pathlib
import re

# External libraries
np.set_printoptions(linewidth=150)
# }}}

# {{{ Read meta info defined in first mining script
artmeta = "data/ncbi/art/artmeta.csv"
artmeta = pd.read_csv("data/ncbi/art/artmeta.csv", sep="|").T[0]
# }}}

# {{{ artmeta variables
KLengthsInitial = np.vectorize(np.int)(artmeta["K.lengths.initial"].split(";"))
n_neighbors = np.vectorize(np.int)(artmeta["n_neighbors"].split(";"))
pseudo_count = artmeta['pseudo_count']
groups = ('c', 'a')
folder = 'data/ncbi/sklearn'
if not os.path.isdir(folder):
    pathlib.Path(folder).mkdir(parents=True,
                               exist_ok=True)
SequencingDepths = np.vectorize(int)(artmeta["sequencing.depths"].split(";"))
KLengthsART = np.vectorize(int)(artmeta["K.lengths.art"].split(";"))
SequencingTechnology = artmeta["sequencing.technology"]
artdir = artmeta["artdir"]
# }}}

# {{{ Set other iterable vars
Rs = [1, 2]
errs = ["_errFree", ""]
# }}}

__author__ = 'Rasmus Liland'
__email__ = 'jens.rasmus.liland@nmbu.no'

# {{{ Function: Accuracy
def Accuracy(tp, tn, fp, fn):
    # accuracy (ACC)
    return((tn+tp)/(tp+tn+fp+fn))
# }}}

# {{{ Function: TPR
def TPR(tp, fn):
    # sensitivity, recall, hit rate, or true positive rate (TPR)
    return(tp/(tp+fn))
# }}}

# {{{ Function: TNR
def TNR(tn, fp):
    # specificity, selectivity or true negative rate (TNR)
    return(tn/(tn+fp))
# }}}

# {{{ Function: PPV
def PPV(tp, fp):
    # precision or positive predictive value (PPV)
    return(tp/(tp+fp))
# }}}

# {{{ Function: FOne
def FOne(tp, fp, fn):
    # F1 score
    return(2*tp/(2*tp + fp + fn))
# }}}

# {{{ Function: MCC
def MCC(tp, tn, fp, fn):
    # Matthews correlation coefficient
    return((tp*tn - fp*fn)/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))**.5)
# }}}

# {{{ Function: df_acc_tpr_tnr_ppv
def df_acc_tpr_tnr_ppv(df):
    fun = lambda x: \
        pd.Series([Accuracy(tp=x["tp"], tn=x["tn"],
                            fp=x["fp"], fn=x["fn"]),
                   TPR(tp=x["tp"], fn=x["fn"]),
                   TNR(tn=x["tn"], fp=x["fp"]),
                   PPV(tp=x["tp"], fp=x["fp"]),
                   FOne(tp=x["tp"], fp=x["fp"], fn=x["fn"]),
                   MCC(tp=x["tp"], tn=x["tn"],
                       fp=x["fp"], fn=x["fn"])],
                  index=["acc", "tpr", "tnr", "ppv",
                         "fone", "mcc"])
    return(df.apply(fun, axis=1))
# }}}

# {{{ Function: ConfusionLabelLine
def ConfusionLabelLine(predict, y, positive):
    out = np.repeat("TP", len(predict))
    out[np.logical_and(predict==y, predict!=positive)] = "TN"
    out[np.logical_and(predict!=y, predict==positive)] = "FP"
    out[np.logical_and(predict!=y, predict!=positive)] = "FN"
    return(out)
# }}}

# {{{ Function: ConfusionStacked
def ConfusionStacked(df, cn_predict, positive, methods):
    x = df[cn_predict].values.T
    y = df["Y"].values
    tp = (np.logical_and(x==y, y==positive)).sum(axis=1)
    tn = (np.logical_and(x==y, x!=positive)).sum(axis=1)
    fp = (np.logical_and(x!=y, x==positive)).sum(axis=1)
    fn = (np.logical_and(x!=y, x!=positive)).sum(axis=1)
    columns = ["tp", "tn", "fp", "fn"]
    return(pd.DataFrame(np.array([tp, tn, fp, fn]).T,
                        columns=columns))
# }}}

# {{{ Function: ConfusionProba
def ConfusionProba(path, p_idx, pp_idx, positive):
    csv = pd.read_csv(path)
    f = np.vectorize(lambda y: np.float(y.split(";")[1]))
    ff = lambda predict: \
        ConfusionLabelLine(predict=predict, y=csv["Y"],
                           positive=positive)
    confusion = csv.loc[:,p_idx].apply(ff, axis=0)
    return(pd.DataFrame([f(csv.loc[:,pp_idx].values.flatten()),
                         confusion.values.flatten()]).T)
# }}}

# {{{ Function: ReplicateVars
def ReplicateVars(V):
    out, v = V[0]
    methods = V[1]
    dat = pd.DataFrame(np.tile(v, out.shape[0]).reshape(out.shape[0], v.shape[0]))
    meth = np.tile(methods, np.int(out.shape[0]/methods.shape[0]))
    return(pd.DataFrame(np.column_stack((dat, meth))))
# }}}

# {{{ Function: ConfusionStackedART
def ConfusionStackedART(path):
    return(ConfusionStacked(df=pd.read_csv(path),
                            cn_predict=cn_predict,
                            positive=positive,
                            methods=methods))
# }}}

# {{{ Function: MapScores
def MapScores(paths, v, n_jobs):
    with mp.Pool(n_jobs) as pool:
        c = pool.map(ConfusionStackedART, paths)
        c = pd.concat(c, axis=0)
        c.index = np.arange(c.shape[0])
        idx = np.invert(v.duplicated().values)
        m = lambda x: c.loc[(v==x).sum(axis=1)==v.shape[1],].sum(axis=0)
        c = v.iloc[idx,:].apply(m, axis=1)
        art_anova = pd.concat([v.iloc[idx,:], c, df_acc_tpr_tnr_ppv(c)], axis=1)
        return(art_anova)
# }}}

if __name__ == '__main__':
    # {{{ Load response using lzma and pandas
    path = 'data/ncbi/location.csv.xz'
    """ All data files, e.g. response and features, are
        compressed using lzma to read less from disk
        because cpu to memory are faster than reading
        from disk """
    location = pd.read_csv(filepath_or_buffer=path)
    Y = location["folder"].values
    path = 'data/ncbi/artlocation.csv.xz'
    artlocation = pd.read_csv(filepath_or_buffer=path)
    path = 'data/ncbi/nanolocation.csv.xz'
    nanolocation = pd.read_csv(filepath_or_buffer=path)
    # }}}

    # {{{ mesg strings
    mesgresults = 'data/ncbi/sklearn/{}_K={}.csv.xz'
    mesgwatches = 'data/ncbi/sklearn/{}_K={}_time.csv.xz'
    mesgresultsart = 'data/ncbi/sklearn/xtrain_{}_K={}_{}{}_{}.csv.xz'
    mesgwatchesart = 'data/ncbi/sklearn/xtrain_{}_K={}_{}{}_{}_time.csv.xz'
    mesgresultsnano = 'data/ncbi/sklearn/nano_{}_K={}_{}.csv.xz'
    mesgwatchesnano = 'data/ncbi/sklearn/nano_{}_K={}_{}_time.csv.xz'
    mesgPathXtest = "{}/{}/D={}/xtest_{}{}_{}_K={}.csv.xz"

    initialanovasavepath = 'data/ncbi/sklearn/initialanova.csv.xz'
    artanovasavepath = 'data/ncbi/sklearn/artanova.csv.xz'
    nanoanovasavepath = 'data/ncbi/sklearn/nanoanova.csv.xz'

    initialprobasavepath = 'data/ncbi/sklearn/initialproba.csv'
    artprobasavepath = 'data/ncbi/sklearn/artproba.csv'
    nanoprobasavepath = 'data/ncbi/sklearn/nanoproba.csv'
    # }}}

    # {{{ find consensus of genomes existing between locations
    iterator = it.product(np.unique(nanolocation["system"]),
                          groups, KLengthsART)
    v = pd.DataFrame(iterator, columns=["genome", "group", "K"])
    fun = lambda x: mesgresultsnano.format(x["group"],
                                           x["K"],
                                           x["genome"])
    paths = v.apply(fun, axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    nanogenomes = v.loc[exists,"genome"].values

    v = pd.DataFrame(it.product(np.unique(artlocation["system"]),
                                groups, errs, KLengthsART),
                     columns=["genome", "group", "err", "K"])
    #### Sequencing depth; redundancy of coverage
    D = 10
    fun = lambda x: mesgresultsart.format(x["group"],
                                          x["K"],
                                          D,
                                          x["err"],
                                          x["genome"])
    paths = v.apply(fun, axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    artgenomes = v.loc[exists,"genome"].values

    fullsystem = np.unique(location["system"])
    a = np.isin(fullsystem, nanogenomes)
    b = np.isin(fullsystem, artgenomes)
    c = np.logical_and(a, b)

    ### Filtering location based on which genomes exists in 
    #   artlocation and nanolocation, and then using the 50 
    #   first of those genomes further down this pipe.
    location = location.iloc[np.isin(location["system"],
                                     fullsystem[c]),:]
    idx = np.isin(location["system"],
                  np.unique(location["system"])[0:50])
    location = location.iloc[idx,:]
    # }}}

    # {{{ location derived variables
    genomes = np.unique(location["system"])
    accessions = location["accession"].values
    Yart = location["folder"].values
    # }}}

    positive = "plasmids"
    n_jobs = artmeta["n_jobs"]
    D = SequencingDepths[0]
    # {{{ n_jobs can be specified as first argv
    if len(sys.argv) > 1:
        n_jobs = np.int(sys.argv[1])
    # }}}

    """
    methodtime = time.time()
    # {{{ df_anova
    # {{{ methods
    K = KLengthsInitial[0]
    df = pd.read_csv(mesgresults.format(groups[0], K))
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    v = pd.DataFrame(it.product(KLengthsInitial, groups))
    # {{{ paths
    fun = lambda x: mesgresults.format(x[1], x[0])
    paths = v.apply(fun, axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    # }}}
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    oldlen = v.shape[0]

    vv = pd.DataFrame(np.reshape(np.tile(v.values, len(methods)),
                                 (-1, v.shape[1])),
                      columns=["K", "group"])
    vv["Method"] = np.tile(methods, oldlen)
    vv["Canonical"] = \
      np.vectorize(lambda x: "Yes" if x=="c" else "No")(vv["group"])
    vv.index = np.arange(vv.shape[0])
    df_anova = MapScores(paths=paths, v=vv, n_jobs=n_jobs).drop(["group"], axis=1)
    df_anova.to_csv(initialanovasavepath, index=False, compression='xz')
    # }}}
    print("df_anova compilation time: {} seconds".format(time.time() - methodtime))

    methodtime = time.time()
    # {{{ ART anova data
    v = pd.DataFrame(it.product(genomes, groups, errs, KLengthsART),
                     columns=["genome", "group", "err", "K"])
    # {{{ paths
    paths = v.apply(lambda x:
                    mesgresultsart.format(x["group"],
                                          x["K"],
                                          D,
                                          x["err"],
                                          x["genome"]),
                    axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    oldlen = v.shape[0]
    # }}}

    # {{{ methods
    df = pd.read_csv(paths[0])
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    vv = pd.DataFrame(np.reshape(np.tile(v.values, len(methods)),
                                 (-1,v.shape[1])), columns=v.columns)
    vv["Method"] = np.tile(methods, oldlen)
    vv.index = np.arange(vv.shape[0])
    art_anova = MapScores(paths=paths,
                          v=vv.drop(["genome"], axis=1),
                          n_jobs=n_jobs)
    art_anova["Err"] = \
      np.vectorize(lambda x: "No" if x=="_errFree" else "Yes")(art_anova["err"])
    art_anova["Canonical"] = \
      np.vectorize(lambda x: "Yes" if x=="c" else "No")(art_anova["group"])
    art_anova = art_anova.drop(["err", "group"], axis=1)
    art_anova.to_csv(artanovasavepath, index=False, compression='xz')
    # }}}
    print("art_anova compilation time: {} seconds".format(time.time() - methodtime))

    methodtime = time.time()
    # {{{ nano anova table
    v = pd.DataFrame(it.product(genomes, groups, KLengthsART),
                     columns=["genome", "group", "K"])
    # {{{ paths
    paths = v.apply(lambda x:
                    mesgresultsnano.format(x["group"],
                                           x["K"],
                                           x["genome"]),
                    axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    oldlen = v.shape[0]
    # }}}

    # {{{ methods
    df = pd.read_csv(paths[0])
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    vv = pd.DataFrame(np.reshape(np.tile(v.values, len(methods)),
                                 (-1,v.shape[1])), columns=v.columns)
    vv["Method"] = np.tile(methods, oldlen)
    vv.index = np.arange(vv.shape[0])
    nano_anova = MapScores(paths=paths,
                           v=vv.drop(["genome"], axis=1),
                           n_jobs=n_jobs)
    nano_anova["Canonical"] = \
      np.vectorize(lambda x: "Yes" if x=="c" else "No")(nano_anova["group"])
    nano_anova = nano_anova.drop(["group"], axis=1)
    nano_anova.to_csv(nanoanovasavepath, index=False, compression='xz')
    # }}}
    print("nano_anova compilation time: {} seconds".format(time.time() - methodtime))

    methodtime = time.time()
    # {{{ nano proba plot table
    v = pd.DataFrame(it.product(genomes, groups, KLengthsART),
                     columns=["genome", "group", "K"])
    # {{{ paths
    paths = v.apply(lambda x:
                    mesgresultsnano.format(x["group"],
                                           x["K"],
                                           x["genome"]),
                    axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    oldlen = v.shape[0]
    # }}}

    # {{{ methods
    df = pd.read_csv(paths[0])
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    # {{{ vparts
    partlen = 50
    partlens = \
      np.concatenate([[0], np.repeat(partlen, np.int(len(paths)/partlen)),
                [len(paths) % partlen]])
    partlens = np.cumsum(partlens)
    vparts = [v.iloc[partlens[i]:partlens[i+1],] for i in np.arange(partlens.shape[0]-1)]
    # }}}

    header = True
    with mp.Pool(n_jobs) as pool:
        for v in vparts:
            linetime = time.time()
            # {{{ paths
            paths = v.apply(lambda x:
                            mesgresultsnano.format(x["group"],
                                                   x["K"],
                                                   x["genome"]),
                            axis=1).values
            # }}}
            iterable = it.product(paths, [p_idx], [pp_idx], [positive])
            out = pool.starmap(ConfusionProba, iterable)
            vvv = pd.concat(pool.map(ReplicateVars,
                                     it.product(zip(out, v.values), [methods])), axis=0)
            vvv.columns = ["Genome", "group", "K", "Method"]
            vvv.index = np.arange(vvv.shape[0])
            out = pd.concat(out, axis=0)
            out.index = np.arange(out.shape[0])
            out.columns = ["Proba", "Label"]
            out = pd.concat([vvv, out], axis=1)
            out["Canonical"] = \
              np.vectorize(lambda x: "Yes" if x=="c"
                           else "No")(out["group"])
            out = out.drop(["group"], axis=1)
            mode = "w" if header else "a"
            out.to_csv(nanoprobasavepath,
                       index=False,
                       mode=mode,
                       header=header)
            header = False
            msg = "Wrote {} rows to {} after {} seconds ... "
            print(msg.format(out.shape[0], nanoprobasavepath, time.time() - linetime))
        # {{{ Compress outpath
        compressiontime = time.time()
        msg = "xz -zf -T {} -9 {}"
        msg = msg.format(n_jobs, nanoprobasavepath)
        os.system(msg)
        msg = "Done compressing {} after {} seconds."
        print(msg.format(nanoprobasavepath, time.time()-compressiontime))
        # }}}
    # }}}
    print("nano_proba compilation time: {} seconds".format(time.time() - methodtime))

    methodtime = time.time()
    # {{{ art proba plot table
    v = pd.DataFrame(it.product(genomes, groups, errs, KLengthsART),
                     columns=["genome", "group", "err", "K"])
    # {{{ paths
    paths = v.apply(lambda x:
                    mesgresultsart.format(x["group"],
                                          x["K"],
                                          D,
                                          x["err"],
                                          x["genome"]),
                    axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    # }}}

    # {{{ methods
    df = pd.read_csv(paths[0])
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    # {{{ vparts
    partlen = 50
    partlens = \
      np.concatenate([[0], np.repeat(partlen, np.int(len(paths)/partlen)),
                [len(paths) % partlen]])
    partlens = np.cumsum(partlens)
    vparts = [v.iloc[partlens[i]:partlens[i+1],] for i in np.arange(partlens.shape[0]-1)]
    # }}}

    header = True
    with mp.Pool(n_jobs) as pool:
        for v in vparts:
            linetime = time.time()
            # {{{ paths
            paths = v.apply(lambda x:
                            mesgresultsart.format(x["group"],
                                                  x["K"],
                                                  D,
                                                  x["err"],
                                                  x["genome"]),
                            axis=1).values
            # }}}
            iterable = it.product(paths, [p_idx], [pp_idx], [positive])
            out = pool.starmap(ConfusionProba, iterable)
            iterable = it.product(zip(out, v.values), [methods])
            vvv = pd.concat(pool.map(ReplicateVars, iterable), axis=0)
            vvv.columns = ["Genome", "group", "Err", "K", "Method"]
            vvv.index = np.arange(vvv.shape[0])
            out = pd.concat(out, axis=0)
            out.index = np.arange(out.shape[0])
            out.columns = ["Proba", "Label"]
            out = pd.concat([vvv, out], axis=1)
            out["Canonical"] = \
              np.vectorize(lambda x: "Yes" if x=="c"
                           else "No")(out["group"])
            out["Err"] = \
              np.vectorize(lambda x: "No" if x=="_errFree"
                           else "Yes")(out["Err"])
            out = out.drop(["group"], axis=1)
            mode = "w" if header else "a"
            out.to_csv(artprobasavepath,
                       index=False,
                       mode=mode,
                       header=header)
            header = False
            msg = "Wrote {} rows to {} after {} seconds ... "
            print(msg.format(out.shape[0], artprobasavepath, time.time() - linetime))
        # {{{ Compress outpath
        compressiontime = time.time()
        msg = "xz -zf -T {} -9 {}"
        msg = msg.format(n_jobs, artprobasavepath)
        os.system(msg)
        msg = "Done compressing {} after {} seconds."
        print(msg.format(artprobasavepath, time.time()-compressiontime))
        # }}}

    # }}}
    print("art_proba compilation time: {} seconds".format(time.time() - methodtime))

    methodtime = time.time()
    # {{{ initial proba plot table
    v = pd.DataFrame(it.product(KLengthsInitial, groups))
    # {{{ paths
    paths = v.apply(lambda x: mesgresults.format(x[1], x[0]), axis=1).values
    exists = np.vectorize(os.path.isfile)(paths)
    paths = paths[exists]
    v = v.loc[exists,:]
    v.index = np.arange(v.shape[0])
    # }}}

    # {{{ methods
    df = pd.read_csv(paths[0])
    cn = df.columns.values
    p = np.array(list(filter(lambda x: "predict" in x, cn)))
    pp = np.array(list(filter(lambda x: "predict_proba" in x, cn)))
    pp_idx = np.isin(cn, pp)
    p_idx = np.logical_xor(np.isin(cn, p), pp_idx)
    cn_predict = cn[p_idx]
    cn_predict_proba = cn[pp_idx]
    methodsub = np.vectorize(lambda x: re.sub(".predict.*|skl.", "", x))
    methods = methodsub(cn_predict)
    # }}}

    header = True
    with mp.Pool(n_jobs) as pool:
        iterable = it.product(paths, [p_idx], [pp_idx], [positive])
        out = pool.starmap(ConfusionProba, iterable)
        iterable = it.product(zip(out, v.values), [methods])
        vvv = pd.concat(pool.map(ReplicateVars,
                                 iterable), axis=0)
        vvv.columns = ["K", "group", "Method"]
        vvv.index = np.arange(vvv.shape[0])
        out = pd.concat(out, axis=0)
        out.index = np.arange(out.shape[0])
        out.columns = ["Proba", "Label"]
        out = pd.concat([vvv, out], axis=1)
        out["Canonical"] = \
          np.vectorize(lambda x: "Yes" if x=="c"
                       else "No")(out["group"])
        out = out.drop(["group"], axis=1)
        out.to_csv(initialprobasavepath,
                   index=False,
                   header=header)
        # {{{ Compress outpath
        compressiontime = time.time()
        msg = "xz -zf -T {} -9 {}"
        msg = msg.format(n_jobs, initialprobasavepath)
        os.system(msg)
        msg = "Done compressing {} after {} seconds."
        print(msg.format(initialprobasavepath, time.time()-compressiontime))
        # }}}
    # }}}
    print("initial_proba compilation time: {} seconds".format(time.time() - methodtime))
    """

    print('ok')
