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

# Base libraries
import time
import os
import sys
import pathlib

# External libraries
np.set_printoptions(linewidth=190)
# }}}

__author__ = 'Rasmus Liland'
__email__ = 'jens.rasmus.liland@nmbu.no'

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

# {{{ Options: Classifiers
n_countjobs = artmeta["n_countjobs"]
n_jobs = artmeta["n_jobs"]  # My laptop has four threads ...
n_jobs = 1
knn_p = 1  # Power parameter for the Minkowski metric
knn_metrics = ["minkowski", "euclidean", "manhattan", "chebyshev"]
n_estimators = 100  # RandomForest. 100 is the default value.
# }}}


# {{{ Function: BinaryLogRatio
def BinaryLogRatio(x, y):
    return(np.sum(np.absolute(np.log2(y/x))))
# }}}


# {{{ Function: LoadPseudoX
def LoadPseudoX(group,
                K,
                pseudo_count):
    xtrainpath = 'data/ncbi/counts/{}_K={}.csv.xz'
    xtrainpath = xtrainpath.format(group, K)
    X = pd.read_csv(filepath_or_buffer=xtrainpath)

    # ## Adding pseudo-count of 0.5 to all values of
    #    the count matrix before computing the relative
    #    K-mer frequencies, to make sure no value in
    #    the matrix contains 0, to be able to compute
    #    the KNN binary log-ratio distance later.
    #    Pseudo-count of 0.5 was chosen based on 0.5
    #    being in-between the smallest increase of 1
    #    from 0 to 1, for any new K-mer to be counted.
    X = X + pseudo_count

    # ## Compute relative frequencies from the count
    #    matrix by dividing the count matrix by the row
    #    sum in the column (each row is divided by its
    #    sum.
    X = X.div(X.sum(axis=1), axis=0)
    return(X)
# }}}


# {{{ Function: ReadXtest
def ReadXtest(artdir,
              SequencingTechnology,
              D,
              accession,
              group,
              K,
              err):
    PathXtest = "{}/{}/D={}/xtest_{}{}_{}_K={}.csv.xz"
    PathXtest = PathXtest.format(artdir,
                                 SequencingTechnology,
                                 D,
                                 accession,
                                 err,
                                 group,
                                 K)
    return(pd.read_csv(PathXtest))
# }}}


# {{{ Function: IsPlasmid
IsPlasmid = np.vectorize(lambda x: 1 if x == "plasmids" else 0)
# }}}


# {{{ Function: ClassificationRoundTrain
def ClassificationRoundTrain(Y,
                             Y_num,
                             X,
                             ClassToNum,
                             NumToClass,
                             CvFolds,
                             n_neighbors,
                             knn_metrics,
                             n_estimators,
                             pseudo_count):
    # {{{ Prepare classifiers and result headers
    forest = ens.RandomForestClassifier(n_jobs=n_jobs,
                                        n_estimators=n_estimators)
    results = pd.DataFrame(data={'Y': Y})
    stopwatches = {}
    MethodKNN = [('skl.{}NN.{}'.format(knnK, m), m, knnK)
                 for knnK, m in it.product(n_neighbors, knn_metrics)]
    MethodKNN.extend([('skl.{}NN.BinaryLogRatio'.format(knnK),
                       BinaryLogRatio, knnK)
                      for knnK in n_neighbors])

    MethodIdx = list(map(lambda x: x[0], MethodKNN))
    MethodIdx.extend(['skl.RandomForest'])
    # }}}

    # {{{ Cross-validation loop

    # {{{ Prepare cross-validation dataset folds in cv loop
    for i in MethodIdx:
        results['{}.predict'.format(i)] = ''
        results['{}.predict_proba'.format(i)] = ''
        stopwatches[i] = 0
    # }}}

    # {{{ Loop through cv folds
    for test_idx in CvFolds.values:
        train_idx = np.invert(test_idx)
        ytrain = Y[train_idx]
        # ytest = Y[test_idx]
        ynumtrain = Y_num[train_idx]
        # ynumtest = Y_num[test_idx]
        xtrain = X.values[train_idx]
        xtest = X.values[test_idx]
    # }}}

        # {{{ KNN
        for methodname, metric, knnK in MethodKNN:
            knn = nei.KNeighborsClassifier(n_neighbors=knnK,
                                           p=knn_p,
                                           metric=metric,
                                           n_jobs=n_jobs)
            methodtime = time.time()
            knn.fit(xtrain, ytrain)
            arg = knn.predict(xtest)
            results['{}.predict'.format(methodname)][test_idx] = arg
            argproba = flattenproba(knn.predict_proba(xtest))
            results['{}.predict_proba'.format(methodname)][test_idx] = argproba
            stopwatches[methodname] += time.time() - methodtime
        # }}}

        # {{{ Random forest
        methodname = 'skl.RandomForest'
        methodtime = time.time()
        forest.fit(xtrain, ynumtrain)
        arg = np.where(forest.predict(xtest) > 0.5,
                       NumToClass[1], NumToClass[0])
        results['{}.predict'.format(methodname)][test_idx] = arg
        argproba = flattenproba(forest.predict_proba(xtest))
        results['{}.predict_proba'.format(methodname)][test_idx] = argproba
        stopwatches[methodname] += time.time() - methodtime
        # }}}

    # }}}

    watches = pd.DataFrame(stopwatches, index=[0])
    return([results, watches])

# }}}


# {{{ Function: flattenproba
def flattenproba(x):
    return(np.array(list(map(lambda y: ";".join(np.vectorize(np.str)(y)), x))))
# }}}


# {{{ Function: ClassificationRoundTest
def ClassificationRoundTest(Y_train,
                            Y_train_num,
                            X_train,
                            Y_test,
                            Y_test_num,
                            X_test,
                            K,
                            ClassToNum,
                            NumToClass,
                            n_neighbors,
                            knn_metrics,
                            n_estimators):
    # {{{ Prepare classifiers and result headers
    forest = ens.RandomForestClassifier(n_jobs=n_jobs,
                                        n_estimators=n_estimators)
    results = pd.DataFrame(data={'Y': Y_test})
    stopwatches = {}
    MethodKNN = [('skl.{}NN.{}'.format(knnK, m), m, knnK)
                 for knnK, m in it.product(n_neighbors, knn_metrics)]
    MethodKNN.extend([('skl.{}NN.BinaryLogRatio'.format(knnK),
                       BinaryLogRatio, knnK)
                      for knnK in n_neighbors])
    MethodIdx = list(map(lambda x: x[0], MethodKNN))
    MethodIdx.extend(['skl.RandomForest'])
    # }}}

    # {{{ Prepare results and stopwatches dataset folds in cv loop
    for i in MethodIdx:
        results['{}.predict'.format(i)] = ''
        results['{}.predict_proba'.format(i)] = ''
        stopwatches[i] = 0
    # }}}

    # {{{ KNN
    for methodname, metric, knnK in MethodKNN:
        knn = nei.KNeighborsClassifier(n_neighbors=knnK,
                                       p=knn_p,
                                       metric=metric,
                                       n_jobs=n_jobs)
        methodtime = time.time()
        knn.fit(X_train, Y_train)
        arg = knn.predict(X_test)
        results['{}.predict'.format(methodname)] = arg
        argproba = flattenproba(knn.predict_proba(X_test))
        results['{}.predict_proba'.format(methodname)] = argproba
        stopwatches[methodname] += time.time() - methodtime
    # }}}

    # {{{ Random forest
    methodname = 'skl.RandomForest'
    methodtime = time.time()
    forest.fit(X_train, Y_train_num)
    arg = np.where(forest.predict(X_test) > 0.5,
                   NumToClass[1], NumToClass[0])
    results['{}.predict'.format(methodname)] = arg
    argproba = flattenproba(forest.predict_proba(X_test))
    results['{}.predict_proba'.format(methodname)] = argproba
    stopwatches[methodname] += time.time() - methodtime
    # }}}

    watches = pd.DataFrame(stopwatches, index=[0])
    return([results, watches])
# }}}


# {{{ Function: string_chunks
def string_chunks(txt, n):
    length = len(txt)
    average = length/n
    last = 0
    asint = np.vectorize(np.int)
    lst = []
    while last < length:
        lst.append(asint([np.ceil(last), np.ceil(last + average - 1)]))
        last = last + average
    return(pd.DataFrame(lst, columns=["start", "stop"]))
# }}}


# {{{ ReadXSampleN
def ReadXSampleN(x, n):
    return(pd.read_csv(PathXtest(x)).sample(n, replace=False))
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
    # }}}

    # {{{ IsPlasmid column on location and artlocation
    location["IsPlasmid"] = IsPlasmid(location["folder"])
    artlocation["IsPlasmid"] = IsPlasmid(artlocation["folder"])
    # }}}

    # {{{ artlocation derived variables
    genomes = np.unique(artlocation["system"])
    accessions = artlocation["accession"].values
    Yart = artlocation["folder"].values
    # }}}

    # {{{ When a countjob index is the argument being
    #     passed as sys.argv[1], the genome list is
    #     being sliced accordingly.
    if len(sys.argv) > 1:
        countjob_idx = np.int(sys.argv[1]) - 1
        chunks = string_chunks(txt=genomes, n=n_countjobs)
        FUN = lambda x: genomes[np.arange(x[0], x[1]+1)]
        chunks = chunks.apply(FUN, axis=1)
        genomes = chunks[countjob_idx]
    # }}}

    # {{{ genome based cv folds: Similar to what the
    #     genome.chunks function in the knnfun.R
    #     library does: in this folds matrix, the
    #     accessions in the row direction are folded
    #     into the genomes they belong to.
    CvFolds = lambda x: (location["system"] == x).values
    CvFolds = pd.DataFrame(map(CvFolds, np.unique(location["system"])))
    # }}}

    # {{{ Numeric response for RandomForestClassifier

    # ## Define clear dictionaries on what is what
    ClassToNum = {'chromosomes': 0, 'plasmids': 1}
    NumToClass = dict((v, k) for k, v in ClassToNum.items())

    # ## Create dict of numeric responses
    Y_num = np.vectorize(lambda x: ClassToNum[x])(Y)
    # }}}

    # {{{ mesg strings
    mesgresults = 'data/ncbi/sklearn/{}_K={}.csv.xz'
    mesgwatches = 'data/ncbi/sklearn/{}_K={}_time.csv.xz'
    mesgresultsart = 'data/ncbi/sklearn/xtrain_{}_K={}_{}{}_{}.csv.xz'
    mesgwatchesart = 'data/ncbi/sklearn/xtrain_{}_K={}_{}{}_{}_time.csv.xz'
    mesgPathXtest = "{}/{}/D={}/xtest_{}{}_{}_K={}.csv.xz"
    # }}}

    # {{{ Classify cross-validated xtrain data
    #
    #     Integer word length, K.  As K increases
    #     number of features increses, thus runtime are
    #     expected to increase.
#     for K, group in [(2, groups[0])]:
    for K, group in it.product(KLengthsInitial, groups):
        savepathresults = mesgresults.format(group, K)
        savepathwatches = mesgwatches.format(group, K)
        if not os.path.isfile(savepathresults):
            X = LoadPseudoX(group, K, pseudo_count)

            # {{{ Start to time this process
            mesg = '[xtrain] Classifying {} {}-mers ... '
            print(mesg.format(group, K), end='')
            ptm = time.time()  # NB: This will surely not include computing time for threads!
            # }}}

            # {{{ ClassificationRoundTrain
            results, watches = \
                ClassificationRoundTrain(Y=Y,
                                         Y_num=Y_num, X=X,
                                         ClassToNum=ClassToNum,
                                         NumToClass=NumToClass,
                                         CvFolds=CvFolds,
                                         n_neighbors=n_neighbors,
                                         knn_metrics=knn_metrics,
                                         n_estimators=n_estimators,
                                         pseudo_count=pseudo_count)
            # }}}

            # {{{ End of timing this process
            mesg = 'done after {:.2f} seconds.'
            print(mesg.format(time.time() - ptm))
            # }}}

            # {{{ Write classification results to xz compressed
            #     csv file
            results.to_csv(savepathresults, index=False, compression='xz')
            watches.to_csv(savepathwatches, index=False, compression='xz')
            # }}}

    # }}}

    print('ok')
