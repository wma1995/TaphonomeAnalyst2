import numpy as np
import pandas as pd
from pathlib import Path

from typing import Union, List, Any
from pandas.io.parsers import read_csv as _read_csv
from pandas.io.parsers import read_table as _read_txt
import os
import h5py
import shutil
import logging
import warnings
from glob import glob
from shutil import rmtree
from numba import njit, prange
import dask
import dask.array as da
import dask.dataframe as dd
from numpy.random.mtrand import dirichlet
from pandas import DataFrame as DF


def read_txt(file_name: str, T: bool = True, verbose: bool = False, **kwargs):
    '''
    Read general delimited file into DataFrame.

    This a wrapper around pandas' read_table function which adds
    optional parsing of lineage information, and sets some default
    parameter values.

    Note:
    By default the data is transposed!
    To avoid this behavior set the parameter 'T' to False.

    Parameters
    ----------
    file : string
        Path to input file.
    T : bool (default True)
        Indicated whether the produced DataFrame will be transposed.
    verbose : bool (default True)
        Indicated whether to print to screen the parsed table stats.

    Returns
    -------
    table : DataFrame
        Parsed table.
    '''
    # Check file
    file_name = Path(file_name)
    if '.txt' in file_name.name:
        temp = _read_txt(file_name, **kwargs)

    elif '.csv' in file_name.name:
        temp = _read_csv(file_name, **kwargs)
    else:
        raise IOError("ERROR - The file cannot be read.")

    # Validation of the minimum size required-Transposed Matrix
    ncol = min(temp.shape[0], 3)
    nrow = min(temp.shape[1], 3)
    assert (ncol <= 3), "The data size is insufficient to apply the algorithm!"
    assert (nrow <= 3), "The data size is insufficient to apply the algorithm!"

    if T:
        temp = temp.T
        s = """\nFinished parsing table.\nTable dimensions, num_rows: {0} & num_colums: {1}\n""" \
            .format(temp.shape[0], temp.shape[1])
        s += '**** Data has been transposed! ****'
    else:
        s = """\nFinished parsing table.\nTable dimensions, num_rows: {0} & num_colums: {1}\n""" \
            .format(temp.shape[0], temp.shape[1])

    # Verbose
    if verbose:
        print(s)
        return temp
    else:
        return temp


def write_txt(frame: Union[pd.DataFrame, np.ndarray], file_name: Union[str, Path], T: bool = True, **kwargs):
    '''
    Write frame to txt file.

    This a wrapper around pandas' to_csv function which adds
    optional writing of lineage information, and sets some default
    parameter values.

    Note:
    By default the data is transposed!
    To avoid this behavior set the parameter 'T' to False.

    Parameters
    ----------
    file : string
        Path to input file.
    T : bool (default True)
        Indicated whether the produced DataFrame will be transposed.
    '''
    if isinstance(frame, np.ndarray):
        frame = pd.DataFrame(frame)

    file_name = Path(file_name)

    # index
    if T:
        frame.T.to_csv(file_name, **kwargs)
    else:
        frame.to_csv(file_name, **kwargs)


def clean_previous_file(name_file: str) -> None:
    file_name = Path(name_file)
    if file_name.is_file():
        file_name.unlink()


# compute
def normalize(frame: Union[np.ndarray, pd.DataFrame], axis: int = 1):
    '''
    Normalize counts by sample total.

    Parameters
    ----------
    axis : {0, 1}
        0 : normalize each column
        1 : normalize each row

    Returns new instance of same class as input frame.
    '''
    # To do for axis=0
    if axis == 0:
        return frame / frame.sum(axis=0, keepdims=True)
    else:
        # to do for the axis=1
        return frame / frame.sum(axis=1, keepdims=True)


def to_fractions(frame: Union[np.ndarray, pd.DataFrame], method: str = 'dirichlet',
                 p_counts: int = 1, axis: int = 1):
    '''
    Covert counts to fraction using given method.

    Parameters
    ----------
    method : string {'dirichlet' (default) | 'normalize' | 'pseudo'}
        dirichlet - randomly draw from the corresponding posterior
                    Dirichlet distribution with a uniform prior.
                    That is, for a vector of counts C,
                    draw the fractions from Dirichlet(C+1).
        normalize - simply divide each row by its sum.
        pseudo    - add given pseudo count (defualt 1) to each count and
                    do simple normalization.
    p_counts : int/float (default 1)
        The value of the pseudo counts to add to all counts.
        Used only if method is dirichlet
    axis : {0 | 1}
        0 : normalize each column.
        1 : normalize each row.

    Returns
    -------
    fracs: frame/array
        Estimated component fractions.
        Returns new instance of same class as input frame.
    '''
    # Validation
    if isinstance(frame, pd.DataFrame):
        frame = frame.values

    # Define Dirichlet Funtion
    def dirichlet_fun(x):
        a = x + int(p_counts)
        return dirichlet(a)

    # normalize case
    if method == 'normalize':
        fracs = normalize(frame, axis)
        return fracs

    # Dirichlet Case
    elif method == 'dirichlet':
        fracs = np.apply_along_axis(dirichlet_fun, axis, frame)
        return fracs
    else:
        logging.info('Unsupported method "%s"' % method)
        raise ValueError('Unsupported method "%s"' % method)


@njit()
def Mesh(a: int):
    '''simple version of :
    https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html
    '''
    n = len(a)
    A1 = np.repeat(a, n)
    A1 = np.reshape(A1, (-1, n))
    A2 = A1.copy()
    return A1.T, A2


def basis_var(Var_mat, M, V_min: float = 1e-4):
    '''
    Estimate the variances of the basis of the compositional data x.
    Assumes that the correlations are sparse (mean correlation is small).
    The element of V_mat are refered to as t_ij in the SparCC paper.
    '''

    if isinstance(Var_mat, np.ndarray):
        Var_mat = da.from_array(Var_mat)

    if isinstance(M, np.ndarray):
        M = da.from_array(M)

    V_min = V_min
    V_vec = Var_mat.sum(axis=1).compute()
    V_base = da.linalg.solve(M, V_vec)
    basis_variance = da.where(V_base <= 0, V_min, V_base).compute()
    return basis_variance


def C_from_V(Var_mat, V_base):
    '''
    Given the estimated basis variances and observed fractions variation matrix,
    compute the basis correlation & covaraince matrices.
    '''

    Vi, Vj = Mesh(V_base)
    Cov_base = 0.5 * (Vi + Vj - Var_mat)
    C_base = Cov_base / np.sqrt(Vi) / np.sqrt(Vj)
    return C_base, Cov_base


def new_excluded_pair(C: Any, previously_excluded: List = [], th: float = 0.1):
    '''
    Find component pair with highest correlation among pairs that
    weren't previously excluded.
    Return the i,j of pair if it's correlaiton >= than th.
    Otherwise return None.
    '''
    C_temp = np.triu(np.abs(C), 1).copy()  # work only on upper triangle, excluding diagonal

    if len(previously_excluded) > 0:
        C_temp[tuple(zip(*previously_excluded))] = 0

    a = np.unravel_index(np.argmax(C_temp), C_temp.shape)
    cmax = C_temp[a]

    if cmax > th:
        return a
    else:
        return None


@njit(parallel=True)
def variation_mat(frame):
    '''
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    Slower version to be used in case the fast version runs out of memory.
    '''
    k = frame.shape[1]
    V = np.zeros((k, k))

    for i in range(k - 1):
        for j in prange(i + 1, k):
            v = np.log(frame[:, i] / frame[:, j])
            v = np.var(v)
            # Changed-DL ddof=0, set ddof to divide by (n-1),
            # rather than n, thus getting an unbiased estimator (rather than the ML one).
            V[i, j] = v
            V[j, i] = v
    return V


def clr(frame: Union[da.Array, dd.DataFrame, np.ndarray], centrality: str = 'mean', axis: int = 1):
    '''
    Do the central log-ratio (clr) transformation of frame.
    'centraility' is the metric of central tendency to divide by
    after taking the logarithm.

    Parameters
    ----------
    centrality : 'mean' (default) | 'median'
    axis : {0, 1}
        0 : transform each column
        1 : transform each row
    '''
    if isinstance(frame, np.ndarray):
        frame = da.from_array(frame)

    if isinstance(frame, dd.DataFrame):
        frame = frame.to_dask_array()

    frame_temp = da.log(frame)
    if centrality == 'mean':
        v = da.mean(frame_temp, axis=axis, keepdims=True)
        R = frame_temp - v
        return R
    else:
        # centrality is 'median':
        v = da.median(frame_temp, axis=axis, keepdims=True)
        R = frame - v
        return R


def run_clr(frame: np.ndarray):
    '''CLR estimation in the matrix'''
    z = clr(frame)
    Cov_base = da.cov(z, rowvar=0)
    C_base = da.corrcoef(z, rowvar=0)

    return dask.compute(C_base, Cov_base)


def run_sparcc(frame, th: float = 0.1, x_iter: int = 10):
    '''
    Estimate the correlations of the basis of the compositional data f.
    Assumes that the correlations are sparse (mean correlation is small).
    '''
    ## observed log-ratio variances
    Var_mat = variation_mat(frame)
    Var_mat_temp = Var_mat.copy()

    ## Make matrix from eqs. 13 of SparCC paper such that: t_i = M * Basis_Varainces
    D = frame.shape[1]  # number of components
    M = np.ones((D, D)) + np.diag([D - 2] * D)

    ## get approx. basis variances and from them basis covariances/correlations
    V_base = basis_var(Var_mat_temp, M)
    C_base, Cov_base = C_from_V(Var_mat, V_base)

    ## Refine by excluding strongly correlated pairs
    excluded_pairs = []
    excluded_comp = np.array([])

    for xi in range(x_iter):
        # search for new pair to exclude
        to_exclude = new_excluded_pair(C=C_base, th=th, previously_excluded=excluded_pairs)

        if to_exclude is None:  # terminate if no new pairs to exclude
            break
        # exclude pair
        excluded_pairs.append(to_exclude)
        i, j = to_exclude
        M[i, j] -= 1
        M[j, i] -= 1
        M[i, i] -= 1
        M[j, j] -= 1
        # inds = zip(*excluded_pairs)

        inda, indb = np.transpose(excluded_pairs)
        Var_mat_temp[inda, indb] = 0
        Var_mat_temp.T[inda, indb] = 0

        # search for new components to exclude
        nexcluded = np.bincount(np.ravel(excluded_pairs))  # number of excluded pairs for each component
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp = np.where(nexcluded >= D - 3)[0]
        excluded_comp_new = set(excluded_comp) - excluded_comp_prev

        if len(excluded_comp_new) > 0:
            # check if enough components left
            if len(excluded_comp) > D - 4:
                warnings.warn('Too many component excluded. Returning clr result.')
                return run_clr(frame)
            for xcomp in excluded_comp_new:
                Var_mat_temp[xcomp, :] = 0
                Var_mat_temp[:, xcomp] = 0
                M[xcomp, :] = 0
                M[:, xcomp] = 0
                M[xcomp, xcomp] = 1
        # run another sparcc iteration
        V_base = basis_var(Var_mat_temp, M)
        C_base, Cov_base = C_from_V(Var_mat, V_base)

        # set excluded components infered values to nans
        for xcomp in excluded_comp:
            V_base[xcomp] = np.nan
            C_base[xcomp, :] = np.nan
            C_base[:, xcomp] = np.nan
            Cov_base[xcomp, :] = np.nan
            Cov_base[:, xcomp] = np.nan
    return C_base, Cov_base


# <4
def basic_corr(frame, method: str = 'sparcc', th: float = 0.1, x_iter: int = 10):
    '''
    Compute the basis correlations between all components of
    the compositional data f.

    Parameters
    ----------
    frame : array_like
        2D array of relative abundances.
        Columns are counts, rows are samples.
    method : str, optional (default 'SparCC')
        The algorithm to use for computing correlation.
        Supported values: SparCC, clr, pearson, spearman, kendall
        Note that the pearson, spearman, kendall methods are not
        altered to account for the fact that the data is compositional,
        and are provided to facilitate comparisons to
        the clr and sparcc methods.
    th : float,default 0.1
        Exclusion threshold for SparCC,the valid values are 0.0<th<1.0
    x_iter : int,default 10
        Number of exclusion iterations for SparCC.

    Returns
    -------
    C_base: array
        Estimated basis correlation matrix.
    Cov_base: array
        Estimated basis covariance matrix.

    '''
    # Check th
    assert (th > 0 and th < 1.0), "The value must be between 0 and 1"

    method = method.lower()

    k = frame.shape[1]
    ## compute basis variances & correlations
    # if k < 4:
    #     logging.info('Can not detect correlations between compositions of <4 components (%d given)' % k)
    #     raise ValueError('Can not detect correlations between compositions of <4 components (%d given)' % k)
    if method == 'clr':
        C_base, Cov_base = run_clr(frame)
    elif method == 'sparcc':
        C_base, Cov_base = run_sparcc(frame, th=th, x_iter=x_iter)
        tol = 1e-3  # tolerance for correlation range
        if np.max(np.abs(C_base)) > 1 + tol:
            warnings.warn('Sparcity assumption violated. Returning clr result.')
            C_base, Cov_base = run_clr(frame)
    else:
        raise ValueError('Unsupported basis correlation method: "%s"' % method)
    return C_base, Cov_base


def main_alg(frame, method: str = 'sparcc',
             th: float = 0.1,
             x_iter: int = 10,
             n_iter: int = 20,
             norm: str = 'dirichlet',
             log: bool = True,
             path_subdir_cor: str = './',
             path_subdir_cov: str = './',
             savedir: str = 'data',
             verbose: bool = False):
    '''
    The main function to organize the execution of the algorithm and the
    processing of temporary files in hdf5 format.

    Parameters
    ----------
    frame : array_like
        2D array of relative abundances.
        Columns are counts, rows are samples.
    method : str, optional (default 'SparCC')
        The algorithm to use for computing correlation.
        Supported values: SparCC, clr, pearson, spearman, kendall
        Note that the pearson, spearman, kendall methods are not
        altered to account for the fact that the data is compositional,
        and are provided to facilitate comparisons to
        the clr and sparcc methods.
    th : float,default 0.1
        Exclusion threshold for SparCC,the valid values are 0.0<th<1.0
    x_iter : int,default 10
        Number of exclusion iterations for SparCC.
    n_iter : int,default 20
        Number of estimation iteration to average over.
    norm : str,(dirichlet|norm),defualt: dirichlet
        Method used to normalize the counts to fractions.
    log : bool, default True
        log-transform fraction? used if method ~= SparCC/CLR
    path_subdir_cor:str,default './'
        Folder path for the temporary correlation estimates file.
    path_subdir_cov:str,default './'
        Folder path for the temporary covariance estimates file
    verbose : bool, default True

    Returns
    -------
    C_base: array
        Estimated basis correlation matrix.
    Cov_base: array
        Estimated basis covariance matrix.

    '''

    if method in ['sparcc', 'clr']:
        for i in range(n_iter):
            if verbose: print('\tRunning iteration ' + str(i))
            logging.info("Running iteration {}".format(i))
            fracs = to_fractions(frame, method=norm)
            cor_sparse, cov_sparse = basic_corr(fracs, method=method, th=th, x_iter=x_iter)
            var_cov = np.diag(cov_sparse)
            # Create files

            file_name_cor = path_subdir_cor + '/cor_{:08d}.hdf5'.format(i)
            file_name_cov = path_subdir_cov + '/cov_{:08d}.hdf5'.format(i)

            h5f_cor = h5py.File(file_name_cor, 'w')
            h5f_cov = h5py.File(file_name_cov, 'w')
            h5f_cor.create_dataset('dataset', data=cor_sparse, shape=cor_sparse.shape)
            h5f_cov.create_dataset('dataset', data=var_cov, shape=var_cov.shape)
            h5f_cor.close()
            h5f_cov.close()

        logging.info("Processing the files from data directory")
        filenames_cor = sorted([f for f in glob(f'{savedir}' + "**/corr_files/*", recursive=True)])
        filenames_cov = sorted([f for f in glob(f'{savedir}' + "**/cov_files/*", recursive=True)])
        dsets_cor = [h5py.File(filename, mode='r') for filename in filenames_cor]
        dsets_cov = [h5py.File(filename, mode='r') for filename in filenames_cov]
        arrays_cor = [da.from_array(dset['dataset']) for dset in dsets_cor]
        arrays_cov = [da.from_array(dset['dataset']) for dset in dsets_cov]

        cor_array = da.asarray(arrays_cor)
        cov_array = da.asarray(arrays_cov)
        var_med = da.nanmedian(cov_array, axis=0).compute()
        cor_med = da.nanmedian(cor_array, axis=0).compute()

        var_med = var_med

        x, y = Mesh(var_med)
        cov_med = cor_med * x ** 0.5 * y ** 0.5
        logging.info("The main process has finished")

        return cor_med, cov_med


def clean_data_folder(path_folder: Union[str, Path]) -> Any:
    path_folder = Path(path_folder)

    if path_folder.is_dir():
        rmtree(path_folder.resolve(), ignore_errors=True)
    else:
        raise EOFError("This is not a directory")


# bootstrapping
def permute_w_replacement(frame: Union[pd.DataFrame, np.ndarray], axis=0):
    '''
    Permute the frame values across the given axis.
    Create simulated dataset were the counts of each component (column)
    in each sample (row), are randomly sampled from the all the
    counts of that component in all samples.

    Parameters
    ----------
    frame : Numpy Array
        Frame to permute.
    axis : {0, 1}
        - 0 - Permute row values across columns
        - 1 - Permute column values across rows

    Returns
    -------
    Permuted DataFrame (new instance).
    '''

    if isinstance(frame, pd.DataFrame):
        frame = frame.values

    fp = lambda x: np.random.permutation(x)

    if axis == 0:
        # Along the columns
        aux = np.apply_along_axis(fp, 0, frame)
        aux = np.apply_along_axis(fp, 1, aux)
        return aux

    elif axis == 1:
        # Along the rows
        aux = np.apply_along_axis(fp, 1, frame)
        aux = np.apply_along_axis(fp, 0, aux)
        return aux


def make_bootstraps(counts: pd.DataFrame, nperm: int, perm_template: str, outpath: str = './', iprint: int = 0):
    '''
    Make n simulated datasets used to get pseudo p-values.
    Simulated datasets are generated by assigning each OTU in each sample
    an abundance that is randomly drawn (w. replacement) from the
    abundances of the OTU in all samples.
    Simulated datasets are either written out as txt files.

    Parameters
    ----------
    counts : DataFrame
        Inferred correlations whose p-values are to be computed.
    nperm : int
        Number of permutations to produce.
    perm_template : str
        Template for the permuted data file names.
        Should not include the path, which is specified using the
        outpath parameter.
        The iteration number is indicated with a "#".
        For example: 'permuted/counts.permuted_#.txt'
    outpath : str (default './')
        The path to which permuted data will be written.
        If not provided files will be written to the cwd.
    iprint : int (default = 0)
        The interval at which iteration number is printed out.
        If iprint<=0 no printouts are made.
    '''
    if not os.path.exists(outpath): os.makedirs(outpath)

    for i in range(nperm):
        if iprint > 0:
            if not i % iprint: print(i)
        # New Matrix
        counts_perm = permute_w_replacement(counts, axis=1)

        outfile = outpath + perm_template.replace('#', '%d' % i)
        # The output is written
        write_txt(counts_perm, outfile, index=True)


# print
def MakeBoots(counts_file: Union[str, Path], nperm: int,
              perm_template: str, outpath: str = './'):
    '''
    Make n simulated datasets used to get pseudo p-values.
    Simulated datasets are generated by assigning each OTU in each sample
    an abundance that is randomly drawn (w. replacement) from the
    abundances of the OTU in all samples.
    Simulated datasets are either written out as csv files.
    '''
    if perm_template is None:
        perm_template = counts_file + '.permuted_#.csv'
    ## read counts data
    counts = read_txt(counts_file, index_col=0)
    if counts.shape[0] == 0:
        print('A problem has occurred with the file, it will be resolved.')

        try:
            counts = read_txt(counts_file, sep=',', index_col=0)
        except IOError as IOE:
            raise (IOE)
    assert counts.shape[0] != 0, "ERROR!"

    ## make permutated data
    # print("Permutations")
    make_bootstraps(counts, nperm, perm_template, outpath=outpath)


# pvalues
def compare2sided(perm, real):
    return np.abs(perm) >= np.abs(real)


def compare1sided(perm, real):
    inds_abs = compare2sided(perm, real)
    inds_sign = np.sign(perm) == np.sign(real)
    return inds_abs & inds_sign


def get_pvalues(cor: DF, perm_template: str, nperm: int,
                test_type: str = 'two_sided', iprint: int = 0):
    '''
    Compute pseudo p-vals from a set correlations obtained from permuted data'
    Pseudo p-vals are the percentage of times a correlation at least
    as extreme as the "real" one was observed in simulated datasets.

    Files containing the permuted correlations should be named with a
    consistent template, and these file names cannot contain any "#" characters.

    Parameters
    ----------
    cor : DataFrame
        Inferred correlations whose p-values are to be computed.
    perm_template : str
        The template used for naming the correlation files of the
        permuted data. The iteration number is indicated with a "#".
        For example: 'permuted/cor.sparcc.permuted_#.txt'
    nperm : int
        Number of permutations available.
    test_type : 'two_sided' (default) | 'one_sided'
        two-sided  = considering only the correlation magnitude.
        one-sided  = accounting for the sign of correlations.
    iprint : int (default = 0)
        The interval at which iteration number is printed out.
        If iprint<=0 no printouts are made.

    Returns
    -------
    p_vals: frame
        Computed pseudo p-values.
    '''
    # Definition of the type of test
    if test_type == 'two_sided':
        cmpfun = compare2sided
    elif test_type == 'one_sided':
        cmpfun = compare1sided
    else:
        raise ValueError('unsupported test type "%s"' % test_type)

    # DataFrame
    n_sig = DF(np.zeros(cor.shape),
               index=cor.index,
               columns=cor.columns)

    for i in range(nperm):
        if iprint > 0:
            if not i % iprint: print(i)
        permfile = perm_template.replace('#', '%d' % i)
        cor_perm = read_txt(permfile, index_col=0).values  # Read each file
        n_sig[cmpfun(cor_perm, cor)] += 1

    p_vals = 1. * n_sig / nperm
    p_vals.values[np.diag_indices_from(p_vals.values)] = 1

    return p_vals


def PseudoPv(cor_file: Union[str, Path], perm_template: str, nperm: int,
             test_type: str = 'two_sided', outfile: Any = None):
    '''
    Compute pseudo p-vals from a set correlations obtained from permuted data'
    Pseudo p-vals are the percentage of times a correlation at least
    as extreme as the "real" one was observed in simulated datasets.

    Files containing the permuted correlations should be named with a
    consistent template, and these file names cannot contain any "#" characters.
    '''
    cor = read_txt(cor_file, verbose=False, index_col=0)
    if cor.shape[0] == 0:
        print('A problem has occurred with the file, it will be resolved.')

    assert cor.shape[0] != 0, "ERROR!"

    print(f"Shape of Corr:{cor.shape}")

    p_vals = get_pvalues(cor, perm_template, nperm, test_type)
    if outfile is None:
        outfile = cor_file + '.nperm_%d.pvals' % nperm

    write_txt(p_vals, outfile)


class SparCC_MicNet:

    def __init__(self,
                 name: str = 'experiment_sparCC',
                 method: str = 'sparcc',
                 low_abundance: bool = False,
                 n_iteractions: int = 2,
                 x_iteractions: int = 2,
                 threshold: float = 0.1,
                 normalization: str = 'dirichlet',
                 log_transform: bool = True,
                 save_corr_file: Union[str, Path] = 'sparcc/example/cor_sparcc.csv',
                 save_cov_file: Union[str, Any] = None,
                 num_simulate_data: int = 3,
                 perm_template: str = 'permutation_#.csv',
                 outpath: Union[str, Path] = 'sparcc/example/pvals/',
                 type_pvalues: str = 'one_sided',
                 outfile_pvals: Union[str, Path] = f'sparcc/example/pvals/pvals_',
                 name_output_file: str = 'sparcc_output'
                 ) -> None:

        self.name = name
        self.method = method
        self.n_iteractions = n_iteractions
        self.x_iteractions = x_iteractions
        self.threshold = threshold
        self.normalization = normalization
        self.log_transform = log_transform
        self.save_corr_file = save_corr_file
        self.save_cov_file = save_cov_file
        self.num_simulate_data = num_simulate_data
        self.perm_template = perm_template
        self.outpath = outpath
        self.type_pvalues = type_pvalues
        self.outfile_pvals = outfile_pvals + type_pvalues + '.csv'
        self.name_output_file = name_output_file
        self.low_abundance = low_abundance

    # print
    def compute(self, data_input: Union[str, pd.DataFrame] = 'sparcc/example/fake_data.txt',
                save_corr_file: str = 'sparcc/example/cor_sparcc.csv'):

        self._preprocess()

        if type(data_input) == str:
            self.data_input = data_input

        else:
            data_input.to_csv('temp_sample_output.csv')
            self.data_input = 'temp_sample_output.csv'

        # Load the file
        try:
            L1 = read_txt(self.data_input, index_col=0)
        except IOError as IOE:
            raise (IOE)

        if L1.shape[0] == 0:
            try:
                L1 = read_txt(self.data_input, sep=',')
            except IOError as IOE:
                raise (IOE)
        assert L1.shape[0] != 0, "ERROR!"

        # SparCC Algorithm
        cor, cov = main_alg(frame=L1,
                            method=self.method,
                            norm=self.normalization,
                            n_iter=self.n_iteractions,
                            verbose=False,
                            log=self.log_transform,
                            th=self.threshold,
                            x_iter=self.x_iteractions,
                            path_subdir_cor=self.path_corr_file,
                            path_subdir_cov=self.path_cov_file,
                            savedir=self.savedir)

        # print("Shape of Correlation Matrix:", cor.shape)
        # print("Shape of Covariance Matrix:", cov.shape)

        # Save Correlation
        write_txt(frame=cor, file_name=save_corr_file)

        clean_data_folder(path_folder=self.path_corr_file)
        clean_data_folder(path_folder=self.path_cov_file)

    def bootstrapping(self, data_input: Union[str, pd.DataFrame] = 'sparcc/example/fake_data.txt'):

        if not hasattr(self, 'data_input'):
            if type(data_input) == str:
                self.data_input = data_input

            else:
                data_input.to_csv('temp_sample_output.csv')
                self.data_input = 'temp_sample_output.csv'

        MakeBoots(counts_file=self.data_input, nperm=self.num_simulate_data,
                  perm_template=self.perm_template, outpath=self.outpath)

    def pvalues(self, save_corr_file: str = 'sparcc/example/cor_sparcc.csv', template: str = ''):

        PseudoPv(save_corr_file, template,
                 self.num_simulate_data, self.type_pvalues,
                 self.outfile_pvals)

    # print
    def run_all(self, data_input: Union[str, pd.DataFrame] = 'sparcc/example/fake_data.txt',
                save_corr_file='Cor_SparCC.csv'):

        clean_previous_file(save_corr_file)
        clean_previous_file('temp_sample_output.csv')
        clean_previous_file(self.outfile_pvals)

        data_input = self._validation_format(data_input)
        data_input = self._filter_otus(data_input, self.low_abundance)

        self.compute(data_input, save_corr_file)
        self.bootstrapping(data_input='')
        file_pvals = str(self.outpath) + str(self.perm_template)
        corr_file = str(self.outpath) + 'perm_cor_#.csv'

        # Iteraction in files
        for i in range(int(self.num_simulate_data)):
            # print('#' * 100)
            # print(f'Iteration: {str(i)}')
            file_pvals1 = file_pvals.replace('#', str(i))
            corr_file1 = corr_file.replace('#', str(i))

            self.compute(data_input=file_pvals1, save_corr_file=corr_file1)

        # Clean Files: perm_corr
        folder = Path(self.outpath)
        List_files_rm = list(folder.glob(self.perm_template.replace('#', '*')))
        for f in List_files_rm:
            f.unlink()

        self.pvalues(save_corr_file, template=corr_file)

        # Clean perm_cor_* files
        List_files_rm = list(folder.glob('perm_cor_*.csv'))
        for f in List_files_rm:
            f.unlink()

        setattr(self, 'save_corr_file', save_corr_file)

        # move log
        # if Path(save_corr_file).parent.is_dir():
        #     output=Path(save_corr_file).parent/f'{self.name}.log'
        #     source=Path(f'sparcc/example/{self.name}.log')
        #     shutil.move(source,output)

        # self.save_corr_file=save_corr_file

    def _preprocess(self):
        # Define Temp Folder
        setattr(self, 'savedir', './sparcc/data')

        # self._check_save_files()

        if os.path.exists(self.savedir) and os.path.isdir(self.savedir):
            try:
                shutil.rmtree("./sparcc/data")

            except OSError as e:
                print("Error: {0}:{1}".format('./temp_files/*', e.strerror))

            self.path_corr_file = os.path.join(str(self.savedir), 'corr_files')
            self.path_cov_file = os.path.join(str(self.savedir), 'cov_files')

            if not os.path.exists(self.path_corr_file):
                os.makedirs(self.path_corr_file)

            if not os.path.exists(self.path_cov_file):
                os.makedirs(self.path_cov_file)

        else:
            os.makedirs(self.savedir)
            self.path_corr_file = os.path.join(str(self.savedir), 'corr_files')
            os.makedirs(self.path_corr_file)
            self.path_cov_file = os.path.join(str(self.savedir), 'cov_files')
            os.makedirs(self.path_cov_file)

    def _validation_format(self, frame: pd.DataFrame) -> pd.DataFrame:

        # OTUS Values
        if all((frame.iloc[:, 0:].dtypes != 'object').values):
            values_df = frame.iloc[:, 0:].values.copy()
        else:
            raise ValueError('There is something wrong with the OTUS values')

        DF_out = pd.DataFrame(index=frame.index, data=values_df)
        return DF_out

    # <2
    def _filter_otus(self, frame: pd.DataFrame, low_abundance: bool = False) -> pd.DataFrame:

        # Remove singletons
        # frame = frame.loc[(frame != 0).sum(axis=1) >= 2, :].copy()
        # Remove low abudance < 5
        if low_abundance == True:
            frame = frame.loc[frame.sum(axis=1) > 5, :].copy()
        else:
            pass

        self._Index_col = frame.index

        return frame

    def __repr__(self) -> str:
        return f"SparCC with n_iteractions:{str(self.n_iteractions)} and x_iteractions: {str(self.x_iteractions)}"

    def __str__(self) -> str:
        return f"SparCC with n_iteractions: {str(self.n_iteractions)} and x_iteractions: {str(self.x_iteractions)}"
