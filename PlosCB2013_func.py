#This script is used to generate the set of representative positive Boolean function up to n=7
#Functions are generated level by level. Functions of a given level l have l terms/clauses in their DNF/CNF
#For n=6 everything can be run locally, however, for n=7 one needs a cluster or a very powerful machine
import numpy as np
import h5py
from itertools import permutations

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
    1-D arrays to form the cartesian product of.
    out : ndarray
    Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
    2-D array of shape (M, len(arrays)) containing cartesian products
    formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def iptgen(n=4,ex=[2]):
    """
    Generate an input subspace containing vectors with a number of active inputs listed in ex.

    Parameters
    ----------
    n : an integer, the input dimentionality
    ex : a list, to select input subsets composed of input vectors containing the same number of 1s

    Returns
    -------
    A list of arrays, the size depend on n and ex. There are combination(n,ex[0])+combination(n,ex[1])+... arrays
    each array have is of length n. Arrays are in organized from smaller to higher binary number.

    Examples
    --------
    >>> iptgen(n=4, ex=[2])
    array([[0,0,1,1],
           [0,1,0,1],
           [0,1,1,0],
           [1,0,0,1],
           [1 0 1 0],
           [1 1,0,0]])

    """
    ipt = cartesian(np.repeat([[0,1]],n,axis=0))
    tr = np.zeros(ipt.shape[0])
    for i in ex:
        tr += np.sum(ipt,axis=1) == i #This trace vector enable to to pick the desired vectors
    return np.repeat(ipt,tr>=1,axis=0)

def gsubs(n):
    """

    Built the subsummation matrix for a given number of input variables n

    Parameters
    ----------
    n : an integer
        number of input variables

    Returns
    -------
    subs : a numpy Boolean array of size 2**n x 2**n
        telling if the two literals are subsuming or not

    Examples
    --------
    >>> subsummation_mat = gsubs(2)

    """
    ipt = range(2**n)
    subs = []
    for c1 in ipt:
        line = []
        for c2 in ipt:
            if c1|c2==c1 or c1|c2==c2 or c1>c2: #Test for subsummation and order
                line.append(False) # Do not concatenate if they subsum
            else:
                line.append(True) # Concatenate their conjunction otherwise
        subs.append(line)
    return np.array(subs,dtype=np.bool)

def gperms(n):
    """

    Generate all the possible permutations of the input vectors you can obtain
    by permuting the label of the input lines.

    Parameters
    ----------
    n : an integer
        number of input variables

    Returns
    -------
    perms : a list of integer lists
        All possible permutations of the input vector set when labels are permuted.

    Examples
    --------
    >>> gperms(n=2)
    [[0,1,2,3],[0,2,1,3]]

    Comments
    --------
    The first list is the non-permuted input vector set.
    The number of lists is equal to n factorial, a vectorized form would be faster
    but it is okay for n=7.

    """
    ipt = [list(iptc) for iptc in iptgen(n,range(n+1))]
    permipt = [0 for i in range(len(ipt))]
    perms = list()

    for per in permutations(range(len(ipt[0]))):
        #Recontruct a new permuted input space
        permipt = [[ic[i] for i in per] for ic in ipt]
        #Recording this permuted input space into a list of list
        perms.append([int(''.join([str(i) for i in ic]),base=2) for ic in permipt])

    return perms

def nxtl(n, fprev, subs, perms):
    """

    Generate the set of positive representative Boolean functions of level mth
    from the one of level m-1th. A function is a tuple of integers which are the
    term (resp. clause) of the Disjunctive Normal Form (resp Conjunctive NF)


    Parameters
    ----------
    n : an integer
        number of input variables
    fprev : a list of tuples
        minimal positive DNFs (CNFs) from the previous level
    subsu : an Boolean numpy array
        prebuilt to determine if the two literals subsum
    perms : a list of list
        possible permutations


    Returns
    -------
    fnext : a set of tuples
        minimal positives DNF of the next level

    Examples
    --------
    >>> nxtl(3, [(i,) for i in range(2**3)], subs=gsubs(3), perms=gperms(3))
    """
    literals = [(i,) for i in range(2**n)] #Built the list of all possible literals
    fnext = set()
    for f1 in fprev: # Iterate though all DNF of level m-1
        for f2 in literals[max(f1):]: # select only the prime superior to the highest literal
            i = 0
            while subs[(f1[i],) + f2] != 0: #Test if any literal in the DNF subsum the prime
                i += 1 # go to the next literal in the DNF
                if i == len(f1): # if no literal in the DNF subsum the prime
                    fc = f1+f2 # concatenate the DNF with the prime
                    #Create of all possible children of fc
                    f = [tuple(sorted([pc[i] for i in fc])) for pc in perms]
                    #Add the child which is the "smallest" tuple
                    fnext.add(min(f))
                    break
    return fnext

def wrap(n, lmax, filename="record.hdf5"):
    """

    Wrap gsubs, gperms, and nxtl to generate all positive monotone function for a given number of variables

    Parameters
    ----------
    lmax : an integer
        level for which the generation stops
    filename : a char
        name of the hdf5 file

    Returns
    -------
    Nothing

    Examples
    --------
    >>> wrap(3, 2)

    """
    phi=set([(2**i-1,) for i in range(1,n+1)]) #Building level one
    perms = gperms(n)
    subsu = gsubs(n)

    gname =  "n" + str(n)
    dname = "/" + gname + "/level"
    print 2 # Counting the two level0 functions
    for i in range(1,lmax+1):
        hdf = h5py.File(filename,"a")
        print len(phi)
        rec = [tuple(fcur) for fcur in phi]
        rec.sort()
        if dname + str(i) not in hdf and phi:
            hdf.create_dataset(dname + str(i), data=np.array(rec), dtype=np.int8)
        hdf.close()
        phi = nxtl(n,phi,subsu,perms)

if __name__ == '__main__':
    wrap(3,3)
    wrap(4,6)
    wrap(5,10)
    #wrap(6,20)
