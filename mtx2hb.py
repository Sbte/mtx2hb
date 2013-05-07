''' Small script to convert Matrix Market format matrices to
    Harwell-Boeing format matrices

    This script takes as first command-line argument the path to the
    matrix and the next arguments can be multiple files containing
    right-hand sides.

    Sven Baars, 2013.'''

import sys
import math
import os

class MatrixType:
    General = 1
    Symmetric = 2
    Array = 3

class Matrix(list):
    def __init__(self, m, n, type_, name=''):
        list.__init__(self)
        self.m = m
        self.n = n
        self.type = type_
        self.name = name

    nnz = property(len)

def read_mtx_matrix(name):
    matrix_in_file = open(name, 'r')
    lines = matrix_in_file.readlines()
    matrix_in_file.close()

    type_ = 0
    if lines[0] == '%%MatrixMarket matrix coordinate real general\n':
        type_ = MatrixType.General
    elif lines[0] == '%%MatrixMarket matrix coordinate real symmetric\n':
        type_ = MatrixType.Symmetric
    elif lines[0] == '%%MatrixMarket matrix array real general\n':
        type_ = MatrixType.Array
    else:
        print 'Invalid format for', name
        print lines[0]
        sys.exit(0)

    # Read mtx format
    if type_ == MatrixType.Array:
        first = True
        for line in lines:
            if line.startswith('%'):
                continue
            if first:
                m, n = [int(i) for i in line.split()]
                matrix = Matrix(m, n, type_, name)
                first = False
                continue
            v = line.strip()
            matrix.append([len(matrix), 1, float(v)])
    else:
        first = True
        for line in lines:
            if line.startswith('%'):
                continue
            if first:
                m, n, nnz = [int(i) for i in line.split()]
                matrix = Matrix(m, n, type_, name)
                first = False
                continue
            i, j, v = line.split()
            matrix.append([int(i), int(j), float(v)])
            if type_ == MatrixType.Symmetric and i != j:
                matrix.append([int(j), int(i), float(v)])

    return matrix

def coocsr(matrix):
    n = matrix.m
    nnz = matrix.nnz

    nnzrow = [0]*n
    for i in range(nnz):
        nnzrow[matrix[i][0]-1] += 1

    ia = [1]*(n+1)
    for i in range(n):
        ia[i+1] = ia[i] + nnzrow[i]

    ja = [0]*nnz
    a  = [0]*nnz
    ia2 = list(ia)
    for i in range(nnz):
        i2 = ia2[matrix[i][0]-1]-1
        ja[i2] = matrix[i][1]
        a[i2] = matrix[i][2]
        ia2[matrix[i][0]-1] += 1

    return (ia, ja, a)

def coocsc(matrix):
    n = matrix.n
    nnz = matrix.nnz

    nnzcol = [0]*n
    for i in range(nnz):
        nnzcol[matrix[i][1]-1] += 1

    ja = [1]*(n+1)
    for i in range(n):
        ja[i+1] = ja[i] + nnzcol[i]

    ia = [0]*nnz
    a  = [0]*nnz
    ja2 = list(ja)
    for i in range(nnz):
        i2 = ja2[matrix[i][1]-1]-1
        ia[i2] = matrix[i][0]
        a[i2] = matrix[i][2]
        ja2[matrix[i][1]-1] += 1

    return (ia, ja, a)

def write_hb_matrix(matrix, rhs=None):
    if rhs is None:
        rhs = []

    ia, ja, a = coocsc(matrix)
    rhsa = [j for i in rhs for j in coocsc(i)[2]]

    m = matrix.m
    n = matrix.n
    nnz = matrix.nnz

    # Default format
    title = '1U'
    ifmt = 8

    if matrix.type == MatrixType.Symmetric:
        type_ = 'RSA'
    else:
        type_ = 'RUA'

    l = int(math.ceil(math.log10(0.1 + nnz + 1) ) + 1)

    # Compute column pointer format
    ptr_len = l
    ptr_nperline = min(80 // ptr_len, n+1)
    ptrcrd = n // ptr_nperline + 1
    ptrfmt = '(%dI%d)' % (ptr_nperline, l)

    # Compute row index format
    ind_len = l
    ind_nperline = min(80 // ind_len, nnz)
    indcrd = (nnz-1) // ind_nperline + 1
    indfmt = '(%dI%d)' % (ind_nperline, l)

    # Compute values and rhs format (same format for both)
    valcrd = 0
    rhscrd = 0
    if ifmt >= 100:
        ihead = ifmt // 100
        ifmt = ifmt - 100 * ihead
        l = ihead + ifmt + 2 #?????
        val_nperline = 80 // l
        c_len = l # 80 / nperline
        c_valfmt = '%%%d.%df' % (c_len, ifmt)
        valfmt = '%dF%d.%d' % (val_nperline, l, ifmt)
    else:
        l = ifmt + 8
        val_nperline = 80 // l
        c_len = l # 80 / nperline
        c_valfmt = '%%%d.%dE' % (c_len, ifmt)
        valfmt = '%dE%d.%d' % (val_nperline, l, ifmt)

    valcrd = (nnz-1) // val_nperline + 1
    valfmt = '('+valfmt+')'

    nrhs = len(rhs)
    if nrhs >= 1:
        rhscrd = (nrhs*m - 1) // val_nperline + 1

    totcrd = ptrcrd + indcrd + valcrd + rhscrd

    # Now write 4-line header
    outfile = open(matrix.name.replace('.mtx', '.hb'), 'w')

    # Line 1
    t = title
    t += ' ' * (72-len(title))

    key = os.path.splitext(os.path.split(matrix.name)[-1])[0]
    if len(key) > 8:
        t += key[:8]
    else:
        t += key
        t += ' ' * (8-len(key))

    t += '\n'
    outfile.write(t)

    # Line 2
    outfile.write('%14d%14d%14d%14d%14d%10s\n' % (totcrd, ptrcrd, indcrd, valcrd, rhscrd, ''))

    # Line 3
    t = type_
    t += ' ' * (14-len(type_))
    outfile.write('%14s%14i%14i%14i%14i%10s\n' % (t, m, n, nnz, nrhs, ''))

    # Line 4
    outfile.write('%16s%16s%20s%20s%8s\n' % (ptrfmt, indfmt, valfmt, valfmt, ''))

    # column pointers
    t = (('%%%dd' % ptr_len) * ptr_nperline + '\n') * ((n+1) // ptr_nperline)
    k = (n+1) % ptr_nperline
    if k > 0:
        t += ('%%%dd' % ptr_len) * k + '\n'
    outfile.write(t % tuple(ja))

    # row indices
    t = (('%%%dd' % ind_len) * ind_nperline + '\n') * (nnz // ind_nperline)
    k = nnz % ind_nperline
    if k > 0:
        t += ('%%%dd' % ptr_len) * k + '\n'
    outfile.write(t % tuple(ia))

    # numerical values of nonzero elements of the matrix
    t = (c_valfmt * val_nperline + '\n') * (nnz // val_nperline)
    k = nnz % val_nperline
    if k > 0:
        t += c_valfmt * k + '\n'
    outfile.write(t % tuple(a))

    # numerical values of right hand side
    t = (c_valfmt * val_nperline + '\n') * ((m * nrhs) // val_nperline)
    k = (m * nrhs) % val_nperline
    if k > 0:
        t += c_valfmt * k + '\n'
    outfile.write(t % tuple(rhsa))

    outfile.close()

def main():
    if len(sys.argv) < 2:
        print 'No matrix passed'
        sys.exit(0)

    name = sys.argv[1]
    matrix = read_mtx_matrix(name)
    rhs = []
    if len(sys.argv) > 2:
        for i in sys.argv[2:]:
            rhs.append(read_mtx_matrix(i))
    write_hb_matrix(matrix, rhs)

if __name__ == '__main__':
    main()
