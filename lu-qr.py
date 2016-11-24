from numpy import *
from scipy.linalg import *
_author_ = 'Zhu Peihao'


def create_matrix():
    """
    :return: Generate the default matrix A
    """
    A = array([[0, -20, -14], [3, 27, -4], [4, 11, -2]], float)
    return A


def gram_schmidt(a):
    m, n = shape(a)
    if linalg.matrix_rank(a) != n:
        print "\nThe rank is not n "
    else:
        R = zeros((m, n), float)
        Q = zeros((m, m), float)
        for i in range(0, n):
            Q[:, i] = a[:, i]
            for j in range(0, i):
                R[j, i] = dot(Q[:, j], a[:, i])
                Q[:, i] -= R[j, i] * Q[:, j]
            R[i, i] = norm(Q[:, i])
            Q[:, i] = Q[:, i] / R[i, i]
        print "A =", '\n', a
        print "Q = ", '\n', Q
        print "R = ", '\n', R
        return Q, R


def householder_reduction(a):
    m, n = shape(a)
    if linalg.matrix_rank(a) != n:
        print "\nThe rank is not n "
    else:
        A = a.copy()
        Q = eye(m, m)
        for i in range(n - 1):
            R = eye(m, m)
            u = a[i: m, i].copy()
            u[0] -= norm(a[i:, i])
            R[i:m, i:m] = eye(m - i, m - i) - 2 * array(mat(u).T * mat(u)) / dot(u, u)

            print 'R',R
            print 'a',a
            print 'temp',u
            #a = array(mat(R)*mat(a))
            a[i:m, i:n] = dot(R[i:m, i:m], a[i:m, i:n])
            print 'change a',a
            Q = dot(Q, R.T)
        R_final = around(a,4)        # Decimal point for accuracy
        print "A =", '\n', A
        print "Q = ", '\n', Q
        print "R = ", '\n', R_final
        return Q, R_final


def lu_factorization(a):
    """
    :param a: Input matrix
    :return: Decomposed matrices L, U, P
    """
    a = mat(a)
    m, n = shape(a)  # Calculate the number of rows and columns of a matrix
    if m != n:
        print "\nThe matrix is not square"
    elif linalg.matrix_rank(a) != n:
        print "\nThis is not a nonsingular matrix "
    else:
        L = eye(n, dtype=float)     # Initialize matrix L, U, P
        U = zeros((n, n), float)
        P = zeros((n, n), float)
        b = arange(1, n + 1, dtype=float)       # b is used to record the order of exchange
        c = c_[a, b]
        for i in range(0, n - 1):       # Partial Pivoting
            temp = i
            for j in range(i + 1, n):       # Select and exchange the row that corresponds to  the max absolute value
                if abs(c[j, i]) > abs(c[temp, i]):
                    temp = j
            temp_row = c[temp].copy()
            c[temp] = c[i]
            c[i] = temp_row
            for j in range(i + 1, n):       # Eliminate the following lines
                if abs(c[j, i]) > 0:
                    temp_ratio = c[j, i] / c[i, i]
                    c[j, i: n] = c[j, i: n] - temp_ratio * c[i, i:n]
                    c[j, i] = temp_ratio
        for i in range(1, n):       # Calculate L
            for j in range(0, i):
                L[i, j] = c[i, j]
        for i in range(0, n):       # Calculate U and P
            for j in range(i, n):
                U[i, j] = c[i, j]
            P[i, int(c[i, n]) - 1] = 1
        print "A =", '\n', a
        print "L = ", '\n', L
        print "U = ", '\n', U
        print "P = ", '\n', P
        return L, U, P


if '__main__' == __name__:
    A = create_matrix()
    #lu_factorization(A)
    #gram_schmidt(A)
    householder_reduction(A)