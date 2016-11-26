from numpy import *
from scipy.linalg import *
_author_ = 'Zhu Peihao'


class MatrixFactorization:

    'lu_factorization, gram_schmidt, householder_reduction, givens_reduction'
    matrix = array([[0, -20, -14], [3, 27, -4], [4, 11, -2]], float)

    def __init__(self, a=matrix.copy()):
        self.A = a.copy()

    def lu_factorization(self):
        """
            :param a: Input matrix
            :return: Decomposed matrices L, U, P
        """
        a = mat(self.A.copy())
        m, n = shape(a)  # Calculate the number of rows and columns of a matrix
        if m != n:
            print "\nThe matrix is not square"
        elif linalg.matrix_rank(a) != n:
            print "\nThis is not a nonsingular matrix "
        else:
            L = eye(n, dtype=float)  # Initialize matrix L, U, P
            U = zeros((n, n), float)
            P = zeros((n, n), float)
            b = arange(1, n + 1, dtype=float)  # b is used to record the order of exchange
            c = c_[a, b]
            for i in range(0, n - 1):  # Partial Pivoting
                temp = i
                for j in range(i + 1, n):  # Select and exchange the row that corresponds to  the max absolute value
                    if abs(c[j, i]) > abs(c[temp, i]):
                        temp = j
                temp_row = c[temp].copy()
                c[temp] = c[i]
                c[i] = temp_row
                for j in range(i + 1, n):  # Eliminate the following lines
                    if abs(c[j, i]) > 0:
                        temp_ratio = c[j, i] / c[i, i]
                        c[j, i: n] = c[j, i: n] - temp_ratio * c[i, i:n]
                        c[j, i] = temp_ratio
            for i in range(1, n):  # Calculate L
                for j in range(0, i):
                    L[i, j] = c[i, j]
            for i in range(0, n):  # Calculate U and P
                for j in range(i, n):
                    U[i, j] = c[i, j]
                P[i, int(c[i, n]) - 1] = 1
            return L, U, P

    def gram_schmidt(self):
        """
            :param a: Input matrix
            :return: Decomposed matrices Q, R
        """
        a = self.A.copy()
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
            return Q, R

    def householder_reduction(self):
        """
            :param a: Input matrix
            :return: Decomposed matrices Q, R
        """
        a = self.A.copy()
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
                # a = array(mat(R)*mat(a))
                a[i:m, i:n] = dot(R[i:m, i:m], a[i:m, i:n])
                Q = dot(Q, R.T)
            R_final = around(a, 4)  # Decimal point for accuracy
            return Q, R_final

    def givens_reduction(self):
        """
            :param a: Input matrix
            :return: Decomposed matrices Q, R
        """
        a = self.A.copy()
        m, n = shape(a)
        if linalg.matrix_rank(a) != n:
            print "\nThe rank is not n "
        else:
            A = a.copy()
            Q = eye(m, m)
            for i in range(n - 1):
                for j in range(i + 1, m):
                    P = eye(m, m)
                    c = a[i, i] / sqrt(power(a[i, i], 2) + power(a[j, i], 2))
                    s = a[j, i] / sqrt(power(a[i, i], 2) + power(a[j, i], 2))
                    P[i, i] = c
                    P[i, j] = s
                    P[j, i] = -s
                    P[j, j] = c
                    a = dot(P, a)
                    Q = dot(Q, P.T)
            R = around(a, 2)
            return Q, R


def create_matrix():
    """
    :return: Generate the default matrix A
    """
    A = array([[0, -20, -14], [3, 27, -4], [4, 11, -2]], float)
    return A


def main():
    A = create_matrix()
    a = MatrixFactorization(A)
    while True:
        flag1 = input('1:LU  2:gram_schmidt  3:householder_reduction  4:givens_reduction\n')
        if flag1 == 1:
            L, U, P = a.lu_factorization()
            print 'A = ', '\n', A, '\n', 'L = ', '\n', L, '\n', 'U = ', '\n', U, '\n', 'P = ', '\n', P, '\n'
        elif flag1 == 2:
            Q, R = a.gram_schmidt()
            print 'A = ', '\n', A, '\n', 'Q = ', '\n', Q, '\n','R = ', '\n', R, '\n'
        elif flag1 == 3:
            Q, R = a.householder_reduction()
            print 'A = ', '\n', A, '\n', 'Q = ', '\n', Q, '\n','R = ', '\n', R, '\n'
        elif flag1 == 4:
            Q, R = a.givens_reduction()
            print 'A = ', '\n', A, '\n', 'Q = ', '\n', Q, '\n','R = ', '\n', R, '\n'
        else:
            print ' no correct input'

        flag2 = input('continue: 1, end: 2\n')
        if flag2 == 2:
            break

if '__main__' == __name__:
    main()