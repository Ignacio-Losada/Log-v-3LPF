import logv3lpf.linpf as linpf
import opendssdirect as dss
from logv3lpf.DSSParser import DSScase
import numpy as np
from logv3lpf.math_functions import polar_to_cartesian
import datetime

class case(object):
    def __init__(self, file,sourcebus,refvm,refva):
        DSScase.__init__(self,file,sourcebus,refvm,refva)
        DSScase.parse_DSS_file(self)
        print("Calculating Line Matrices")
        linpf.get_line_matrices(self)
        print("Calculating Transformer Matrices")
        linpf.get_transformer_matrices(self)
        print("Calculating Regulator Matrices")
        linpf.get_regulator_matrices(self)
        print("Calculating Base Matrices")
        linpf.calculate_base_matrices(self)
        case.sourcebus = sourcebus
        linpf.check_logv3lpf_performance(self)
    
    def run_case(self,algorithm):
        """
        Compute the eigenvalues and right eigenvectors of a square array.
        Parameters
        ----------
        a : (..., M, M) array
            Matrices for which the eigenvalues and right eigenvectors will
            be computed
        Returns
        -------
        w : (..., M) array
            The eigenvalues, each repeated according to its multiplicity.
            The eigenvalues are not necessarily ordered. The resulting
            array will be of complex type, unless the imaginary part is
            zero in which case it will be cast to a real type. When `a`
            is real the resulting eigenvalues will be real (0 imaginary
            part) or occur in conjugate pairs
        v : (..., M, M) array
            The normalized (unit "length") eigenvectors, such that the
            column ``v[:,i]`` is the eigenvector corresponding to the
            eigenvalue ``w[i]``.
        Raises
        ------
        LinAlgError
            If the eigenvalue computation does not converge.
        See Also
        --------
        eigvals : eigenvalues of a non-symmetric array.
        eigh : eigenvalues and eigenvectors of a real symmetric or complex
            Hermitian (conjugate symmetric) array.
        eigvalsh : eigenvalues of a real symmetric or complex Hermitian
                (conjugate symmetric) array.
        scipy.linalg.eig : Similar function in SciPy that also solves the
                        generalized eigenvalue problem.
        scipy.linalg.schur : Best choice for unitary and other non-Hermitian
                            normal matrices.
        Notes
        -----
        .. versionadded:: 1.8.0
        Broadcasting rules apply, see the `numpy.linalg` documentation for
        details.
        This is implemented using the ``_geev`` LAPACK routines which compute
        the eigenvalues and eigenvectors of general square arrays.
        The number `w` is an eigenvalue of `a` if there exists a vector
        `v` such that ``a @ v = w * v``. Thus, the arrays `a`, `w`, and
        `v` satisfy the equations ``a @ v[:,i] = w[i] * v[:,i]``
        for :math:`i \\in \\{0,...,M-1\\}`.
        The array `v` of eigenvectors may not be of maximum rank, that is, some
        of the columns may be linearly dependent, although round-off error may
        obscure that fact. If the eigenvalues are all different, then theoretically
        the eigenvectors are linearly independent and `a` can be diagonalized by
        a similarity transformation using `v`, i.e, ``inv(v) @ a @ v`` is diagonal.
        For non-Hermitian normal matrices the SciPy function `scipy.linalg.schur`
        is preferred because the matrix `v` is guaranteed to be unitary, which is
        not the case when using `eig`. The Schur factorization produces an
        upper triangular matrix rather than a diagonal matrix, but for normal
        matrices only the diagonal of the upper triangular matrix is needed, the
        rest is roundoff error.
        Finally, it is emphasized that `v` consists of the *right* (as in
        right-hand side) eigenvectors of `a`.  A vector `y` satisfying
        ``y.T @ a = z * y.T`` for some number `z` is called a *left*
        eigenvector of `a`, and, in general, the left and right eigenvectors
        of a matrix are not necessarily the (perhaps conjugate) transposes
        of each other.
        References
        ----------
        G. Strang, *Linear Algebra and Its Applications*, 2nd Ed., Orlando, FL,
        Academic Press, Inc., 1980, Various pp.
        Examples
        --------
        >>> from numpy import linalg as LA
        (Almost) trivial example with real e-values and e-vectors.
        >>> w, v = LA.eig(np.diag((1, 2, 3)))
        >>> w; v
        array([1., 2., 3.])
        array([[1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.]])
        Real matrix possessing complex e-values and e-vectors; note that the
        e-values are complex conjugates of each other.
        >>> w, v = LA.eig(np.array([[1, -1], [1, 1]]))
        >>> w; v
        array([1.+1.j, 1.-1.j])
        array([[0.70710678+0.j        , 0.70710678-0.j        ],
            [0.        -0.70710678j, 0.        +0.70710678j]])
        Complex-valued matrix with real e-values (but complex-valued e-vectors);
        note that ``a.conj().T == a``, i.e., `a` is Hermitian.
        >>> a = np.array([[1, 1j], [-1j, 1]])
        >>> w, v = LA.eig(a)
        >>> w; v
        array([2.+0.j, 0.+0.j])
        array([[ 0.        +0.70710678j,  0.70710678+0.j        ], # may vary
            [ 0.70710678+0.j        , -0.        +0.70710678j]])
        Be careful about round-off error!
        >>> a = np.array([[1 + 1e-9, 0], [0, 1 - 1e-9]])
        >>> # Theor. e-values are 1 +/- 1e-9
        >>> w, v = LA.eig(a)
        >>> w; v
        array([1., 1.])
        array([[1., 0.],
            [0., 1.]])
        """
        if algorithm=="logv3lpf":
            print("Solving base case logv3lpf...")
            case.run_logv3lpf(self)
        elif algorithm=="NR":
            case.run_newton_raphson(self)
        elif algorithm=="GS":
            case.run_gauss_seidel(self)
        elif algorithm=="FBS":
            case.run_forward_backward_sweep(self)
        elif algorithm=="openDSS":
            print("Solving DSS case in OpenDSS...")
            dss.Solution.Solve()
            print("Storing solution from OpenDSS...")
            DSScase.process_openDSS_solution(self)
        else:
            raise Error('Algorithm not implemented. Options are logv3lpf, NR, GS, FBS. See documentation for details')
    
    def run_logv3lpf(self):
        linpf.rank_k_correction_solve(self,True)
    # def update_logv3lpf(self):
    #     linpf.update_logv3lpf(self,True)

    def run_newton_raphson(self):
        return False
    def run_gauss_seidel(self):
        return False
    def run_forward_backward_sweep(self):
        return False

    def set_load_kvar(self,name,value):
        try:
            self.loads.loc[self.loads.name==name,"kvar"] = value
        except ValueError:
            print("Load {:} does not exist!".format(name))

    def set_load_kW(self,name,value):
        try:
            self.loads.loc[self.loads.name==name,"kW"] = value
        except ValueError:
            print("Load {:} does not exist!".format(name))

    def set_capacitor_kvar(self,name,value):
        try:
            self.capacitors.loc[self.capacitors.name==name,"kvar"] = value
        except ValueError:
            print("Capacitor {:} does not exist!".format(name))

    def set_transformer_tap(self,name,value,bus):
        try:
            self.transformers.loc[self.transformers.name==name, bus] = value
            if name in list(self.regcontrols.name):
                self.regcontrols.loc[self.regcontrols.name==name,"tap"] = value
        except ValueError:
            print("Transformer {:} does not exist or {:} is not a valid bus!".format(name.bus))

    def get_line_current(self,name,pu):
        idx = self.lines[self.lines.name == name].index[0]
        fbus,tbus = self.lines.fbus[idx],self.lines.tbus[idx]
        fphase,tphase = self.lines.fphase[idx],self.lines.tphase[idx]
        fidx = [self.bus_phases[fbus].index(i) for i in fphase]
        tidx = [self.bus_phases[tbus].index(i) for i in tphase]
        un = np.log(np.array([self.results["logv3lpf"]["vm"][fbus][i] for i in fidx]))
        um = np.log(np.array([self.results["logv3lpf"]["vm"][tbus][i] for i in tidx]))
        θn = np.log(np.array([self.results["logv3lpf"]["va"][fbus][i]-case.base[fbus][i] for i in fidx]))
        θm = np.log(np.array([self.results["logv3lpf"]["va"][tbus][i]-case.base[tbus][i] for i in tidx]))
        u = np.hstack([un,um]).reshape(-1,1)
        θ = np.hstack([θn,θm]).reshape(-1,1)
        y = case.yl[case.lines["nphases"][:idx].sum()*2:case.lines["nphases"][:idx].sum()*2+case.lines["nphases"][idx]*2]

        s = case.Yul[idx]@u + case.Yθl[idx]@θ + y

        if pu:
            return s
        else:
            scaling = np.diag([self.base[fbus]["IBase"]]*len(fidx)+[self.base[tbus]["IBase"]]*len(tidx))
            return scaling@I


        


