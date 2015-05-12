#
#
#  multired.py 
#
# 
#  Copyright (C) 2015 Vincenzo (Enzo) Nicosia <katolaz@yahoo.it>
# 
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.  
# 
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  long with this program.  If not, see  <http://www.gnu.org/licenses/>.
# 
# 
#
#  This module provides the class multiplex_red, which implements the
#  algorithm for the structural reduction of multi-layer networks
#  based on the Von Neumann entropy and Quantum Jensen-Shannon
#  divergence of graphs. 
#
#  If you use this code please cite the paper:
#
#    M. De Domenico, V. Nicosia, A. Arenas, V. Latora, 
#    "Structural reducibility of multilayer networks" 
#    Nat. Commun. 6, 6864 (2015) doi:10.1038/ncomms7864
#
# --------------------------------------------
#  
#  -- 2015/04/23 -- release 0.1 
#  -- 2015/05/11 -- release 0.1.1 -- removed the last full matrices 
# 



import sys
import math 
import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.linalg import eigh, eig
import copy
from scipy.cluster.hierarchy import linkage, dendrogram

has_matplotlib = False

try:
    import matplotlib
    has_matplotlib = True

except ImportError:
    has_matplotlib = False


class XLogx_fit:
    def __init__(self, degree, npoints= 100, xmax=1):
        if xmax > 1:
            xmax = 1
        self.degree = degree
        x = np.linspace(0, xmax, npoints)
        y = [i * math.log(i) for i in x[1:]]
        y.insert(0, 0)
        self.fit = np.polyfit(x, y, degree)

    def __getitem__ (self, index):
        if index <= self.degree:
            return self.fit[index]
        else:
            print "Error!!! Index %d is larger than the degree of the fitting polynomial (%d)" \
                % (index, degree)
            sys.exit(-1)


class layer:
    def __init__ (self, layerfile= None, matrix=None):
        self.N = 0
        self.num_layer = -1
        self.fname = layerfile
        self.adj_matr = None
        self.laplacian = None
        self.resc_laplacian = None
        self.entropy = None
        self.entropy_approx = None
        self._ii = []
        self._jj = []
        self._ww = []
        self._matrix_called = False
        if layerfile != None:
            try:
                min_N = 10e10
                with open(layerfile, "r") as lines:
                    for l in lines:
                        if l[0] == '#':
                            continue
                        elems = l.strip(" \n").split(" ")
                        s = int(elems[0])
                        d = int(elems[1])
                        self._ii.append(s)
                        self._jj.append(d)
                        #if s > self.N:
                        #    self.N = s
                        #if d > self.N:
                        #    self.N = d
                        if s < min_N:
                            min_N = s
                        if  d < min_N:
                            min_N = d
                        if len(elems) >2 : ## A weight is specified 
                            val = [float(x) if "e" in x or "." in x else int(x) for x in [elems[2]]][0]
                            self._ww.append(float(val))
                        else:
                            self._ww.append(int(1))
                m1 = max(self._ii)
                m2 = max(self._jj)
                if m1 > m2:
                    self.N = m1
                else:
                    self.N = m2
                
            except (IOError):
                print "Unable to find/open file %s -- Exiting!!!" % layerfile
                exit(-2)
        elif matrix != None:
            self.adj_matr = copy.copy(matrix)
            self.N, _x = matrix.shape 
            K = self.adj_matr.sum(0).reshape((1, self.N)).tolist()[0]
            D = csr_matrix((K, (range(self.N), range(self.N)) ), shape=(self.N, self.N))
            self.laplacian = csr_matrix(D - self.adj_matr)
            K = self.laplacian.diagonal().sum()
            self.resc_laplacian = csr_matrix(self.laplacian / K)
            self._matrix_called = True
        else:
            print "The given matrix is BLANK"
    def make_matrices(self, N):
        self.N = N 
        self.adj_matr = csr_matrix((self._ww, (self._ii, self._jj)), shape=(self.N, self.N))
        self.adj_matr = self.adj_matr + self.adj_matr.transpose()
        K = self.adj_matr.sum(0).reshape((1, self.N)).tolist()[0]
        D = csr_matrix((K, (range(self.N), range(self.N)) ), shape=(self.N, self.N))
        self.laplacian = csr_matrix(D - self.adj_matr)
        K = self.laplacian.diagonal().sum()
        self.resc_laplacian = csr_matrix(self.laplacian / K)
        self._matrix_called = True

    def dump_info(self):
        N, M = self.adj_matr.shape
        K = self.adj_matr.nnz
        sys.stderr.write("Layer File: %s\nNodes: %d Edges: %d\nEntropy: %g Approx. Entropy: %g\n" % \
                             (self.fname, N, K, self.entropy, self.entropy_approx) )

    def compute_VN_entropy(self):
        eigvals = eigh(self.resc_laplacian.todense())

        self.entropy = 0
        for l_i in eigvals[0]:
            if (l_i > 10e-20):
                self.entropy -= l_i * math.log (l_i)


    def compute_VN_entropy_approx(self, poly):
        p = poly.degree
        h = - poly[p] * self.N
        M = csr_matrix(eye(self.N))
        for i in range(p-1, -1, -1):
            M = M *  self.resc_laplacian
            h += - poly[i] * sum(M.diagonal())
        self.entropy_approx = h

    def aggregate(self, other_layer):
        if self.adj_matr != None:
            self.adj_matr = self.adj_matr + other_layer.adj_matr
        else:
            self.adj_matr = copy.copy(other_layer.adj_matr)
        K = self.adj_matr.sum(0).reshape((1, self.N)).tolist()[0]
        D = csr_matrix((K, (range(self.N), range(self.N)) ), shape=(self.N, self.N))
        self.laplacian = csr_matrix(D - self.adj_matr)
        K = self.laplacian.diagonal().sum()
        self.resc_laplacian = csr_matrix(self.laplacian / K)
        self._matrix_called = True

    def dump_laplacian(self):
        print self.laplacian

class multiplex_red:
    
    def __init__ (self, multiplexfile, directed = None, fit_degree=10, verbose=False):
        self.layers = []
        self.N = 0
        self.M = 0
        self.entropy = 0
        self.entropy_approx = 0
        self.JSD = None
        self.JSD_approx = None
        self.Z = None
        self.Z_approx = None
        self.aggr = None
        self.q_vals = None
        self.q_vals_approx = None
        self.fit_degree = fit_degree
        self.poly = XLogx_fit(self.fit_degree)
        self.verb = verbose
        self.cuts = None
        self.cuts_approx = None
        try:
            with open(multiplexfile, "r") as lines:
                for l in lines:
                    if (self.verb):
                        sys.stderr.write("Loading layer %d from file %s" % (len(self.layers), l))
                    A = layer(l.strip(" \n"))
                    #if A.N > self.N:
                    #    self.N = A.N+1
                    self.layers.append(A)
            n = 0
            N = max([x.N for x in self.layers])
            self.N = N + 1
            for l in self.layers:
                l.make_matrices(self.N)
                l.num_layer = n
                n += 1
            self.M = len(self.layers)
        except ( IOError):
            print "Unable to find/open file %s -- Exiting!!!" % layer_file
            exit(-2)

    def dump_info(self):
        i = 0
        for l in self.layers:
            sys.stderr.write("--------\nLayer: %d\n" % i)
            l.dump_info()
            i += 1


    def compute_aggregated(self):
        self.aggr = copy.copy(self.layers[0])
        self.aggr.entropy = 0
        self.aggr.entropy_approx = 0
        for l in self.layers[1:]:
            self.aggr.aggregate(l)

    def compute_layer_entropies(self):
        for l in self.layers:
            l.compute_VN_entropy()

    def compute_layer_entropies_approx(self):
        for l in self.layers:
            l.compute_VN_entropy_approx(self.poly)


    def compute_multiplex_entropy(self, force_compute=False):
        ### The entropy of a multiplex is defined as the sum of the entropies of its layers
        for l in self.layers:
            if l.entropy == None:
                l.compute_VN_entropy()
                self.entropy += l.entropy

    def compute_multiplex_entropy_approx(self, force_compute=False):
        ### The entropy of a multiplex is defined as the sum of the entropies of its layers
        for l in self.layers:
            if l.entropy_approx == None:
                l.compute_VN_entropy_approx(self.poly)
            self.entropy_approx += l.entropy_approx

    def compute_JSD_matrix(self):
        if (self.verb):
            sys.stderr.write("Computing JSD matrix\n")
        self.JSD = np.zeros((self.M, self.M))
        for i in range(len(self.layers)):
            for j in range(i+1, len(self.layers)):
                li = self.layers[i]
                lj = self.layers[j]
                if not li.entropy:
                    li.compute_VN_entropy()
                if not lj.entropy:
                    lj.compute_VN_entropy()
                # m_sigma = (li.resc_laplacian + lj.resc_laplacian)/2.0
                # m_sigma_entropy = mr.compute_VN_entropy_LR(m_sigma)
                m_sigma_matr = (li.adj_matr + lj.adj_matr)/2.0
                m_sigma = layer(matrix=m_sigma_matr)
                m_sigma.compute_VN_entropy()
                d = m_sigma.entropy - 0.5 * (li.entropy + lj.entropy)
                d = math.sqrt(d)
                self.JSD[i][j] = d
                self.JSD[j][i] = d
        pass

    def compute_JSD_matrix_approx(self):
        if (self.verb):
            sys.stderr.write("Computing JSD matrix (approx)\n")
        self.JSD_approx = np.zeros((self.M, self.M))
        for i in range(len(self.layers)):
            for j in range(i+1, len(self.layers)):
                li = self.layers[i]
                lj = self.layers[j]
                if not li.entropy_approx:
                    li.compute_VN_entropy_approx(self.poly)
                if not lj.entropy_approx:
                    lj.compute_VN_entropy_approx(self.poly)
                m_sigma_matr = (li.adj_matr + lj.adj_matr)/2.0
                m_sigma = layer(matrix=m_sigma_matr)
                m_sigma.compute_VN_entropy_approx(self.poly)
                d = m_sigma.entropy_approx - 0.5 * (li.entropy_approx + lj.entropy_approx)
                d = math.sqrt(d)
                self.JSD_approx[i][j] = d
                self.JSD_approx[j][i] = d
                
    def dump_JSD(self, force_compute=False):
        if self.JSD == None:
            if force_compute:
                self.compute_JSD_matrix()
            else:
                print "Error!!! call to dump_JSD but JSD matrix has not been computed!!!"
                sys.exit(1)
        idx = 0
        for i in range(self.len):
            for j in range(i+1, self.len):
                print i, j, self.JSD[idx]
                idx += 1

    def dump_JSD_approx(self, force_compute=False):
        if self.JSD_approx == None:
            if force_compute:
                self.compute_JSD_matrix_approx()
            else:
                print "Error!!! call to dump_JSD_approx but JSD approximate matrix has not been computed!!!"
                sys.exit(1)
        idx = 0
        for i in range(self.M):
            for j in range(i+1, self.M):
                print i, j, self.JSD_approx[idx]
                idx += 1
    

    def reduce(self, method="ward"):
        if (self.verb):
            sys.stderr.write("Performing '%s' reduction\n" % method)
        if self.JSD == None:
            self.compute_JSD_matrix()
        self.Z = linkage(self.JSD, method=method)
        return self.Z

    def reduce_approx(self, method="ward"):
        if (self.verb):
            sys.stderr.write("Performing '%s' reduction (approx)\n" % method)
        if self.JSD_approx == None:
            self.compute_JSD_matrix_approx()
        self.Z_approx = linkage(self.JSD_approx, method=method)
        return self.Z_approx
    
    def get_linkage(self):
        return self.Z

    def get_linkage_approx(self):
        return self.Z_approx

    def __compute_q(self, layers):
        H_avg = 0
        if not self.aggr:
            self.compute_aggregated()
            self.aggr.compute_VN_entropy()
        for l in layers:
            if not l.entropy:
                l.compute_VN_entropy()
            H_avg += l.entropy
        H_avg /= len(layers)
        q = 1.0 - H_avg / self.aggr.entropy
        return q

    def get_q_profile(self):
        mylayers = copy.copy(self.layers)
        rem_layers = copy.copy(self.layers)
        q_vals = []
        if self.Z == None:
            self.reduce()
        q = self.__compute_q(rem_layers)
        q_vals.append(q)
        n = len(self.layers)
        for l1, l2, _d, _x in self.Z:
            l_new = layer(matrix=mylayers[int(l1)].adj_matr)
            l_new.num_layer = n 
            n += 1 
            l_new.aggregate(mylayers[int(l2)])
            rem_layers.remove(mylayers[int(l1)])
            rem_layers.remove(mylayers[int(l2)])
            rem_layers.append(l_new)
            mylayers.append(l_new)
            q = self.__compute_q(rem_layers)
            q_vals.append(q)
        self.q_vals = q_vals
        return q_vals
        pass

    
    def __compute_q_approx(self, layers):
        H_avg = 0
        if not self.aggr:
            self.compute_aggregated()
            self.aggr.compute_VN_entropy_approx(self.poly)
        for l in layers:
            if not l.entropy_approx:
                l.compute_VN_entropy_approx(self.poly)
            H_avg += l.entropy_approx
        H_avg /= len(layers)
        q = 1.0 - H_avg / self.aggr.entropy_approx
        return q

    def get_q_profile_approx(self):
        mylayers = copy.copy(self.layers)
        rem_layers = copy.copy(self.layers)
        q_vals = []
        if self.Z_approx == None:
            self.reduce_approx()
        q = self.__compute_q_approx(rem_layers)
        q_vals.append(q)
        n = len(self.layers)
        for l1, l2, _d, _x in self.Z_approx:
            l_new = layer(matrix=mylayers[int(l1)].adj_matr)
            l_new.num_layer = n 
            n += 1 
            l_new.aggregate(mylayers[int(l2)])
            rem_layers.remove(mylayers[int(l1)])
            rem_layers.remove(mylayers[int(l2)])
            rem_layers.append(l_new)
            mylayers.append(l_new)
            q = self.__compute_q_approx(rem_layers)
            q_vals.append(q)
        self.q_vals_approx = q_vals
        return q_vals
        
    def compute_partitions(self):
        if (self.verb):
            sys.stderr.write("Getting partitions...\n")
        if self.Z == None:
            self.reduce()
        if self.q_vals == None:
            self.get_q_profile()
        sets = {}
        M = len(self.layers)
        for i in range(len(self.layers)):
            sets[i] = [i]
        best_pos = self.q_vals.index(max(self.q_vals))
        j = 0
        cur_part = sets.values()
        self.cuts = [copy.deepcopy(cur_part)]
        while j < M-1:
            l1, l2, _x, _y = self.Z[j]
            l1 = int(l1)
            l2 = int(l2)
            val = sets[l1]
            val.extend(sets[l2])
            sets[M+j] = val
            r1 = cur_part.index(sets[l1])
            cur_part.pop(r1)
            r2 = cur_part.index(sets[l2])
            cur_part.pop(r2)
            cur_part.append(val)
            j += 1
            self.cuts.append(copy.deepcopy(cur_part))
        self.cuts.append(copy.deepcopy(cur_part))
        return zip(self.q_vals, self.cuts)


    def compute_partitions_approx(self):
        if (self.verb):
            sys.stderr.write("Getting partitions (approx)...\n")
        if self.Z_approx == None:
            self.reduce_approx()
        if self.q_vals_approx == None:
            self.get_q_profile_approx()
        sets = {}
        M = len(self.layers)
        for i in range(len(self.layers)):
            sets[i] = [i]
        best_pos = self.q_vals_approx.index(max(self.q_vals_approx))
        j = 0
        cur_part = sets.values()
        self.cuts_approx = [copy.deepcopy(cur_part)]
        while j < M-1:
            l1, l2, _x, _y = self.Z_approx[j]
            l1 = int(l1)
            l2 = int(l2)
            val = sets[l1]
            val.extend(sets[l2])
            sets[M+j] = val
            r1 = cur_part.index(sets[l1])
            cur_part.pop(r1)
            r2 = cur_part.index(sets[l2])
            cur_part.pop(r2)
            cur_part.append(val)
            j += 1
            self.cuts_approx.append(copy.deepcopy(cur_part))
        self.cuts_approx.append(copy.deepcopy(cur_part))
        return zip(self.q_vals_approx, self.cuts_approx)

    def draw_dendrogram(self, force = False):
        if not has_matplotlib:
            sys.stderr.write("No matplotlib module found in draw_dendrogram...Exiting!!!\n")
            sys.exit(3)
        if self.Z == None:
            if not force:
                sys.stderr.write("Please call reduce() first or specify 'force=True'")
            else:
                self.reduce()
        dendrogram(self.Z, no_plot=False)
        matplotlib.pyplot.draw()
        matplotlib.pyplot.show()

    def draw_dendrogram_approx(self, force = False):
        if not has_matplotlib:
            sys.stderr.write("No matplotlib module found in draw_dendrogram_approx...Exiting!!!\n")
            sys.exit(3)
        if self.Z_approx == None:
            if not force:
                sys.stderr.write("Please call reduce_approx() first or specify 'force=True'")
            else:
                self.reduce_approx()
        dendrogram(self.Z_approx, no_plot=False)
        matplotlib.pyplot.draw()
        matplotlib.pyplot.show()

    def dump_partitions(self):
        part = zip(self.q_vals, self.cuts)
        for q, p in part:
            print q, "->", p

    def dump_partitions_approx(self):
        part = zip(self.q_vals_approx, self.cuts_approx)
        for q, p in part:
            print q, "->", p




