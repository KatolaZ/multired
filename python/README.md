# multired
/**
*
* Copyright (C) 2015 Vincenzo (Enzo) Nicosia <katolaz@yahoo.it>
*
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.  
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* long with this program.  If not, see  <http://www.gnu.org/licenses/>.
*
*/

This is multired-0.1.


This is a Python implementation of the algorithm for structural
reduction of multi-layer networks based on the Von Neumann Entropy 
and on the Quantum Jensen-Shannon divergence of graphs, as explained 
in:

  M. De Domenico. V. Nicosia, A. Arenas, V. Latora
  "Structural reducibility of multilayer networks", 
  Nat. Commun. 6, 6864 (2015) doi:10.1038/ncomms7864

If you happen to find any use of this code please do not forget to
cite that paper ;-)


--------------------
       INFO
--------------------

The module "multired.py" provides the class "multiplex_red", which
implements the algorithm to reduce a multilayer network described in
the paper cited above. 

In order to use it, you just need to 

 import multired as mr

in your python script and create a multiplex_red object. Please make
sure that "multired.py" is in PYTHONPATH. The constructor requires as
its first argument the path of a file which in turn contains a list of
files (one for each line) where the graph of each layer is to be
found.

The class provides one set of methods which perform the exact
evaluation of the Von Neumann entropy, and another set of methods
(those whose name end with the suffix "_approx") which rely on a
polynomial approximation of the Von Neumann entropy. By default the
approximation is based on a 10th order polynomial fit of x log(x) in
[0,1], but the order of the polynomial can be set through the
parameter "fit_degree" of the constructor.

Several sample scripts can be found in the "test/" directory. You also
find a sample data set in the folder "sample_data/". That is the 4-layer 
Noordin Top Terrorist multiplex network, originally provided in:

  N. Roberts, S. F. Everton, Roberts and Everton "Terrorist
  Data: Noordin Top Terrorist Network" (Subset) (2011).

and extensively studied in:

  F. Battiston, V. Nicosia, V. Latora,
  "Structural measures for multiplex networks",
  Phys. Rev. E 89, 032804 (2014).

Please consider citing those papers if you use that data set in a
scientific work. 

--------------------
    DEPENDENCIES
--------------------

The only strict dependencies are a recent version of Python, Numpy and
SciPy. The methods "draw_dendrogram" and "draw_dendrogram_approx" will
work only if matplotlib is installed.

The module has been tested on a Debian GNU/Linux system, using:

 - Python 2.7.8,  
 - SciPy 0.13.3
 - Numpy 1.8.2
 - matplotlib 1.3.1

but it will almost surely work on other platforms and/or with other
versions of those packages. If you would like to report a working
configuration, just email me (the address is at the beginning of this
file).
