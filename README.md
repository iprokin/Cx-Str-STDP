This is the readme for the model
\
"Endocannabinoid dynamics gate spike-timing dependent depression and potentiation"
\
Yihui Cui, Ilya Prokin, Hao Xu, Bruno Delord, Stephane Genet, Laurent Venance, Hugues Berry
\
<https://elifesciences.org/content/5/e13185>
\
A copy of the code can be found here <https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=187605>.


Abstract
-------------------

Synaptic plasticity is a cardinal cellular mechanism for learning and
memory. The endocannabinoid (eCB) system has emerged as a pivotal
pathway for synaptic plasticity because of its widely characterized
ability to depress synaptic transmission on short- and long-term scales.
Recent reports indicate that eCBs also mediate potentiation of the
synapse. However it is not known how eCB signaling may support such
bidirectionality. Here, we combined electrophysiology experiments with
biophysical modeling to question the mechanisms of eCB bidirectionality
in spike-timing dependent plasticity at corticostriatal synapses. We
demonstrate that STDP outcome is controlled by eCB levels and dynamics:
prolonged and moderate levels of eCB lead to eCB-mediated long-term
depression (eCB-tLTD) while short and large eCB transients produce
eCB-mediated long-term potentiation (eCB-tLTP). Moreover, we show that
eCB-tLTD requires active calcineurin whereas eCB-tLTP necessitates the
activity of presynaptic PKA. Therefore, just like neurotransmitters
glutamate or GABA, eCB form a bidirectional system to encode learning
and memory.


Technical note
-------------------

This package was designed to keep its dependencies to minimum.
However, it requires:

- python2
- gcc
- gfortran
- NumPy
- SciPy (optional, for Gaussian smoothing of STDP curve)
- matplotlib (for plotting)

If you have matplotlib, your default python is python2 and you have
installed NumPy from source or with pip, there is a high chance that all
dependencies are already satisfied.

The equations of RHS of the model are implemented in FORTRAN95 using
FORTRAN77 legacy code of ODEPACK for the numerical integration of ODEs.
The FORTRAN code compiles to python module with gfortran by f2py (part
of NumPy). The module could be imported by python as any
other normal module "import solve_py".
To install dependencies and to compile the module
there are two helper bash scripts.

Note on manual installation:
If you prefer manual install of NumPy, SciPy and matplotlib, refer to
<http://www.scipy.org/install.html>.
To manualy install gcc and gfortran, refer to
<https://gcc.gnu.org/install/binaries.html> and
<https://gcc.gnu.org/wiki/GFortranBinaries>.

This package was tested on Ubuntu, Fedora, Arch Linux and OSX (with
homebrew). Using windows, you will have to manually install
dependencies and manually compile module for python (refer to commands
in file "install_dependencies.sh" and "compile_py_module.sh")


Usage note
-------------------

1) If you have dependencies installed skip to 2). If not, either
manually install dependencies e.g. with your package manager, or to do it
automatically use the script "install_dependencies.sh". The script will
try to detect your system's package manager and suggest packages to
install. To use it, call:\

    ```{.sh}
    bash install_dependencies.sh
    ```

2) To compile module call:

    ```{.sh}
    bash compile_py_module.sh
    ```

3) To produce figures call:

    ```{.sh}
    python2 run_me.py
    ```

4) After 3), newly created folder "pngs" should contain figures with
time traces of model's variables and STDP curves.

5) (optional) modify parameters in "paramets.py" and repeat from 3)


Disclaimer
-------------------

odepack code was taken from <https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html> and split to separate files with [f77split](https://people.sc.fsu.edu/~jburkardt/c_src/f77split/f77split.html).
