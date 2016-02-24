
# Topics in Computational Economics

[John Stachurski](http://johnstachurski.net/)

This is the home page of ECON-GA 3002, a PhD level course on computational economics to be held at [NYU](http://econ.as.nyu.edu/page/home) in the spring semester of 2016.  

(Note: This document is preliminary and still under development)

Semi-Random quote

>   All this technology  carries risk. There is no faster way for a trading
>   firm to destroy itself  than to deploy a piece of trading software that
>   makes a bad decision over and over in a tight loop. Part of Jane Street's
>   reaction to these  technological risks was to put a very strong focus on
>   building software  that was easily understood--software that was readable.
>
>   -- Yaron  Minsky, Jane Street


Table of Contents:

* [News](#news)
* [References](#references)
* [Prerequisites](#prerequisites)
* [Syllabus](#syllabus)
    * [Part I: Programming](#part-i-programming)
    * [Part II: Comp Econ Foundations](#part-ii-comp-econ-foundations)
    * [Part III: Applications](#part-iii-applications)
* [Assessment](#assessment)
* [Additional Resources](#additional-resources)

## News

Please note that the lecture room has changed to **room 5-75 in the Stern Building**.

The time is unchanged: Friday 9am--11am

Please be sure to bring your laptop 



## References

* http://quant-econ.net/
* Secondary / Useful / Related / Recommended texts
    * Kendall Atkinson and Weimin Han (2009). *Theoretical Numerical Analysis* (3rd ed)
    * Ward Cheney (2001).  *Analysis for Applied Mathematics* 
    * Nancy Stokey and Robert Lucas Jr. (1989) *Recursive Methods in Economic Dynamics*
    * John Stachurski (2009).  *Economic Dynamics: Theory and Computation* 


## Prerequisites

I assume that you have

* At least a bit of programming experience
    * E.g., some experience writing Matlab code or similar
* Econ PhD level quantitative skills, including some familiarity with
    * Linear algebra
    * Basic analysis (sequences, limits, continuity, etc.)
    * Dynamics (diff equations, finite Markov chains, AR(1) processes, etc.)

If you would like to prepare for the course before hand please consider 

* Installing [Linux](http://www.ubuntu.com/desktop) on a [VM](https://www.virtualbox.org/wiki/Linux_Downloads) or in a bootable partition on your laptop 
    * Backup your data first!
    * Help available in the first class
* Build up your [Linux skills](http://manuals.bioinformatics.ucr.edu/home/linux-basics) (and
  [profit](http://www.eweek.com/it-management/demand-for-linux-skills-growing-faster-than-talent-pool-report.html)) 
* Do some exercises in real analysis if you are rusty
    * [These notes](http://math.louisville.edu/~lee/ira/IntroRealAnal.pdf) look like about the right level
* Read the first 3 chapters of [RMT](https://mitpress.mit.edu/books/recursive-macroeconomic-theory-1) if you don't know any Markov chain theory or dynamic programming


## Syllabus 

Below is a sketch of the syllabus for the course.  The details are still
subject to some change.

### Part I: Programming 


#### Introduction

* Scientific programming environments --- what do we want?
    * Speed?
    * Productivity?
    * [Fun?](https://xkcd.com/353/)
* Why [Python](https://www.python.org/)?  And what is it anyway?
    * Background
        * http://quant-econ.net/py/about_py.html
        * http://www.galvanize.com/blog/2015/10/01/bill-and-melinda-gates-foundation-saves-lives-with-python/
    * Philosophy
        * http://legacy.python.org/dev/peps/pep-0020/
        * https://gist.github.com/sloria/7001839
    * The [second best language for everything](http://blog.mikiobraun.de/2013/11/how-python-became-the-language-of-choice-for-data-science.html)
        * https://github.com/jstac/backup_scripts
* What's Julia?
    * http://julialang.org/blog/2012/02/why-we-created-julia/
    * http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html
* Open Source 
    * Examples of how contributions improve on the standard library
        * http://docs.python-requests.org/en/latest/
        * https://python-programming.courses/general/better-date-and-time-handling-with-arrow/
    * Open science
        * http://www.openscience.org/blog/?p=269
        * https://opensource.com/resources/open-science
        * http://www.nature.com/news/interactive-notebooks-sharing-the-code-1.16261
        * http://devblogs.nvidia.com/parallelforall/open-reproducible-computational-chemistry-python-cuda/
    * How can open source produce **better** software than firms acting alone?
        * https://github.com/
        * https://www.moore.org/newsroom/press-releases/2015/07/07/$6m-for-uc-berkeley-and-cal-poly-to-expand-and-enhance-open-source-software-for-scientific-computing-and-data-science
        * https://www.continuum.io/


#### Coding Foundations

* UNIX and the UNIX shell
    * http://swcarpentry.github.io/shell-novice/
* Editing = Vim
    * https://danielmiessler.com/study/vim/
    * https://realpython.com/blog/python/vim-and-python-a-match-made-in-heaven/
    * http://vim-adventures.com/
    * http://www.openvim.com/
* Tmux
    * http://tangosource.com/blog/a-tmux-crash-course-tips-and-tweaks/
* Version control
    * https://github.com/swcarpentry/git-novice
    * http://gitimmersion.com/
    * http://luisbg.blogalia.com//historias/76017 --- Git cheatsheet
* General software engineering skills
    * http://software-carpentry.org/
* Speed and Efficiency
    * Hardware
    * Interpreted / JIT compiled / AOT compiled
    * Vectorized code
* C and Fortran
    * [GSL](http://www.gnu.org/software/gsl/)
    * http://computationalmodelling.bitbucket.org/tools/FORTRAN.html
* Test driven development: 
    * http://code.tutsplus.com/tutorials/beginning-test-driven-development-in-python--net-30137



#### Core Python

* [Getting started](http://quant-econ.net/py/getting_started.html)
    * The REPLs: Python and IPython shells
    * Jupyter
    * The beauty of introspection on the fly
* Basic syntax
    * http://quant-econ.net/py/python_by_example.html
    * http://quant-econ.net/py/python_essentials.html
* OOP. It's like structs with lazy evaluation
    * http://quant-econ.net/py/python_oop.html
    * http://quant-econ.net/py/python_foundations.html
    * http://quant-econ.net/py/python_advanced_features.html
* Python style
    * https://blog.hartleybrody.com/python-style-guide/
    * https://google.github.io/styleguide/pyguide.html
    * https://www.python.org/dev/peps/pep-0008/
* Other general Python resources
    * https://leanpub.com/intermediatepython/read
    * http://nbviewer.ipython.org/github/rajathkumarmp/Python-Lectures/blob/master/01.ipynb
    * http://book.pythontips.com/en/latest/
* Debugging 
    * http://www.scipy-lectures.org/advanced/debugging/


#### Scientific Python I: SciPy and Friends 

* General Resources
    * https://github.com/jrjohansson/scientific-python-lectures
    * http://bender.astro.sunysb.edu/classes/python-science/
    * http://computationalmodelling.bitbucket.org/tools/
* [NumPy and SciPy](http://www.scipy.org/)
    * http://quant-econ.net/py/numpy.html
    * http://quant-econ.net/py/scipy.html
* [Jupyter](http://jupyter.org/)
    * http://nbviewer.ipython.org/
    * http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb
    * http://mindtrove.info/#nb-extensions
    * https://plot.ly/ipython-notebooks/ipython-notebook-tutorial/
    * https://github.com/nicolaskruchten/pyconca/blob/master/jupyter_magic.ipynb
* [Matplotlib](http://matplotlib.org/)
    * http://quant-econ.net/py/matplotlib.html
    * http://nbviewer.ipython.org/github/clbarnes/plotstyles/blob/master/plotstyles.ipynb


#### Scientific Python II: The Ecosystem

* [Pandas](http://pandas.pydata.org/)
    * http://geoffboeing.com/2015/11/landscape-us-rents/
* [Numba](http://numba.pydata.org/) and other JIT compilers 
    * http://blog.pyston.org/2015/11/03/102/
    * http://nbviewer.ipython.org/github/postelrich/fin_examples/blob/master/cva/cva1.ipynb
    * https://www.ibm.com/developerworks/community/blogs/jfp/entry/A_Comparison_Of_C_Julia_Python_Numba_Cython_Scipy_and_BLAS_on_LU_Factorization?lang=en
* AOT compilers
    * [Cython](http://cython.org/)
    * Others (Nuitka?)
* Visualization
    * [Plotly](https://plot.ly/), Bokeh
* Statistics and machine learning
  * https://www.youtube.com/watch?v=5W715nfJNJw
  * PyMC, [Statsmodels](http://statsmodels.sourceforge.net/)
  * http://scikit-learn.org/stable/related_projects.html
  * https://www.youtube.com/watch?v=L7R4HUQ-eQ0&feature=youtu.be
  * [Seaborn](http://stanford.edu/~mwaskom/software/seaborn/)
* Parallel processing 
    * http://www.admin-magazine.com/HPC/Articles/Parallel-Python-with-Joblib
    * http://www.davekuhlman.org/python_multiprocessing_01.html
    * http://ufora.github.io/ufora/
* Blaze 
* Wrappers
    * https://github.com/wjakob/pybind11
    * f2py and related solutions (https://www.euroscipy.org/2015/schedule/presentation/58/)
* [NetworkX](https://networkx.github.io/)
* [Sympy](http://www.sympy.org/en/index.html)
    * http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/IPython%20Kernel/SymPy.ipynb
* Webscraping 
    * http://shop.oreilly.com/product/0636920034391.do
    * http://robertwdempsey.com/simple-python-web-scraper-get-pricing-data/



#### Julia


* General
    * http://julialang.org/
    * http://www.slideshare.net/acidflask/an-introduction-to-julia
    * http://doodlingindata.com/2015/08/11/writing-good-julia-functions/
    * http://computationalmodelling.bitbucket.org/tools/Julia.html
* Libraries
    * [QuantEcon.jl](https://github.com/QuantEcon/QuantEcon.jl)
    * [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)
    * [Gadfly](http://dcjones.github.io/Gadfly.jl/)



### Part II: Comp Econ Foundations


#### Markov Dynamics I: Finite State

* Asymptotics
* The Dobrushin coefficient
* A simple coupling argument
* Code from QuantEcon
* Applications


#### Functional Analysis

* A dash of measure and integration
* Metric / Banach / Hilbert space
    * Space of bounded functions (cbS is a closed subset)
    * The Lp spaces
* Banach contraction mapping theorem
    * Blackwell's sufficient condition
* Orthogonal projections
* Neumann series lemma
* Applications
    * The Lucas 78 asset pricing paper


#### Markov Dynamics II: General State

* General state spaces
    * Feller chains, Boundedness in prob
    * Monotone methods
* LLN and CLT
* Look ahead method
    * examples in lae_extension?
    * examples in poverty traps survey?
* Applications 
    * ARCH, AZ, STAR, MCMC, etc.


#### Solving Forward Looking Models

* L2 methods
* Asset Pricing


#### Dynamic Programming

* Fundamental theory
    * The principle of optimality
    * VFI
    * Howard's policy iteration algorithm
* Approximation
    * Preserving the contraction property
    * MC for integrals
* Weighted sup norm approach


### Part III: Applications


#### DP II: Applications and Extensions

* The Coleman operator
    * [The income fluctuation problem](http://quant-econ.net/py/ifp.html)
    * Benhabib wealth distribution paper, heavy tails
* Recursive and risk sensitive preferences
    * [Stochastic Optimal Growth Model with Risk Sensitive Preferences](http://arxiv.org/abs/1509.05638)
* Other (see TE paper, monotone LLN)


#### Optimal Stopping

* Reservation rule operator
    *  Theory 
    *  Applications


#### Coase's Theory of the Firm

* Theory 
* Implementation



## Assessment

See lecture 1 slides.

### Notes on Class Presentations

All students enrolled in the course must give a 20 minute presentation.
The presentation can be on your class project or on a code library or
algorithm in Julia or Python that you find interesting.  Here are some
suggestions:

* Profiling (see, e.g., [this link](http://pynash.org/2013/03/06/timing-and-profiling.html) or [this one](https://zapier.com/engineering/profiling-python-boss/))
* [scikit-learn](http://scikit-learn.org/stable/) (a machine learning library)
* Unit tests (see, e.g., [here](http://docs.python-guide.org/en/latest/writing/tests/) or [here](https://www.jeffknupp.com/blog/2013/12/09/improve-your-python-understanding-unit-testing/))
* Alternative plotting libraries and their strengths / weaknesses
* [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) (a well-written Julia library)
* Some features of vim or vim plug-in(s) that you find particularly useful
* Techniques for parallel processing
* Interfacing with C and Fortran code in either Python or Julia


### Notes on the Class Project

You should discuss your class project at least briefly with me before you
start.  I am flexible about topics and mainly concerned with quality.

#### Structure of the Project

A completed class project is a GitHub repository containing

* Code
* A Jupyter notebook that pulls all the code together and runs it
* A PDF document that provides analysis and reports results
    * like a short research paper

Good projects demonstrate proficiency with 

* Python or Julia
* Good programming style
* Ideally, the techical material discussed during the course

#### Random Ideas

Here are some very random ideas that I'll add to over the semester.  The links
are to papers, code or discussions of algorithms, quantitative work, etc. that could
be implemented / replicated / improved using Python or Julia.  Feel free to use or ignore.  (Ideally you
will find your own topic according to your own interests.  Please discuss your
topic with me either way).

* [Computing equilibria in dynamic games](https://www.andrew.cmu.edu/user/sevin/sevin/Research_files/Supergame_March_2015_KJ.pdf)
* [Heterogeneous agents in continuous time](http://www.princeton.edu/~moll/HACTproject.htm)
* [Computing Nash equilibria](https://en.wikipedia.org/wiki/Lemke%E2%80%93Howson_algorithm)
* [The stable marriage problem](https://en.wikipedia.org/wiki/Stable_marriage_problem)
* [Krusell-Smith](https://ideas.repec.org/c/dge/qmrbcd/180.html)
* [Krusell-Smith II](http://www.econ.yale.edu/smith/code.htm)
* Numbafy everything in random.py (ask me)
* Numbafy some of the optimization / root finding routines from SciPy (ask me)
* [Assorted code / ideas from Dean Corbae](https://sites.google.com/site/deancorbae/teaching)
* [Assorted code / ideas from Karen Kopecky](http://www.karenkopecky.net/)
* [Assorted code / ideas from Chris Carroll](http://www.econ2.jhu.edu/people/ccarroll/)
* [A paper on dynamics by Kiminori Matsuyama](http://faculty.wcas.northwestern.edu/~kmatsu/Revisiting%20the%20model%20of%20credit%20cycles%20with%20Good%20and%20Bad%20Projects-2016-2.pdf)
* [An econ geography paper by Paul Krugman](https://ideas.repec.org/a/eee/eecrev/v37y1993i2-3p293-298.html)
* Routines from Miranda and Fackler's [CompEcon](http://www4.ncsu.edu/~pfackler/compecon/toolbox.html) toolkit and [textbook](http://www4.ncsu.edu/~pfackler/compecon/)



## Additional Resources

* Jupyter
    * https://github.com/bloomberg/bqplot
    * https://cloud.google.com/datalab/
    * http://blog.dominodatalab.com/lesser-known-ways-of-using-notebooks/
    * https://github.com/jupyter/jupyterhub
    * http://mybinder.org/

* Data, machine learning and prediction
    * www.galvanize.com/blog/how-random-forest-modeling-solves-seattles-bikesharing-problem/
    * https://anaconda.org/ikkebr/brazilian-federal-payroll/notebook

* Language comparisons
    * http://sebastianraschka.com/Articles/2014_matlab_vs_numpy.html
    * http://scottsievert.github.io/blog/2015/09/01/matlab-to-python/
    * https://www.ibm.com/developerworks/community/blogs/jfp/entry/Python_Meets_Julia_Micro_Performance?lang=en

* Python, general 
    * https://www.reddit.com/r/Python/comments/3s4j6n/zen_of_python_verse_2/
    * https://github.com/s16h/py-must-watch
    * https://www.reddit.com/r/Python/comments/3m3ll9/where_python_is_used_in_industry_other_than_webdev/
    * http://bruceeckel.github.io/2015/08/29/what-i-do/
    * http://blog.apcelent.com/python-decorator-tutorial-with-example.html
    * http://noeticforce.com/best-free-tutorials-to-learn-python-pdfs-ebooks-online-interactive

Vectorization:
    * http://blog.datascience.com/straightening-loops-how-to-vectorize-data-aggregation-with-pandas-and-numpy/

Good reads
    * http://undsci.berkeley.edu/article/cold_fusion_01
    * https://msdn.microsoft.com/en-us/library/dn568100.aspx


