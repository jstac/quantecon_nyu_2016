# quantecon_nyu_2016
Quantitative Economics

Intro
======

What makes a good scientific programming environment?  Suggestions?

- speed
- productivity
- large network / many libraries
- scalability?
- fun
- Readability: 
    All this technology  carries risk. There is no faster way for a trading
    firm to destroy itself  than to deploy a piece of trading software that
    makes a bad decision over and over in a tight loop. Part of Jane Street's
    reaction to these  technological risks was to put a very strong focus on
    building software  that was easily understood—software that was readable.
    -- Yaron  Minsky, Jane Street



What is Python?
====================

zen of Python
    * and follow ups: https://www.reddit.com/r/Python/comments/3s4j6n/zen_of_python_verse_2/
        - Number of bugs is linear with number of code lines (a plus of using Python)

    * import this, import antigravity

http://www.galvanize.com/blog/2015/10/01/bill-and-melinda-gates-foundation-saves-lives-with-python/

https://blog.hartleybrody.com/python-style-guide/

Benefits of open source
    * Examples of how contributions improve on standard library
        * http://docs.python-requests.org/en/latest/
        * https://python-programming.courses/general/better-date-and-time-handling-with-arrow/
    - open science
        * http://devblogs.nvidia.com/parallelforall/open-reproducible-computational-chemistry-python-cuda/
        * http://www.nature.com/news/interactive-notebooks-sharing-the-code-1.16261
        * http://blog.f1000research.com/2014/11/11/what-is-open-science/



Programming Concepts
=============================


Why OOP?  It's like structs -- clearly useful -- with lazy evaluation

The beauty of introspection / IPython / --- use MarkovChain instance mc.[tab] as an example

Test driven development: http://code.tutsplus.com/tutorials/beginning-test-driven-development-in-python--net-30137



Programming Practice
=============================

* UNIX and the UNIX shell
    http://swcarpentry.github.io/shell-novice/

* Automation and scripting
    https://github.com/jstac/backup_scripts
    http://swcarpentry.github.io/make-novice/

* Editing and Vim
    https://danielmiessler.com/study/vim/
    https://realpython.com/blog/python/vim-and-python-a-match-made-in-heaven/
    http://vim-adventures.com/
    http://www.openvim.com/

* tmux
    http://tangosource.com/blog/a-tmux-crash-course-tips-and-tweaks/

* Version control
    https://github.com/swcarpentry/git-novice
    http://gitimmersion.com/
    http://luisbg.blogalia.com//historias/76017 --- Git cheatsheet

* Debugging 

* Profiling

* General software engineering skills
    * http://software-carpentry.org/



Python Core Language
===========================
 
https://google.github.io/styleguide/pyguide.html

    - note how the discourage global vars

* https://leanpub.com/intermediatepython/read

* http://nbviewer.ipython.org/github/rajathkumarmp/Python-Lectures/blob/master/01.ipynb

* http://book.pythontips.com/en/latest/

* List of Python tutorials:
    http://noeticforce.com/best-free-tutorials-to-learn-python-pdfs-ebooks-online-interactive

* Decorators
    blog.apcelent.com/python-decorator-tutorial-with-example.html
    



Scientific Programming
============================

Julia vs Python and Python optimization resources 
    * https://www.ibm.com/developerworks/community/blogs/jfp/entry/Python_Meets_Julia_Micro_Performance?lang=en

Matlab vs NumPy
    http://sebastianraschka.com/Articles/2014_matlab_vs_numpy.html
    http://scottsievert.github.io/blog/2015/09/01/matlab-to-python/


* Resources
    http://bender.astro.sunysb.edu/classes/python-science/
    http://computationalmodelling.bitbucket.org/tools/


* Speed and Efficiency

    Hardware
    Vectorized code
    Interpreted / JIT compiled / AOT compiled
    Comparison of loops and limits of vectorization: http://koaning.io/posts/julia_for_loops.html

* Use Matsuyama synchronization basins of attraction instead of Julia sets to illustrate



Programming in C
==================

* GSL

* Suggested exercise: use simulation to confirm that OLS slope coeff in linear
  AR1 regression is biased downwards.


Programming Topics and Libraries
====================================


Python
-------


* Core language
    * https://gist.github.com/sloria/7001839

* General scientific Python
    * https://github.com/jrjohansson/scientific-python-lectures

* Jupyter and IPython 
    * http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb
    * http://blog.dominodatalab.com/lesser-known-ways-of-using-notebooks/
    * https://www.moore.org/newsroom/press-releases/2015/07/07/$6m-for-uc-berkeley-and-cal-poly-to-expand-and-enhance-open-source-software-for-scientific-computing-and-data-science
    *  https://github.com/jupyter/jupyterhub
    *  http://mybinder.org/
    * my live notebooks for emet book
    * https://github.com/jupyter/jupyterhub
    * https://github.com/jupyter/nbgrader
    * http://nbviewer.ipython.org/
    * https://github.com/bloomberg/bqplot
    * https://github.com/lambdalisue/jupyter-vim-binding
    * http://mindtrove.info/#nb-extensions
    * https://cloud.google.com/datalab/
    * https://plot.ly/ipython-notebooks/ipython-notebook-tutorial/
    * https://github.com/nicolaskruchten/pyconca/blob/master/jupyter_magic.ipynb

* sympy
    * http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/IPython%20Kernel/SymPy.ipynb

* webscraping: http://robertwdempsey.com/simple-python-web-scraper-get-pricing-data/

* pandas -- matt?
    * http://geoffboeing.com/2015/11/landscape-us-rents/

* Matplotlib, plotly, Bokeh
    * http://nbviewer.ipython.org/github/clbarnes/plotstyles/blob/master/plotstyles.ipynb

* Bayesian stats
  * https://www.youtube.com/watch?v=5W715nfJNJw
  * PyMC

* parallel processing in Python 
    http://ufora.github.io/ufora/

* NumPy and SciPy  -- me?

* scikit learn
    * http://scikit-learn.org/stable/related_projects.html
    * https://www.youtube.com/watch?v=L7R4HUQ-eQ0&feature=youtu.be

* Statsmodels, Patsy

* Blaze? 

* Cython, nuitka and other AOT compilers

* Numba and other JIT compilers 
    * http://blog.pyston.org/2015/11/03/102/
    http://nbviewer.ipython.org/github/postelrich/fin_examples/blob/master/cva/cva1.ipynb

* Wrappers
    * https://github.com/wjakob/pybind11
    * f2py and related solutions (https://www.euroscipy.org/2015/schedule/presentation/58/)

* NetworkX

* Seaborn

* Markov chains and MDPs -- me




Julia
===========


http://www.slideshare.net/acidflask/an-introduction-to-julia
http://computationalmodelling.bitbucket.org/tools/Julia.html

* Distributions.jl -- I'll do it?
* Gadfly

C and Fortran
================

http://computationalmodelling.bitbucket.org/tools/FORTRAN.html



Resources
============

Must watch vids:
https://github.com/s16h/py-must-watch

Julia
http://doodlingindata.com/2015/08/11/writing-good-julia-functions/

Discussion of Python
http://bruceeckel.github.io/2015/08/29/what-i-do/

Vectorization:
http://blog.datascience.com/straightening-loops-how-to-vectorize-data-aggregation-with-pandas-and-numpy/

Discussion of speed:
https://www.reddit.com/r/Python/comments/3m3ll9/where_python_is_used_in_industry_other_than_webdev/




Applications / Replications / Technical Topics
================================================

Foundations
------------------------------

* Metric / Banach / Hilbert space
    - L2
    - Space of bounded functions (cbS is a closed subset)

* Banach contraction mapping theorem
* Neumann series lemma


Markov Dynamics
-----------------

* Finite MCs and Daisuke's code
    - The Dobrushin coefficient
    - A simple coupling argument
* General state spaces
    - Feller chains, Boundedness in prob
    - Monotone methods
* LLN and CLT
* applications like ARCH, poverty traps, STAR, MCMC
* Look ahead method
    - examples in lae_extension?
    - examples in poverty traps survey??


Dynamic Programming
--------------------

* Fundamental theory
    - The principle of optimality
    - VFI
    - Howard's policy iteration algorithm
* Approximation
    - Preserving the contraction property
    - MC for integrals

* Applications (see TE paper, monotone LLN)


Unbounded DP
-------------

* Weighted sup norm approach


The Coleman operator
---------------------

* Benhabib wealth dist heavy tails


Recursive and Risk Sensitive Preferences
---------------------------------------

* Stochastic Optimal Growth Model with Risk Sensitive Preferences N Bäuerle, A Jaśkiewicz - arXiv preprint arXiv:1509.05638, 2015

* recursive preferences


Asset Pricing
--------------

* Lucas method first
* L2 methods


Optimal Stopping
-------------------

* Reservation rule operator
    - Theory 
    - Applications


Coase
-----

* Theory of the firm


APS 
----

* A possibility





Course Structure and Assessment
=================================

** Presentations

- Everyone presents 
    * a library / programming topic 
    * their project, at the end of the course
        - choose length of presentations acc to amount of time remaining

-- aim is to communicate, not to be clever

** Marks for 

* presentations
* participation at presentations
* one project

Good projects show proficiency with techical material discussed above
