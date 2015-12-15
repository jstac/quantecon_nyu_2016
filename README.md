
# QUANTECON_NYU_2016

This is the home page of the course on computational economics I'll be taking
at NYU in the spring semester of 2016.  

Semi-Random quote

>   All this technology  carries risk. There is no faster way for a trading
>   firm to destroy itself  than to deploy a piece of trading software that
>   makes a bad decision over and over in a tight loop. Part of Jane Street's
>   reaction to these  technological risks was to put a very strong focus on
>   building software  that was easily understood--software that was readable.
>   -- Yaron  Minsky, Jane Street

Table of Contents:

* [Syllabus](#syllabus)
    * [Part I: Programming](#part-i-programming)
    * [Part II: Comp Econ Foundations](#part-ii-comp-econ-foundations)
    * [Part III: Applications](#part-iii-applications)
* [Assessment](#assessment)
* [Additional Resources](#additional-resources)



## Syllabus 


### Part I: Programming 


#### Introduction

* Scientific programming environments --- what do we want?

* Why Python?  And what is it anyway?

    * Background
        * http://quant-econ.net/py/about_py.html
        * http://www.galvanize.com/blog/2015/10/01/bill-and-melinda-gates-foundation-saves-lives-with-python/

    * Philosophy
        * import this, import antigravity
        * http://legacy.python.org/dev/peps/pep-0020/
        * https://www.reddit.com/r/Python/comments/3s4j6n/zen_of_python_verse_2/

    * Automation and glue language
        * https://github.com/jstac/backup_scripts


* What's Julia?

    * http://julialang.org/blog/2012/02/why-we-created-julia/
    * http://libertystreeteconomics.newyorkfed.org/2015/12/the-frbny-dsge-model-meets-julia.html


* Open Source 

    * Examples of how contributions improve on standard library

        * http://docs.python-requests.org/en/latest/
        * https://python-programming.courses/general/better-date-and-time-handling-with-arrow/

    * Open science

        * http://devblogs.nvidia.com/parallelforall/open-reproducible-computational-chemistry-python-cuda/
        * http://www.nature.com/news/interactive-notebooks-sharing-the-code-1.16261
        * http://blog.f1000research.com/2014/11/11/what-is-open-science/



#### Coding 101

* UNIX and the UNIX shell
    * http://swcarpentry.github.io/shell-novice/
    * http://swcarpentry.github.io/make-novice/

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

    Hardware
    Vectorized code
    Interpreted / JIT compiled / AOT compiled
    Comparison of loops and limits of vectorization: http://koaning.io/posts/julia_for_loops.html

* C and Fortran
    * GSL
    * http://computationalmodelling.bitbucket.org/tools/FORTRAN.html

* Test driven development: 
    * http://code.tutsplus.com/tutorials/beginning-test-driven-development-in-python--net-30137



#### Core Python

* The REPLs: 
    * Python and IPython shells
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
 
    * https://gist.github.com/sloria/7001839
    * https://leanpub.com/intermediatepython/read
    * http://nbviewer.ipython.org/github/rajathkumarmp/Python-Lectures/blob/master/01.ipynb
    * http://book.pythontips.com/en/latest/
    * http://noeticforce.com/best-free-tutorials-to-learn-python-pdfs-ebooks-online-interactive
    * blog.apcelent.com/python-decorator-tutorial-with-example.html

* Debugging 


#### Scientific Python I: Core Libraries

* General Resources
    * https://github.com/jrjohansson/scientific-python-lectures
    * http://bender.astro.sunysb.edu/classes/python-science/
    * http://computationalmodelling.bitbucket.org/tools/

* NumPy and SciPy  
    * http://quant-econ.net/py/numpy.html
    * http://quant-econ.net/py/scipy.html

* Jupyter 
    * http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb
    * http://blog.dominodatalab.com/lesser-known-ways-of-using-notebooks/
    * https://www.moore.org/newsroom/press-releases/2015/07/07/$6m-for-uc-berkeley-and-cal-poly-to-expand-and-enhance-open-source-software-for-scientific-computing-and-data-science
    * my live notebooks for emet book
    * https://github.com/jupyter/jupyterhub
    * http://nbviewer.ipython.org/
    * http://mindtrove.info/#nb-extensions
    * https://plot.ly/ipython-notebooks/ipython-notebook-tutorial/
    * https://github.com/nicolaskruchten/pyconca/blob/master/jupyter_magic.ipynb


#### Scientific Python II: The Ecosystem

* Sympy
    * http://nbviewer.ipython.org/github/ipython/ipython/blob/master/examples/IPython%20Kernel/SymPy.ipynb

* Webscraping 
    * http://robertwdempsey.com/simple-python-web-scraper-get-pricing-data/

* Pandas 
    * http://geoffboeing.com/2015/11/landscape-us-rents/

* Visualization
    * Matplotlib, plotly, Bokeh
    * http://nbviewer.ipython.org/github/clbarnes/plotstyles/blob/master/plotstyles.ipynb

* Statistics and machine learning
  * https://www.youtube.com/watch?v=5W715nfJNJw
  * PyMC, Statsmodels
  * http://scikit-learn.org/stable/related_projects.html
  * https://www.youtube.com/watch?v=L7R4HUQ-eQ0&feature=youtu.be
  * Seaborn

* Parallel processing 
    http://ufora.github.io/ufora/

* Blaze? 

* AOT compilers
    * Cython
    * Nuitka

* Numba and other JIT compilers 
    * http://blog.pyston.org/2015/11/03/102/
    * http://nbviewer.ipython.org/github/postelrich/fin_examples/blob/master/cva/cva1.ipynb

* Wrappers
    * https://github.com/wjakob/pybind11
    * f2py and related solutions (https://www.euroscipy.org/2015/schedule/presentation/58/)

* NetworkX



#### Julia


* General
    * http://doodlingindata.com/2015/08/11/writing-good-julia-functions/
    * http://www.slideshare.net/acidflask/an-introduction-to-julia
    * http://computationalmodelling.bitbucket.org/tools/Julia.html

* Libraries
    * Distributions.jl 
    * Gadfly





### Part II: Comp Econ Foundations



#### Functional Analysis

* Metric / Banach / Hilbert space
    - L2
    - Space of bounded functions (cbS is a closed subset)

* Banach contraction mapping theorem
* Neumann series lemma


#### Markov Dynamics I: Finite State

* Finite MCs and Daisuke's code
    - The Dobrushin coefficient
    - A simple coupling argument


#### Markov Dynamics II: General State

* General state spaces
    - Feller chains, Boundedness in prob
    - Monotone methods
* LLN and CLT
* applications like ARCH, poverty traps, STAR, MCMC
* Look ahead method
    - examples in lae_extension?
    - examples in poverty traps survey??




#### Dynamic Programming

* Fundamental theory
    * The principle of optimality
    * VFI
    * Howard's policy iteration algorithm

* Approximation
    * Preserving the contraction property
    * MC for integrals

* Applications (see TE paper, monotone LLN)



### Part III: Applications


#### DP II: Extensions

* Weighted sup norm approach

* The Coleman operator
    * Benhabib wealth dist heavy tails

* Recursive and risk sensitive preferences
    * Stochastic Optimal Growth Model with Risk Sensitive Preferences N Bäuerle, A Jaśkiewicz - arXiv preprint arXiv:1509.05638, 2015


#### Solving Forward Looking Models

* L2 methods

* Asset Pricing


#### Optimal Stopping

* Reservation rule operator
    *  Theory 
    *  Applications


#### Coase's Theory of the Firm

* Theory 
* Implementation

#### APS 

* Maybe








## Assessment

The full details are yet to be filled in but the marks will be spread across

* A class project, options for which include

    * Replicate some published research in Python or Julia

    * Investigate a topic of your own interest

* Presentations

    1. A library or programming topic 

    1. Your class project, towards the end of the course

* Homework assignments

* Participation 
    
    * Attending presentations of your classmates


##### Notes

A completed class project is a GitHub repository containing

* Code

* A PDF document that looks like a short research paper describing the
  project, providing analysis and reporting results

Good projects demonstrate proficiency with 

* good programming style

* techical material discussed in the coure




## Additional Resources

* Jupyter
    * https://github.com/bloomberg/bqplot
    * https://cloud.google.com/datalab/
    * https://github.com/lambdalisue/jupyter-vim-binding
    *  https://github.com/jupyter/jupyterhub
    *  http://mybinder.org/

* Language comparisons
    * http://sebastianraschka.com/Articles/2014_matlab_vs_numpy.html
    * http://scottsievert.github.io/blog/2015/09/01/matlab-to-python/
    * https://www.ibm.com/developerworks/community/blogs/jfp/entry/Python_Meets_Julia_Micro_Performance?lang=en

* Python, general 
    * https://github.com/s16h/py-must-watch
    * https://www.reddit.com/r/Python/comments/3m3ll9/where_python_is_used_in_industry_other_than_webdev/
    * http://bruceeckel.github.io/2015/08/29/what-i-do/

Vectorization:
    * http://blog.datascience.com/straightening-loops-how-to-vectorize-data-aggregation-with-pandas-and-numpy/



