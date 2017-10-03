
## Cookbook: Parallelized *STRUCTURE* analyses on unlinked SNPs

As part of the `ipyrad.analysis` toolkit we've created convenience functions for easily distributing *STRUCTURE* analysis jobs on an HPC cluster, and for doing so in a programmatic and reproducible way. Importantly, *our workflow allows you to easily sample different distributions of unlinked SNPs among replicate analyses*, with the final inferred population structure summarized from a distribution of replicates. We also provide some simple interactive plotting functions to make barplots and slightly fancier figure, like below. 

### A note on Jupyter/IPython
This is a Jupyter notebook, a reproducible and executable document. The code in this notebook is Python (2.7), and should be executed either in a jupyter-notebook, like this one, or in an IPython terminal. Execute each cell in order to reproduce our entire analysis. We make use of the `ipyparallel` Python library to distribute *STRUCTURE* jobs across processers in parallel. If that is confusing, see our [tutorial on using ipcluster with jupyter](). The example data set used in this analysis is from the [empirical example ipyrad tutorial](http://ipyrad.readthedocs.io/pedicularis_.html).

### Required software
You can easily install the required software for this notebook locally using `conda` by running the commented code below in a terminal. If you are working on an HPC cluster you **do not need** administrator privileges to install the software in this way, since it is only installed locally.


```python
## conda install ipyrad -c ipyrad
## conda install structure -c ipyrad
## conda install clumpp -c ipyrad
## conda install toytree -c eaton-lab
```

### Import Python libraries


```python
import ipyrad.analysis as ipa      ## ipyrad analysis toolkit
import ipyparallel as ipp          ## parallel processing
import toyplot                     ## plotting library
```

### Parallel cluster setup
Start an `ipcluster` instance in a separate terminal. An easy way to do this in a jupyter-notebook running on an HPC cluster is to go to your Jupyter dashboard, and click [new], and then [terminal], and run '`ipcluster start`' in that terminal. This will start a local cluster on the compute node you are connected to. See our [ipyparallel tutorial] (coming soon) for further details. 


```python
##
## ipcluster start --n=40
##
```


```python
## get parallel client
ipyclient = ipp.Client()
print "Connected to {} cores".format(len(ipyclient))
```

    Connected to 56 cores


### Enter input and output file locations


```python
## the structure formatted file
strfile = "/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/slane1_51.str"

## an optional mapfile, to sample unlinked SNPs
#mapfile = "./analysis-ipyrad/pedic-full_outfiles/pedic-full.snps.map"

## the directory where outfiles should be written
workdir = "/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/"
```

### Create a *Structure* Class object
Structure is kind of an old fashioned program that requires creating quite a few input files to run, which makes it not very convenient to use in a programmatic and reproducible way. To work around this we've created a convenience wrapper object to make it easy to submit Structure jobs and to summarize their results. 


```python
## create a Structure object
struct = ipa.structure(name="structure-test",
                       data=strfile, 
                       #mapfile=mapfile,
                       workdir=workdir)
```

### Set parameter options for this object
Our Structure object will be used to submit jobs to the cluster. It has associated with it a name, a set of input files, and a large number of parameter settings. You can modify the parameters by setting them like below. You can also use tab-completion to see all of the available options, or print them like below. See the [full structure docs here](http://computing.bio.cam.ac.uk/local/doc/structure.pdf) for further details on the function of each parameter. In support of reproducibility, it is good practice to print both the mainparams and extraparams so it is clear which options you used. 


```python
## set mainparams for object
struct.mainparams.burnin = 10000
struct.mainparams.numreps = 100000

## see all mainparams
print struct.mainparams

## see or set extraparams
print struct.extraparams
```

    burnin             10000               
    extracols          0                   
    label              1                   
    locdata            0                   
    mapdistances       0                   
    markernames        0                   
    markovphase        0                   
    missing            -9                  
    notambiguous       -999                
    numreps            100000              
    onerowperind       0                   
    phased             0                   
    phaseinfo          0                   
    phenotype          0                   
    ploidy             2                   
    popdata            0                   
    popflag            0                   
    recessivealleles   0                   
    
    admburnin           500                 
    alpha               1.0                 
    alphamax            10.0                
    alphapriora         1.0                 
    alphapriorb         2.0                 
    alphapropsd         0.025               
    ancestdist          0                   
    ancestpint          0.9                 
    computeprob         1                   
    echodata            0                   
    fpriormean          0.01                
    fpriorsd            0.05                
    freqscorr           1                   
    gensback            2                   
    inferalpha          1                   
    inferlambda         0                   
    intermedsave        0                   
    lambda_             1.0                 
    linkage             0                   
    locispop            0                   
    locprior            0                   
    locpriorinit        1.0                 
    log10rmax           1.0                 
    log10rmin           -4.0                
    log10rpropsd        0.1                 
    log10rstart         -2.0                
    maxlocprior         20.0                
    metrofreq           10                  
    migrprior           0.01                
    noadmix             0                   
    numboxes            1000                
    onefst              0                   
    pfrompopflagonly    0                   
    popalphas           0                   
    popspecificlambda   0                   
    printlambda         1                   
    printlikes          0                   
    printnet            1                   
    printqhat           0                   
    printqsum           1                   
    randomize           0                   
    reporthitrate       0                   
    seed                12345               
    sitebysite          0                   
    startatpopinfo      0                   
    unifprioralpha      1                   
    updatefreq          10000               
    usepopinfo          0                   
    


### Submit jobs to run on the cluster
The function `run()` distributes jobs to run on the cluster and load-balances the parallel workload. It takes a number of arguments. The first, `kpop`, is the number of populations. The second, `nreps`, is the number of replicated runs to perform. Each rep has a different random seed, and if you entered a mapfile for your Structure object then it will subsample unlinked snps independently in each replicate. The `seed` argument can be used to make the replicate analyses reproducible. The `extraparams.seed` parameter will be generated from this for each replicate. And finally, provide it the `ipyclient` object that we created above. The structure object will store an *asynchronous results object* for each job that is submitted so that we can query whether the jobs are finished yet or not. Using a simple for-loop we'll submit 20 replicate jobs to run at four different values of K. 


```python
## a range of K-values to test
tests = [2,3,4,5,6]
```


```python
## submit batches of 20 replicate jobs for each value of K 
for kpop in tests:
    struct.run(
        kpop=kpop, 
        nreps=20, 
        seed=12345,
        ipyclient=ipyclient,
        )
```

    submitted 20 structure jobs [structure-test-K-2]
    submitted 20 structure jobs [structure-test-K-3]
    submitted 20 structure jobs [structure-test-K-4]
    submitted 20 structure jobs [structure-test-K-5]
    submitted 20 structure jobs [structure-test-K-6]


### Track progress until finished
You can check for finished results by using the `get_clumpp_table()` function, which tries to summarize the finished results files. If no results are ready it will simply print a warning message telling you to wait. If you want the notebook to block/wait until all jobs are finished then execute the `wait()` function of the ipyclient object, like below. 


```python
## see submitted jobs (we query first 10 here)
struct.asyncs[:10]
```




    [<AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>,
     <AsyncResult: _call_structure>]




```python
## query a specific job result by index
if struct.asyncs[20].ready():
    print struct.asyncs[20].result()
```

    ('\n\n----------------------------------------------------\nSTRUCTURE by Pritchard, Stephens and Donnelly (2000)\n     and Falush, Stephens and Pritchard (2003)\n       Code by Pritchard, Falush and Hubisz\n             Version 2.3.4 (Jul 2012)\n----------------------------------------------------\n\n\nReading file "/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.mainparams.txt".\nReading file "/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.extraparams.txt".\nReading file "/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.strfile.txt".\nNumber of alleles per locus: min= 2; ave=2.1; max= 4\n\n\n--------------------------------------\n\nFinished initialization; starting MCMC \n100000 iterations + 10000 burnin\n\n\n--------------------------------------------\nOverall proportion of membership of the\nsample in each of the 3 clusters\n\nInferred Clusters\n  1      2      3    \n0.451  0.220  0.329  \n\n--------------------------------------------\n\n Rep#:   Lambda   Alpha      F1    F2    F3     D1,2  D1,3  D2,3    Ln Like  Est Ln P(D)\n10000:    1.00    0.024    0.246 0.524 0.760    0.030 0.055 0.074   --  \n\nBURNIN completed\n Rep#:   Lambda   Alpha      F1    F2    F3     D1,2  D1,3  D2,3    Ln Like  Est Ln P(D)\n20000:    1.00    0.028    0.233 0.543 0.763    0.031 0.056 0.076   -16498    0 \n30000:    1.00    0.034    0.261 0.505 0.771    0.032 0.056 0.075   -16534    0 \n40000:    1.00    0.023    0.237 0.517 0.777    0.031 0.055 0.073   -16508    0 \n50000:    1.00    0.034    0.275 0.530 0.752    0.032 0.055 0.076   -16540    0 \n60000:    1.00    0.033    0.242 0.526 0.779    0.030 0.053 0.074   -16489    0 \n70000:    1.00    0.021    0.248 0.522 0.762    0.030 0.054 0.076   -16450    0 \n80000:    1.00    0.019    0.261 0.523 0.770    0.031 0.056 0.073   -16540    0 \n90000:    1.00    0.030    0.247 0.538 0.770    0.031 0.056 0.074   -16573    0 \n100000:    1.00    0.036    0.255 0.513 0.767    0.031 0.055 0.076   -16493    0 \n110000:    1.00    0.030    0.246 0.548 0.770    0.030 0.054 0.074   -16444    0 \n\nMCMC completed\n\nInferred ancestry of individuals:\n        Label (%Miss) :  Inferred clusters\n  1 p_001s_02   (16)   :  0.000 0.000 1.000 \n  2 p_001s_03   (42)   :  0.000 0.000 0.999 \n  3 p_002.5s_01   (28)   :  0.058 0.001 0.940 \n  4 p_002.5s_04    (9)   :  0.000 0.000 1.000 \n  5 p_002.5s_07   (31)   :  0.000 0.000 1.000 \n  6 p_002s_01   (27)   :  0.000 0.000 1.000 \n  7 p_002s_03   (44)   :  0.000 0.000 1.000 \n  8 p_002s_05   (14)   :  0.000 0.000 1.000 \n  9 p_002s_06   (12)   :  0.000 0.000 1.000 \n 10 p_002s_08   (10)   :  0.000 0.000 1.000 \n 11 p_002s_09   (26)   :  0.000 0.000 1.000 \n 12 p_002s_12   (25)   :  0.000 0.000 1.000 \n 13 p_002s_13   (22)   :  0.000 0.000 1.000 \n 14 p_002s_14   (19)   :  0.000 0.000 1.000 \n 15 p_002s_16   (19)   :  0.000 0.000 0.999 \n 16 p_004.5s_05   (36)   :  0.997 0.003 0.000 \n 17 p_005s_06   (38)   :  0.999 0.001 0.000 \n 18 p_005s_11   (24)   :  0.999 0.001 0.000 \n 19 p_006s_01   (26)   :  0.994 0.005 0.002 \n 20 p_006s_06   (28)   :  0.996 0.004 0.000 \n 21 p_006s_14   (14)   :  0.998 0.001 0.001 \n 22 p_009s_13   (26)   :  0.921 0.078 0.000 \n 23 p_009s_15   (20)   :  0.965 0.034 0.001 \n 24 p_011s_08   (28)   :  0.993 0.007 0.000 \n 25 p_011s_10   (35)   :  0.999 0.001 0.000 \n 26 p_011s_11   (25)   :  0.958 0.041 0.001 \n 27 p_012s_06   (34)   :  0.988 0.000 0.012 \n 28 p_012s_12   (37)   :  0.998 0.001 0.001 \n 29 p_013s_02   (42)   :  0.986 0.008 0.007 \n 30 p_013s_08   (46)   :  0.950 0.016 0.034 \n 31 p_015s_13    (9)   :  0.001 0.999 0.000 \n 32 p_015s_14   (26)   :  0.001 0.999 0.000 \n 33 p_016s_01   (31)   :  0.000 0.999 0.000 \n 34 p_016s_02   (37)   :  0.001 0.998 0.000 \n 35 p_016s_03   (27)   :  0.001 0.999 0.001 \n 36 p_016s_04   (18)   :  0.000 0.999 0.000 \n 37 p_016s_06   (34)   :  0.000 0.999 0.000 \n 38 p_016s_15   (25)   :  0.001 0.999 0.000 \n 39 p_017s_01   (18)   :  0.106 0.894 0.001 \n 40 p_017s_03   (31)   :  0.005 0.994 0.001 \n 41 p_017s_05   (23)   :  0.101 0.898 0.001 \n 42 p_026s_12   (31)   :  0.999 0.000 0.000 \n 43 p_026s_12r   (29)   :  0.998 0.001 0.000 \n 44 p_026s_14   (34)   :  0.999 0.001 0.000 \n 45 p_027s_01   (14)   :  0.956 0.044 0.000 \n 46 p_027s_03   (30)   :  0.999 0.000 0.001 \n 47 p_027s_06   (17)   :  0.999 0.000 0.000 \n 48 p_029s_10   (22)   :  0.009 0.231 0.760 \n 49 p_031s_11   (13)   :  0.000 0.000 0.999 \n 50 p_1027s_12r   (38)   :  0.946 0.053 0.001 \n 51 p_1027s_20   (33)   :  0.992 0.001 0.007 \n\nCommand line arguments:   structure -m /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.mainparams.txt -e /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.extraparams.txt -K 3 -D 771445218 -N 51 -L 1533 -i /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.strfile.txt -o /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/structure-test-K-3-rep-0 \nInput File:    /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/tmp-structure-test-3-0.strfile.txt\nOutput File:   /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/structure-test-K-3-rep-0_f\n\nRun parameters:\n   51 individuals\n   1533 loci\n   3 populations assumed\n   10000 Burn-in period\n   100000 Reps\nRANDOMIZE turned off\n\n--------------------------------------------\nEstimated Ln Prob of Data   = -17339.3\nMean value of ln likelihood = -16495.1\nVariance of ln likelihood   = 1688.4\nMean value of alpha         = 0.0286\n\nMean value of Fst_1         = 0.2499\nMean value of Fst_2         = 0.5277\nMean value of Fst_3         = 0.7694\n\n\n--------------------------------------------\nOverall proportion of membership of the\nsample in each of the 3 clusters\n\nInferred Clusters\n  1      2      3    \n0.449  0.222  0.329  \n\n--------------------------------------------\nFinal results printed to file /uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/slane1_51_outfiles/structure-test-K-3-rep-0_f\n\n', None)



```python
## block/wait until all jobs finished
ipyclient.wait() 
```




    True



### Summarize replicates with CLUMPP
We ran 20 replicates per K-value hypothesis. We now need to concatenate and purmute those results so they can be summarized. For this we use the software clumpp. The default arguments to clumpp are generally good, but you can modify them the same as structure params, by accessing the `.clumppparams` attribute of your structure object. See the [clumpp documentation](https://web.stanford.edu/group/rosenberglab/software/CLUMPP_Manual.pdf) for more details. If you have a large number of samples (>50) you may wish to use the `largeKgreedy` algorithm (m=3) for faster runtimes. Below we run clumpp for each value of K that we ran structure on. You only need to tell the `get_clumpp_table()` function the value of K and it will find all of the result files given the Structure object's `name` and `workdir`.


```python
## set some clumpp params
struct.clumppparams.m = 3               ## use largegreedy algorithm
struct.clumppparams.greedy_option = 2   ## test nrepeat possible orders
struct.clumppparams.repeats = 10000     ## number of repeats
struct.clumppparams
```




    datatype                  0                   
    every_permfile            0                   
    greedy_option             2                   
    indfile                   0                   
    m                         3                   
    miscfile                  0                   
    order_by_run              1                   
    outfile                   0                   
    override_warnings         0                   
    permfile                  0                   
    permutationsfile          0                   
    permuted_datafile         0                   
    popfile                   0                   
    print_every_perm          0                   
    print_permuted_data       0                   
    print_random_inputorder   0                   
    random_inputorderfile     0                   
    repeats                   10000               
    s                         2                   
    w                         1                   




```python
## run clumpp for each value of K
tables = {}
for kpop in tests:
    tables[kpop] = struct.get_clumpp_table(kpop)
```

    mean scores across 20 replicates.
    mean scores across 20 replicates.
    mean scores across 20 replicates.
    mean scores across 20 replicates.
    mean scores across 20 replicates.


### Sort the table order how you like it
This can be useful if, for example, you want to order the names to be in the same order as tips on your phylogeny. 


```python
## custom sorting order
myorder = [
    "32082_przewalskii", 
    "33588_przewalskii",
    "41478_cyathophylloides", 
    "41954_cyathophylloides", 
    "29154_superba",
    "30686_cyathophylla", 
    "33413_thamno", 
    "30556_thamno", 
    "35236_rex", 
    "40578_rex", 
    "35855_rex",
    "39618_rex", 
    "38362_rex",
]

print "custom ordering"
print tables[4].ix[myorder]
```

### A function for adding an interactive hover to our plots
The function automatically parses the table above for you. It can reorder the individuals based on their membership in each group, or based on an input list of ordered names. It returns the table of data as well as a list with information for making interactive hover boxes, which you can see below by hovering over the plots.  


```python
def hover(table):
    hover = []
    for row in range(table.shape[0]):
        stack = []
        for col in range(table.shape[1]):
            label = "Name: {}\nGroup: {}\nProp: {}"\
                .format(table.index[row], 
                        table.columns[col],
                        table.ix[row, col])
            stack.append(label)
        hover.append(stack)
    return list(hover)
```

### Visualize population structure in barplots 
Hover over the plot to see sample names and info in the hover box. 


```python
for kpop in tests:
    ## parse outfile to a table and re-order it
    table = tables[kpop]
    #table = table.ix[myorder]
    
    ## plot barplot w/ hover
    canvas, axes, mark = toyplot.bars(
                            table, 
                            title=hover(table),
                            width=400, 
                            height=200, 
                            yshow=False,                            
                            style={"stroke": toyplot.color.near_black},
                            )

```


<div align="center" class="toyplot" id="tfb6dc92397644cdca63630e162e7a69f"><svg class="toyplot-canvas-Canvas" height="200.0px" id="tb1320e1df85645ed9c0e68aa71f2ac33" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 400.0 200.0" width="400.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="tcab1393713744d08b37c871bec273b8d"><clipPath id="t87104a73f3a14b56ac1f10670fadf045"><rect height="120.0" width="320.0" x="40.0" y="40.0"></rect></clipPath><g clip-path="url(#t87104a73f3a14b56ac1f10670fadf045)"><g class="toyplot-mark-BarMagnitudes" id="t4fbe83ef5b964d17b92e2125c057a26b" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="150.0"><title>Name: p_001s_02
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="150.0"><title>Name: p_001s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="5.919408059194069" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="144.08059194080593"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0592</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="150.0"><title>Name: p_002s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="150.0"><title>Name: p_002s_09
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="150.0"><title>Name: p_002s_14
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="150.0"><title>Name: p_002s_16
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.009999000099988"><title>Name: p_004.5s_05
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.009999000099988"><title>Name: p_005s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.009999000099988"><title>Name: p_005s_11
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.890010998900124" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.109989001099883"><title>Name: p_006s_01
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.009999000099988"><title>Name: p_006s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.7900209979002" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.209979002099793"><title>Name: p_006s_14
Group: 0
Prop: 0.998</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.009999000099988"><title>Name: p_009s_13
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.870012998700133" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.129987001299874"><title>Name: p_009s_15
Group: 0
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.019998000199976"><title>Name: p_011s_08
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009999000099988"><title>Name: p_011s_10
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.890010998900124" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.109989001099883"><title>Name: p_011s_11
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="95.100489951004903" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="54.899510048995104"><title>Name: p_012s_06
Group: 0
Prop: 0.9511</title></rect><rect class="toyplot-Datum" height="99.890010998900124" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.109989001099883"><title>Name: p_012s_12
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="94.570542945705427" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="55.429457054294573"><title>Name: p_013s_02
Group: 0
Prop: 0.9458</title></rect><rect class="toyplot-Datum" height="94.04059594040595" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="55.95940405959405"><title>Name: p_013s_08
Group: 0
Prop: 0.9405</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.009999000099988"><title>Name: p_015s_13
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009999000099988"><title>Name: p_015s_14
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.009999000099988"><title>Name: p_016s_01
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.019998000199976"><title>Name: p_016s_02
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.890010998900124" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.109989001099883"><title>Name: p_016s_03
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.009999000099988"><title>Name: p_016s_04
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.009999000099988"><title>Name: p_016s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.019998000199976"><title>Name: p_016s_15
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.950004999500052" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.049995000499948"><title>Name: p_017s_01
Group: 0
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.900009999000105" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.099990000999902"><title>Name: p_017s_03
Group: 0
Prop: 0.9991</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.009999000099988"><title>Name: p_017s_05
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.009999000099988"><title>Name: p_026s_12r
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.009999000099988"><title>Name: p_026s_14
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.950004999500052" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.049995000499948"><title>Name: p_027s_01
Group: 0
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.770022997700238" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.229977002299762"><title>Name: p_027s_03
Group: 0
Prop: 0.9978</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.009999000099988"><title>Name: p_027s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="26.407359264073591" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="123.59264073592641"><title>Name: p_029s_10
Group: 0
Prop: 0.2641</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="150.0"><title>Name: p_031s_11
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="99.890010998900124" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.109989001099883"><title>Name: p_1027s_12r
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="97.900209979002085" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="52.099790020997908"><title>Name: p_1027s_20
Group: 0
Prop: 0.9791</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.009999000099988"><title>Name: p_001s_03
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="94.07059294070595" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.009999000099988"><title>Name: p_002.5s_01
Group: 1
Prop: 0.9408</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.009999000099988"><title>Name: p_002.5s_04
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009999000099988"><title>Name: p_002s_03
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.009999000099988"><title>Name: p_002s_06
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.009999000099988"><title>Name: p_002s_08
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.009999000099988"><title>Name: p_002s_09
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.009999000099988"><title>Name: p_002s_14
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.009999000099988"><title>Name: p_002s_16
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.009999000099988"><title>Name: p_004.5s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.0"><title>Name: p_005s_06
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.009999000099988"><title>Name: p_005s_11
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.009999000099988"><title>Name: p_006s_01
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.009999000099988"><title>Name: p_006s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.19998000199980481" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.009999000099988"><title>Name: p_006s_14
Group: 1
Prop: 0.002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.009999000099988"><title>Name: p_009s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.11998800119988573" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.009999000099988"><title>Name: p_009s_15
Group: 1
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.009999000099988"><title>Name: p_011s_08
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009999000099988"><title>Name: p_011s_10
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.009999000099988"><title>Name: p_011s_11
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="4.8895110488951161" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.009999000099988"><title>Name: p_012s_06
Group: 1
Prop: 0.0489</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.009999000099988"><title>Name: p_012s_12
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="5.4194580541945854" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.009999000099988"><title>Name: p_013s_02
Group: 1
Prop: 0.0542</title></rect><rect class="toyplot-Datum" height="5.9494050594940617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="50.009999000099988"><title>Name: p_013s_08
Group: 1
Prop: 0.0595</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.009999000099988"><title>Name: p_015s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009999000099988"><title>Name: p_015s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.009999000099988"><title>Name: p_016s_01
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.009999000099988"><title>Name: p_016s_02
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.009999000099988"><title>Name: p_016s_03
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.009999000099988"><title>Name: p_016s_04
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.009999000099988"><title>Name: p_016s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.009999000099988"><title>Name: p_016s_15
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.019998000199976"><title>Name: p_017s_01
Group: 1
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.009999000099988"><title>Name: p_017s_03
Group: 1
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.009999000099988"><title>Name: p_017s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.0"><title>Name: p_026s_12r
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.009999000099988"><title>Name: p_026s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.039996000399959541" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.009999000099988"><title>Name: p_027s_01
Group: 1
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.21997800219977393" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.009999000099988"><title>Name: p_027s_03
Group: 1
Prop: 0.0022</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.009999000099988"><title>Name: p_027s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="73.582641735826428" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="50.009999000099988"><title>Name: p_029s_10
Group: 1
Prop: 0.7359</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.009999000099988"><title>Name: p_031s_11
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.009999000099988"><title>Name: p_1027s_12r
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="2.0897910208979198" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.009999000099988"><title>Name: p_1027s_20
Group: 1
Prop: 0.0209</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="t9a30bce3917f463592154872ef98766a" transform="translate(50.0,150.0)translate(0,10.0)"><line style="" x1="0" x2="300.0" y1="0" y2="0"></line><g><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(2.941176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">0</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(61.764705882352935,6)translate(0,7.5)"><tspan style="font-size:10.0px">10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(120.58823529411765,6)translate(0,7.5)"><tspan style="font-size:10.0px">20</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(179.41176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">30</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(238.23529411764704,6)translate(0,7.5)"><tspan style="font-size:10.0px">40</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(297.05882352941177,6)translate(0,7.5)"><tspan style="font-size:10.0px">50</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-3.0" y2="4.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-6"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1"], "id": "t4fbe83ef5b964d17b92e2125c057a26b", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0592, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.999, 1.0, 0.998, 1.0, 0.9988, 0.9999, 1.0, 0.999, 0.9511, 0.999, 0.9458, 0.9405, 1.0, 1.0, 1.0, 0.9999, 0.999, 1.0, 1.0, 0.9999, 0.9996, 0.9991, 1.0, 1.0, 1.0, 1.0, 0.9996, 0.9978, 1.0, 0.2641, 0.0, 0.999, 0.9791], [1.0, 1.0, 0.9408, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0001, 0.0, 0.001, 0.0, 0.002, 0.0, 0.0012, 0.0001, 0.0, 0.001, 0.0489, 0.001, 0.0542, 0.0595, 0.0, 0.0, 0.0, 0.0001, 0.001, 0.0, 0.0, 0.0001, 0.0003, 0.0009, 0.0, 0.0, 0.0001, 0.0, 0.0004, 0.0022, 0.0, 0.7359, 1.0, 0.001, 0.0209]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#tfb6dc92397644cdca63630e162e7a69f .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "tfb6dc92397644cdca63630e162e7a69f";
            var axes = {"t9a30bce3917f463592154872ef98766a": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 300.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>



<div align="center" class="toyplot" id="t335cb3927008436bbd0cd8caa69991a9"><svg class="toyplot-canvas-Canvas" height="200.0px" id="t340ef6f0a2574e8db9c3940dfd8f13de" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 400.0 200.0" width="400.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="t19b1d52c3f9f4f1192d120222b0df9b8"><clipPath id="tb3f02607ddba4b8cb4f04b086bdd8e0e"><rect height="120.0" width="320.0" x="40.0" y="40.0"></rect></clipPath><g clip-path="url(#tb3f02607ddba4b8cb4f04b086bdd8e0e)"><g class="toyplot-mark-BarMagnitudes" id="t2724d57deeed459988711b8269f36f15" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000099989999"><title>Name: p_001s_02
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.99000099989999"><title>Name: p_001s_03
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="5.8694130586941355" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="144.13058694130586"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0587</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000099989999"><title>Name: p_002s_03
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="149.99000099989999"><title>Name: p_002s_08
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000099989999"><title>Name: p_002s_09
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.99000099989999"><title>Name: p_002s_14
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.98000199980001"><title>Name: p_002s_16
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="95.410458954104598" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="54.589541045895409"><title>Name: p_004.5s_05
Group: 0
Prop: 0.9542</title></rect><rect class="toyplot-Datum" height="91.040895910408949" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="58.959104089591044"><title>Name: p_005s_06
Group: 0
Prop: 0.9105</title></rect><rect class="toyplot-Datum" height="91.090890910908911" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="58.909109089091082"><title>Name: p_005s_11
Group: 0
Prop: 0.911</title></rect><rect class="toyplot-Datum" height="83.371662833716641" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="66.628337166283359"><title>Name: p_006s_01
Group: 0
Prop: 0.8338</title></rect><rect class="toyplot-Datum" height="90.940905909409054" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="59.059094090590946"><title>Name: p_006s_06
Group: 0
Prop: 0.9095</title></rect><rect class="toyplot-Datum" height="92.110788921107883" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="57.88921107889211"><title>Name: p_006s_14
Group: 0
Prop: 0.9212</title></rect><rect class="toyplot-Datum" height="87.551244875512452" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="62.448755124487541"><title>Name: p_009s_13
Group: 0
Prop: 0.8756</title></rect><rect class="toyplot-Datum" height="85.161483851614847" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="64.838516148385153"><title>Name: p_009s_15
Group: 0
Prop: 0.8517</title></rect><rect class="toyplot-Datum" height="90.640935906409368" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="59.359064093590639"><title>Name: p_011s_08
Group: 0
Prop: 0.9065</title></rect><rect class="toyplot-Datum" height="91.000899910009011" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="58.999100089990989"><title>Name: p_011s_10
Group: 0
Prop: 0.9101</title></rect><rect class="toyplot-Datum" height="88.821117888211177" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="61.178882111788823"><title>Name: p_011s_11
Group: 0
Prop: 0.8883</title></rect><rect class="toyplot-Datum" height="94.830516948305174" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="55.169483051694826"><title>Name: p_012s_06
Group: 0
Prop: 0.9484</title></rect><rect class="toyplot-Datum" height="95.630436956304379" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="54.369563043695621"><title>Name: p_012s_12
Group: 0
Prop: 0.9564</title></rect><rect class="toyplot-Datum" height="93.500649935006507" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="56.499350064993493"><title>Name: p_013s_02
Group: 0
Prop: 0.9351</title></rect><rect class="toyplot-Datum" height="59.014098590140989" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="90.985901409859011"><title>Name: p_013s_08
Group: 0
Prop: 0.5902</title></rect><rect class="toyplot-Datum" height="0.5799420057994098" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.42005799420059"><title>Name: p_015s_13
Group: 0
Prop: 0.0058</title></rect><rect class="toyplot-Datum" height="0.53994600539945736" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.46005399460054"><title>Name: p_015s_14
Group: 0
Prop: 0.0054</title></rect><rect class="toyplot-Datum" height="0.35996400359962877" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.64003599640037"><title>Name: p_016s_01
Group: 0
Prop: 0.0036</title></rect><rect class="toyplot-Datum" height="0.15998400159983817" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.84001599840016"><title>Name: p_016s_02
Group: 0
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.40995900409961905" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.59004099590038"><title>Name: p_016s_03
Group: 0
Prop: 0.0041</title></rect><rect class="toyplot-Datum" height="0.23997600239974304" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.76002399760026"><title>Name: p_016s_04
Group: 0
Prop: 0.0024</title></rect><rect class="toyplot-Datum" height="0.34996500349964776" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.65003499650035"><title>Name: p_016s_06
Group: 0
Prop: 0.0035</title></rect><rect class="toyplot-Datum" height="0.41995800419957163" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.58004199580043"><title>Name: p_016s_15
Group: 0
Prop: 0.0042</title></rect><rect class="toyplot-Datum" height="24.207579242075781" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="125.79242075792422"><title>Name: p_017s_01
Group: 0
Prop: 0.2421</title></rect><rect class="toyplot-Datum" height="14.508549145085482" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="135.49145085491452"><title>Name: p_017s_03
Group: 0
Prop: 0.1451</title></rect><rect class="toyplot-Datum" height="18.548145185481445" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="131.45185481451855"><title>Name: p_017s_05
Group: 0
Prop: 0.1855</title></rect><rect class="toyplot-Datum" height="46.175382461753813" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="103.82461753824619"><title>Name: p_026s_12
Group: 0
Prop: 0.4618</title></rect><rect class="toyplot-Datum" height="45.735426457354265" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="104.26457354264573"><title>Name: p_026s_12r
Group: 0
Prop: 0.4574</title></rect><rect class="toyplot-Datum" height="46.235376462353756" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="103.76462353764624"><title>Name: p_026s_14
Group: 0
Prop: 0.4624</title></rect><rect class="toyplot-Datum" height="59.174082591740827" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="90.825917408259173"><title>Name: p_027s_01
Group: 0
Prop: 0.5918</title></rect><rect class="toyplot-Datum" height="68.933106689331069" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="81.066893310668931"><title>Name: p_027s_03
Group: 0
Prop: 0.6894</title></rect><rect class="toyplot-Datum" height="71.852814718528151" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="78.147185281471849"><title>Name: p_027s_06
Group: 0
Prop: 0.7186</title></rect><rect class="toyplot-Datum" height="14.128587141285863" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="135.87141285871414"><title>Name: p_029s_10
Group: 0
Prop: 0.1413</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.99000099989999"><title>Name: p_031s_11
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="48.945105489451066" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="101.05489451054893"><title>Name: p_1027s_12r
Group: 0
Prop: 0.4895</title></rect><rect class="toyplot-Datum" height="56.174382561743826" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="93.825617438256174"><title>Name: p_1027s_20
Group: 0
Prop: 0.5618</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000099989999"><title>Name: p_001s_02
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.9700029997"><title>Name: p_001s_03
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.55994400559941937" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="143.57064293570645"><title>Name: p_002.5s_01
Group: 1
Prop: 0.0056</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="149.99000099989999"><title>Name: p_002.5s_04
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="149.99000099989999"><title>Name: p_002s_01
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.98000199980001"><title>Name: p_002s_03
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="149.99000099989999"><title>Name: p_002s_05
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000099989999"><title>Name: p_002s_06
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="149.98000199980001"><title>Name: p_002s_08
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.98000199980001"><title>Name: p_002s_09
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="149.99000099989999"><title>Name: p_002s_12
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="149.99000099989999"><title>Name: p_002s_13
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.98000199980001"><title>Name: p_002s_14
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.9700029997"><title>Name: p_002s_16
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="4.5695430456954327" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.019998000199976"><title>Name: p_004.5s_05
Group: 1
Prop: 0.0457</title></rect><rect class="toyplot-Datum" height="8.9391060893910677" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.019998000199976"><title>Name: p_005s_06
Group: 1
Prop: 0.0894</title></rect><rect class="toyplot-Datum" height="8.8991100889910939" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.009999000099988"><title>Name: p_005s_11
Group: 1
Prop: 0.089</title></rect><rect class="toyplot-Datum" height="16.428357164283561" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.199980001999798"><title>Name: p_006s_01
Group: 1
Prop: 0.1643</title></rect><rect class="toyplot-Datum" height="9.0390960903909701" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.019998000199976"><title>Name: p_006s_06
Group: 1
Prop: 0.0904</title></rect><rect class="toyplot-Datum" height="7.7692230776922315" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.119988001199879"><title>Name: p_006s_14
Group: 1
Prop: 0.0777</title></rect><rect class="toyplot-Datum" height="12.438756124387552" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.009999000099988"><title>Name: p_009s_13
Group: 1
Prop: 0.1244</title></rect><rect class="toyplot-Datum" height="14.678532146785322" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.159984001599831"><title>Name: p_009s_15
Group: 1
Prop: 0.1468</title></rect><rect class="toyplot-Datum" height="9.3190680931906797" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.03999600039996"><title>Name: p_011s_08
Group: 1
Prop: 0.0932</title></rect><rect class="toyplot-Datum" height="8.979102089791013" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.019998000199976"><title>Name: p_011s_10
Group: 1
Prop: 0.0898</title></rect><rect class="toyplot-Datum" height="11.058894110588945" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.119988001199879"><title>Name: p_011s_11
Group: 1
Prop: 0.1106</title></rect><rect class="toyplot-Datum" height="3.7796220377962157" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="51.38986101389861"><title>Name: p_012s_06
Group: 1
Prop: 0.0378</title></rect><rect class="toyplot-Datum" height="4.1795820417958112" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.18998100189981"><title>Name: p_012s_12
Group: 1
Prop: 0.0418</title></rect><rect class="toyplot-Datum" height="5.479452054794514" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="51.019898010198979"><title>Name: p_013s_02
Group: 1
Prop: 0.0548</title></rect><rect class="toyplot-Datum" height="35.696430356964299" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="55.289471052894712"><title>Name: p_013s_08
Group: 1
Prop: 0.357</title></rect><rect class="toyplot-Datum" height="99.390060993900619" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.029997000299964"><title>Name: p_015s_13
Group: 1
Prop: 0.994</title></rect><rect class="toyplot-Datum" height="99.450054994500562" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009999000099988"><title>Name: p_015s_14
Group: 1
Prop: 0.9946</title></rect><rect class="toyplot-Datum" height="99.620037996200395" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.019998000199976"><title>Name: p_016s_01
Group: 1
Prop: 0.9963</title></rect><rect class="toyplot-Datum" height="99.790020997900214" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.049995000499948"><title>Name: p_016s_02
Group: 1
Prop: 0.998</title></rect><rect class="toyplot-Datum" height="99.450054994500519" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.139986001399862"><title>Name: p_016s_03
Group: 1
Prop: 0.9946</title></rect><rect class="toyplot-Datum" height="99.750024997500276" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.009999000099988"><title>Name: p_016s_04
Group: 1
Prop: 0.9976</title></rect><rect class="toyplot-Datum" height="99.640035996400371" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.009999000099988"><title>Name: p_016s_06
Group: 1
Prop: 0.9965</title></rect><rect class="toyplot-Datum" height="99.560043995600452" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.019998000199976"><title>Name: p_016s_15
Group: 1
Prop: 0.9957</title></rect><rect class="toyplot-Datum" height="75.712428757124286" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.079992000799926"><title>Name: p_017s_01
Group: 1
Prop: 0.7572</title></rect><rect class="toyplot-Datum" height="85.381461853814642" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.109989001099883"><title>Name: p_017s_03
Group: 1
Prop: 0.8539</title></rect><rect class="toyplot-Datum" height="81.341865813418678" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.109989001099883"><title>Name: p_017s_05
Group: 1
Prop: 0.8135</title></rect><rect class="toyplot-Datum" height="53.784621537846228" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.03999600039996"><title>Name: p_026s_12
Group: 1
Prop: 0.5379</title></rect><rect class="toyplot-Datum" height="54.244575542445759" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.019998000199976"><title>Name: p_026s_12r
Group: 1
Prop: 0.5425</title></rect><rect class="toyplot-Datum" height="53.73462653734628" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029997000299964"><title>Name: p_026s_14
Group: 1
Prop: 0.5374</title></rect><rect class="toyplot-Datum" height="40.755924407559256" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.069993000699917"><title>Name: p_027s_01
Group: 1
Prop: 0.4076</title></rect><rect class="toyplot-Datum" height="30.666933306669328" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.399960003999603"><title>Name: p_027s_03
Group: 1
Prop: 0.3067</title></rect><rect class="toyplot-Datum" height="28.077192280771932" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.069993000699917"><title>Name: p_027s_06
Group: 1
Prop: 0.2808</title></rect><rect class="toyplot-Datum" height="12.028797120287976" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="123.84261573842616"><title>Name: p_029s_10
Group: 1
Prop: 0.1203</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.92000799920007"><title>Name: p_031s_11
Group: 1
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="50.894910508949089" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.159984001599845"><title>Name: p_1027s_12r
Group: 1
Prop: 0.509</title></rect><rect class="toyplot-Datum" height="42.015798420157985" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="51.809819018098189"><title>Name: p_1027s_20
Group: 1
Prop: 0.4202</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.970002999700014" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.019998000199976"><title>Name: p_001s_02
Group: 2
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.960003999600019" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.009999000099988"><title>Name: p_001s_03
Group: 2
Prop: 0.9997</title></rect><rect class="toyplot-Datum" height="93.560643935606464" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.009999000099988"><title>Name: p_002.5s_01
Group: 2
Prop: 0.9357</title></rect><rect class="toyplot-Datum" height="99.990000999899991" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.0"><title>Name: p_002.5s_04
Group: 2
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.019998000199976"><title>Name: p_002.5s_07
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999899991" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.0"><title>Name: p_002s_01
Group: 2
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.0"><title>Name: p_002s_03
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999899991" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.0"><title>Name: p_002s_05
Group: 2
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.009999000099988"><title>Name: p_002s_06
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.0"><title>Name: p_002s_08
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.0"><title>Name: p_002s_09
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999899991" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.0"><title>Name: p_002s_12
Group: 2
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999899991" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.0"><title>Name: p_002s_13
Group: 2
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.0"><title>Name: p_002s_14
Group: 2
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.960003999600019" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.009999000099988"><title>Name: p_002s_16
Group: 2
Prop: 0.9997</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.009999000099988"><title>Name: p_004.5s_05
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.0"><title>Name: p_005s_06
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.0"><title>Name: p_005s_11
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.1899810018998096" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.009999000099988"><title>Name: p_006s_01
Group: 2
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.009999000099988"><title>Name: p_006s_06
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.10998900109989052" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.009999000099988"><title>Name: p_006s_14
Group: 2
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.0"><title>Name: p_009s_13
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.14998500149984295" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.009999000099988"><title>Name: p_009s_15
Group: 2
Prop: 0.0015</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.009999000099988"><title>Name: p_011s_08
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009999000099988"><title>Name: p_011s_10
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999902406" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.019998000199976"><title>Name: p_011s_11
Group: 2
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="1.3798620137986219" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.009999000099988"><title>Name: p_012s_06
Group: 2
Prop: 0.0138</title></rect><rect class="toyplot-Datum" height="0.17998200179982149" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.009999000099988"><title>Name: p_012s_12
Group: 2
Prop: 0.0018</title></rect><rect class="toyplot-Datum" height="1.0098990100989909" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.009999000099988"><title>Name: p_013s_02
Group: 2
Prop: 0.0101</title></rect><rect class="toyplot-Datum" height="5.2794720527947092" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="50.009999000100002"><title>Name: p_013s_08
Group: 2
Prop: 0.0528</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.009999000099988"><title>Name: p_015s_13
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.0"><title>Name: p_015s_14
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.009999000099988"><title>Name: p_016s_01
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.039996000399959541" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.009999000099988"><title>Name: p_016s_02
Group: 2
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.12998700129987384" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.009999000099988"><title>Name: p_016s_03
Group: 2
Prop: 0.0013</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.0"><title>Name: p_016s_04
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.0"><title>Name: p_016s_06
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.009999000099988"><title>Name: p_016s_15
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.069993000699938079" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.009999000099988"><title>Name: p_017s_01
Group: 2
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.009999000099988"><title>Name: p_017s_03
Group: 2
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.10998900109988341" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.0"><title>Name: p_017s_05
Group: 2
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.009999000099988"><title>Name: p_026s_12r
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.009999000099988"><title>Name: p_026s_14
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.059994000599928654" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.009999000099988"><title>Name: p_027s_01
Group: 2
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.38996100389961441" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.009999000099988"><title>Name: p_027s_03
Group: 2
Prop: 0.0039</title></rect><rect class="toyplot-Datum" height="0.049995000499940545" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.019998000199976"><title>Name: p_027s_06
Group: 2
Prop: 0.0005</title></rect><rect class="toyplot-Datum" height="73.83261673832618" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="50.009999000099988"><title>Name: p_029s_10
Group: 2
Prop: 0.7384</title></rect><rect class="toyplot-Datum" height="99.920007999200067" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.0"><title>Name: p_031s_11
Group: 2
Prop: 0.9993</title></rect><rect class="toyplot-Datum" height="0.14998500149984295" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.009999000100002"><title>Name: p_1027s_12r
Group: 2
Prop: 0.0015</title></rect><rect class="toyplot-Datum" height="1.8098190180981888" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.0"><title>Name: p_1027s_20
Group: 2
Prop: 0.0181</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="t404813caf4c8454cb33552a0c4abbc10" transform="translate(50.0,150.0)translate(0,10.0)"><line style="" x1="0" x2="300.0" y1="0" y2="0"></line><g><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(2.941176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">0</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(61.764705882352935,6)translate(0,7.5)"><tspan style="font-size:10.0px">10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(120.58823529411765,6)translate(0,7.5)"><tspan style="font-size:10.0px">20</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(179.41176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">30</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(238.23529411764704,6)translate(0,7.5)"><tspan style="font-size:10.0px">40</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(297.05882352941177,6)translate(0,7.5)"><tspan style="font-size:10.0px">50</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-3.0" y2="4.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-6"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1", "magnitude2"], "id": "t2724d57deeed459988711b8269f36f15", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0001, 0.0001, 0.0587, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0001, 0.0001, 0.0, 0.0, 0.0001, 0.0002, 0.9542, 0.9105, 0.911, 0.8338, 0.9095, 0.9212, 0.8756, 0.8517, 0.9065, 0.9101, 0.8883, 0.9484, 0.9564, 0.9351, 0.5902, 0.0058, 0.0054, 0.0036, 0.0016, 0.0041, 0.0024, 0.0035, 0.0042, 0.2421, 0.1451, 0.1855, 0.4618, 0.4574, 0.4624, 0.5918, 0.6894, 0.7186, 0.1413, 0.0001, 0.4895, 0.5618], [0.0, 0.0002, 0.0056, 0.0001, 0.0, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0457, 0.0894, 0.089, 0.1643, 0.0904, 0.0777, 0.1244, 0.1468, 0.0932, 0.0898, 0.1106, 0.0378, 0.0418, 0.0548, 0.357, 0.994, 0.9946, 0.9963, 0.998, 0.9946, 0.9976, 0.9965, 0.9957, 0.7572, 0.8539, 0.8135, 0.5379, 0.5425, 0.5374, 0.4076, 0.3067, 0.2808, 0.1203, 0.0007, 0.509, 0.4202], [0.9998, 0.9997, 0.9357, 1.0, 0.9999, 1.0, 0.9999, 1.0, 0.9999, 0.9999, 0.9999, 1.0, 1.0, 0.9999, 0.9997, 0.0001, 0.0002, 0.0001, 0.0019, 0.0001, 0.0011, 0.0001, 0.0015, 0.0003, 0.0001, 0.001, 0.0138, 0.0018, 0.0101, 0.0528, 0.0002, 0.0001, 0.0001, 0.0004, 0.0013, 0.0001, 0.0001, 0.0001, 0.0007, 0.001, 0.0011, 0.0003, 0.0001, 0.0002, 0.0006, 0.0039, 0.0005, 0.7384, 0.9993, 0.0015, 0.0181]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#t335cb3927008436bbd0cd8caa69991a9 .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "t335cb3927008436bbd0cd8caa69991a9";
            var axes = {"t404813caf4c8454cb33552a0c4abbc10": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 300.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>



<div align="center" class="toyplot" id="ta463b4dfd4374193b9c46ddabe60c6aa"><svg class="toyplot-canvas-Canvas" height="200.0px" id="tb02d56ddb10848e6be8cf18a4bce3bdc" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 400.0 200.0" width="400.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="t5e57a815f4ab402780470bff94a476b1"><clipPath id="t53c8906c7a9c4a4c822239d9a8e648e7"><rect height="120.0" width="320.0" x="40.0" y="40.0"></rect></clipPath><g clip-path="url(#t53c8906c7a9c4a4c822239d9a8e648e7)"><g class="toyplot-mark-BarMagnitudes" id="t7d2b775ff8d545ad9fbdd5536852a520" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000099989999"><title>Name: p_001s_02
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.1099890010998763" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.89001099890012"><title>Name: p_001s_03
Group: 0
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="3.5496450354964679" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="146.45035496450353"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0355</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000099989999"><title>Name: p_002s_03
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="150.0"><title>Name: p_002s_09
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.99000099989999"><title>Name: p_002s_14
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.98000199980001"><title>Name: p_002s_16
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="1.8698130186981246" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="148.13018698130188"><title>Name: p_004.5s_05
Group: 0
Prop: 0.0187</title></rect><rect class="toyplot-Datum" height="1.6498350164983435" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="148.35016498350166"><title>Name: p_005s_06
Group: 0
Prop: 0.0165</title></rect><rect class="toyplot-Datum" height="1.5098490150984674" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="148.49015098490153"><title>Name: p_005s_11
Group: 0
Prop: 0.0151</title></rect><rect class="toyplot-Datum" height="30.606939306069393" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="119.39306069393061"><title>Name: p_006s_01
Group: 0
Prop: 0.3061</title></rect><rect class="toyplot-Datum" height="1.5498450154984482" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="148.45015498450155"><title>Name: p_006s_06
Group: 0
Prop: 0.0155</title></rect><rect class="toyplot-Datum" height="1.3198680131986862" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="148.68013198680131"><title>Name: p_006s_14
Group: 0
Prop: 0.0132</title></rect><rect class="toyplot-Datum" height="1.7998200179982007" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="148.2001799820018"><title>Name: p_009s_13
Group: 0
Prop: 0.018</title></rect><rect class="toyplot-Datum" height="3.8296170382961634" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="146.17038296170384"><title>Name: p_009s_15
Group: 0
Prop: 0.0383</title></rect><rect class="toyplot-Datum" height="1.6998300169983054" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="148.30016998300169"><title>Name: p_011s_08
Group: 0
Prop: 0.017</title></rect><rect class="toyplot-Datum" height="1.7598240175982482" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="148.24017598240175"><title>Name: p_011s_10
Group: 0
Prop: 0.0176</title></rect><rect class="toyplot-Datum" height="1.189881011898791" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="148.81011898810121"><title>Name: p_011s_11
Group: 0
Prop: 0.0119</title></rect><rect class="toyplot-Datum" height="1.159884011598848" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="148.84011598840115"><title>Name: p_012s_06
Group: 0
Prop: 0.0116</title></rect><rect class="toyplot-Datum" height="1.6098390160984195" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="148.39016098390158"><title>Name: p_012s_12
Group: 0
Prop: 0.0161</title></rect><rect class="toyplot-Datum" height="3.5496450354964679" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="146.45035496450353"><title>Name: p_013s_02
Group: 0
Prop: 0.0355</title></rect><rect class="toyplot-Datum" height="60.413958604139594" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="89.586041395860406"><title>Name: p_013s_08
Group: 0
Prop: 0.6042</title></rect><rect class="toyplot-Datum" height="0.029997000299999854" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.9700029997"><title>Name: p_015s_13
Group: 0
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.98000199980001"><title>Name: p_015s_14
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.99000099989999"><title>Name: p_016s_01
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.9000099990001"><title>Name: p_016s_02
Group: 0
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.98000199980001"><title>Name: p_016s_03
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.99000099989999"><title>Name: p_016s_04
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.99000099989999"><title>Name: p_016s_06
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.9000099990001"><title>Name: p_016s_15
Group: 0
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.5699430056994288" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="149.43005699430057"><title>Name: p_017s_01
Group: 0
Prop: 0.0057</title></rect><rect class="toyplot-Datum" height="0.39996000399960963" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="149.60003999600039"><title>Name: p_017s_03
Group: 0
Prop: 0.004</title></rect><rect class="toyplot-Datum" height="1.3198680131986862" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="148.68013198680131"><title>Name: p_017s_05
Group: 0
Prop: 0.0132</title></rect><rect class="toyplot-Datum" height="95.120487951204879" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="54.879512048795121"><title>Name: p_026s_12
Group: 0
Prop: 0.9513</title></rect><rect class="toyplot-Datum" height="94.870512948705141" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="55.129487051294866"><title>Name: p_026s_12r
Group: 0
Prop: 0.9488</title></rect><rect class="toyplot-Datum" height="95.230476952304784" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="54.769523047695216"><title>Name: p_026s_14
Group: 0
Prop: 0.9524</title></rect><rect class="toyplot-Datum" height="77.792220777922211" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="72.207779222077789"><title>Name: p_027s_01
Group: 0
Prop: 0.778</title></rect><rect class="toyplot-Datum" height="95.090490950904893" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="54.909509049095099"><title>Name: p_027s_03
Group: 0
Prop: 0.951</title></rect><rect class="toyplot-Datum" height="94.080591940805903" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="55.91940805919409"><title>Name: p_027s_06
Group: 0
Prop: 0.9409</title></rect><rect class="toyplot-Datum" height="0.32996700329968576" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="149.67003299670031"><title>Name: p_029s_10
Group: 0
Prop: 0.0033</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.9000099990001"><title>Name: p_031s_11
Group: 0
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="74.632536746325371" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="75.367463253674629"><title>Name: p_1027s_12r
Group: 0
Prop: 0.7464</title></rect><rect class="toyplot-Datum" height="94.540545945405455" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="55.459454054594538"><title>Name: p_1027s_20
Group: 0
Prop: 0.9455</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.870012998700147" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.019998000199976"><title>Name: p_001s_03
Group: 1
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="93.63063693630636" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="52.81971802819718"><title>Name: p_002.5s_01
Group: 1
Prop: 0.9364</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.009999000099988"><title>Name: p_002.5s_04
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009999000099988"><title>Name: p_002s_03
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.019998000199976"><title>Name: p_002s_06
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.019998000199976"><title>Name: p_002s_08
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.980001999800024" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.019998000199976"><title>Name: p_002s_09
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.009999000099988"><title>Name: p_002s_14
Group: 1
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.970002999700029" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.009999000099988"><title>Name: p_002s_16
Group: 1
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="148.12018798120189"><title>Name: p_004.5s_05
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="148.35016498350166"><title>Name: p_005s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="148.49015098490153"><title>Name: p_005s_11
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.20997900209978582" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="119.18308169183082"><title>Name: p_006s_01
Group: 1
Prop: 0.0021</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="148.44015598440157"><title>Name: p_006s_06
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999923723" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="148.58014198580139"><title>Name: p_006s_14
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="148.19018098190182"><title>Name: p_009s_13
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.17998200179982859" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="145.99040095990401"><title>Name: p_009s_15
Group: 1
Prop: 0.0018</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="148.29017098290171"><title>Name: p_011s_08
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="148.23017698230177"><title>Name: p_011s_10
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089991000899942719" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="148.72012798720127"><title>Name: p_011s_11
Group: 1
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.60993900609938123" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="148.23017698230177"><title>Name: p_012s_06
Group: 1
Prop: 0.0061</title></rect><rect class="toyplot-Datum" height="0.15998400159980974" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="148.23017698230177"><title>Name: p_012s_12
Group: 1
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.90990900909906713" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="145.54044595540446"><title>Name: p_013s_02
Group: 1
Prop: 0.0091</title></rect><rect class="toyplot-Datum" height="13.888611138886105" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="75.6974302569743"><title>Name: p_013s_08
Group: 1
Prop: 0.1389</title></rect><rect class="toyplot-Datum" height="0.0099990000999525819" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.96000399960005"><title>Name: p_015s_13
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.9700029997"><title>Name: p_015s_14
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.98000199980001"><title>Name: p_016s_01
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.88001199880011"><title>Name: p_016s_02
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.1099890010998763" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.87001299870013"><title>Name: p_016s_03
Group: 1
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.98000199980001"><title>Name: p_016s_04
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.99000099989999"><title>Name: p_016s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.87001299870013"><title>Name: p_016s_15
Group: 1
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.049995000499961861" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="149.38006199380061"><title>Name: p_017s_01
Group: 1
Prop: 0.0005</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="149.5000499950005"><title>Name: p_017s_03
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.15998400159983817" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="148.52014798520148"><title>Name: p_017s_05
Group: 1
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="54.879512048795121"><title>Name: p_026s_12
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="55.129487051294866"><title>Name: p_026s_12r
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="54.769523047695216"><title>Name: p_026s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="72.197780221977808"><title>Name: p_027s_01
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.39996000399960252" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="54.509549045095497"><title>Name: p_027s_03
Group: 1
Prop: 0.004</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="55.829417058294176"><title>Name: p_027s_06
Group: 1
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="73.78262173782619" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="75.887411258874124"><title>Name: p_029s_10
Group: 1
Prop: 0.7379</title></rect><rect class="toyplot-Datum" height="99.880011998800128" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.019998000199976"><title>Name: p_031s_11
Group: 1
Prop: 0.9989</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="75.267473252674733"><title>Name: p_1027s_12r
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.71992800719927885" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="54.739526047395259"><title>Name: p_1027s_20
Group: 1
Prop: 0.0072</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.009999000099988"><title>Name: p_001s_03
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="2.699730026997301" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.119988001199879"><title>Name: p_002.5s_01
Group: 2
Prop: 0.027</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.0"><title>Name: p_002.5s_04
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009999000099988"><title>Name: p_002s_03
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.009999000099988"><title>Name: p_002s_06
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.019998000199976"><title>Name: p_002s_08
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.009999000099988"><title>Name: p_002s_09
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.0"><title>Name: p_002s_14
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.0"><title>Name: p_002s_16
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="97.780221977802228" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.33996600339966"><title>Name: p_004.5s_05
Group: 2
Prop: 0.9779</title></rect><rect class="toyplot-Datum" height="98.270172982701723" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.079992000799926"><title>Name: p_005s_06
Group: 2
Prop: 0.9828</title></rect><rect class="toyplot-Datum" height="98.370162983701647" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.119988001199879"><title>Name: p_005s_11
Group: 2
Prop: 0.9838</title></rect><rect class="toyplot-Datum" height="68.213178682131797" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.969903009699031"><title>Name: p_006s_01
Group: 2
Prop: 0.6822</title></rect><rect class="toyplot-Datum" height="97.960203979602056" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.479952004799522"><title>Name: p_006s_06
Group: 2
Prop: 0.9797</title></rect><rect class="toyplot-Datum" height="98.420157984201552" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.159984001599845"><title>Name: p_006s_14
Group: 2
Prop: 0.9843</title></rect><rect class="toyplot-Datum" height="92.690730926907321" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="55.49945005499449"><title>Name: p_009s_13
Group: 2
Prop: 0.927</title></rect><rect class="toyplot-Datum" height="84.82151784821518" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="61.168883111688828"><title>Name: p_009s_15
Group: 2
Prop: 0.8483</title></rect><rect class="toyplot-Datum" height="98.140185981401871" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.149985001499843"><title>Name: p_011s_08
Group: 2
Prop: 0.9815</title></rect><rect class="toyplot-Datum" height="98.180181981801823" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.049995000499948"><title>Name: p_011s_10
Group: 2
Prop: 0.9819</title></rect><rect class="toyplot-Datum" height="98.150184981501837" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.569943005699422"><title>Name: p_011s_11
Group: 2
Prop: 0.9816</title></rect><rect class="toyplot-Datum" height="98.200179982001799" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.029997000299964"><title>Name: p_012s_06
Group: 2
Prop: 0.9821</title></rect><rect class="toyplot-Datum" height="98.110188981101885" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.119988001199879"><title>Name: p_012s_12
Group: 2
Prop: 0.9812</title></rect><rect class="toyplot-Datum" height="94.950504949505046" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.589941005899412"><title>Name: p_013s_02
Group: 2
Prop: 0.9496</title></rect><rect class="toyplot-Datum" height="7.669233076692322" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="68.028197180281978"><title>Name: p_013s_08
Group: 2
Prop: 0.0767</title></rect><rect class="toyplot-Datum" height="0.11998800119988573" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.84001599840016"><title>Name: p_015s_13
Group: 2
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.069993000699895447" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.9000099990001"><title>Name: p_015s_14
Group: 2
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.9700029997"><title>Name: p_016s_01
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.30996900309969533" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.57004299570042"><title>Name: p_016s_02
Group: 2
Prop: 0.0031</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.85001499850014"><title>Name: p_016s_03
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.029997000299943011" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.95000499950007"><title>Name: p_016s_04
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.029997000299943011" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.96000399960005"><title>Name: p_016s_06
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.19998000199981902" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.67003299670031"><title>Name: p_016s_15
Group: 2
Prop: 0.002</title></rect><rect class="toyplot-Datum" height="18.078192180781912" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="131.3018698130187"><title>Name: p_017s_01
Group: 2
Prop: 0.1808</title></rect><rect class="toyplot-Datum" height="4.009599040095992" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="145.4904509549045"><title>Name: p_017s_03
Group: 2
Prop: 0.0401</title></rect><rect class="toyplot-Datum" height="9.8490150984901561" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="138.67113288671132"><title>Name: p_017s_05
Group: 2
Prop: 0.0985</title></rect><rect class="toyplot-Datum" height="4.6195380461953803" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.259974002599741"><title>Name: p_026s_12
Group: 2
Prop: 0.0462</title></rect><rect class="toyplot-Datum" height="4.7695230476952233" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.359964003599643"><title>Name: p_026s_12r
Group: 2
Prop: 0.0477</title></rect><rect class="toyplot-Datum" height="4.5695430456954185" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.199980001999798"><title>Name: p_026s_14
Group: 2
Prop: 0.0457</title></rect><rect class="toyplot-Datum" height="4.5795420457954208" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="67.618238176182388"><title>Name: p_027s_01
Group: 2
Prop: 0.0458</title></rect><rect class="toyplot-Datum" height="4.4595540445955422" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.049995000499955"><title>Name: p_027s_03
Group: 2
Prop: 0.0446</title></rect><rect class="toyplot-Datum" height="5.6194380561943831" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.209979002099793"><title>Name: p_027s_06
Group: 2
Prop: 0.0562</title></rect><rect class="toyplot-Datum" height="20.577942205779429" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="55.309469053094695"><title>Name: p_029s_10
Group: 2
Prop: 0.2058</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.009999000099988"><title>Name: p_031s_11
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="4.3195680431956873" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="70.947905209479046"><title>Name: p_1027s_12r
Group: 2
Prop: 0.0432</title></rect><rect class="toyplot-Datum" height="4.4295570442955707" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.309969003099688"><title>Name: p_1027s_20
Group: 2
Prop: 0.0443</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.0"><title>Name: p_001s_03
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999902406" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.019998000199976"><title>Name: p_002.5s_01
Group: 3
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.0"><title>Name: p_002.5s_04
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.0"><title>Name: p_002s_03
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.0"><title>Name: p_002s_06
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.019998000199976"><title>Name: p_002s_08
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.0"><title>Name: p_002s_09
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.0"><title>Name: p_002s_14
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.0"><title>Name: p_002s_16
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.33996600339965966" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.0"><title>Name: p_004.5s_05
Group: 3
Prop: 0.0034</title></rect><rect class="toyplot-Datum" height="0.069993000699938079" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.009999000099988"><title>Name: p_005s_06
Group: 3
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.099990000999902406" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.019998000199976"><title>Name: p_005s_11
Group: 3
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.95990400959904321" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.009999000099988"><title>Name: p_006s_01
Group: 3
Prop: 0.0096</title></rect><rect class="toyplot-Datum" height="0.4799520047995216" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.0"><title>Name: p_006s_06
Group: 3
Prop: 0.0048</title></rect><rect class="toyplot-Datum" height="0.13998600139984774" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.019998000199998"><title>Name: p_006s_14
Group: 3
Prop: 0.0014</title></rect><rect class="toyplot-Datum" height="5.4994500549944902" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.0"><title>Name: p_009s_13
Group: 3
Prop: 0.055</title></rect><rect class="toyplot-Datum" height="11.148885111488859" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.019998000199969"><title>Name: p_009s_15
Group: 3
Prop: 0.1115</title></rect><rect class="toyplot-Datum" height="0.13998600139985484" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.009999000099988"><title>Name: p_011s_08
Group: 3
Prop: 0.0014</title></rect><rect class="toyplot-Datum" height="0.039996000399959541" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009999000099988"><title>Name: p_011s_10
Group: 3
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.55994400559943358" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.009999000099988"><title>Name: p_011s_11
Group: 3
Prop: 0.0056</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.009999000099988"><title>Name: p_012s_06
Group: 3
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.10998900109989052" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.009999000099988"><title>Name: p_012s_12
Group: 3
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.5699430056994359" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.019998000199976"><title>Name: p_013s_02
Group: 3
Prop: 0.0057</title></rect><rect class="toyplot-Datum" height="18.01819818018199" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="50.009999000099988"><title>Name: p_013s_08
Group: 3
Prop: 0.1802</title></rect><rect class="toyplot-Datum" height="99.840015998400162" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.0"><title>Name: p_015s_13
Group: 3
Prop: 0.9985</title></rect><rect class="toyplot-Datum" height="99.900009999000105" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.0"><title>Name: p_015s_14
Group: 3
Prop: 0.9991</title></rect><rect class="toyplot-Datum" height="99.9700029997" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.0"><title>Name: p_016s_01
Group: 3
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.560043995600438" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.009999000099988"><title>Name: p_016s_02
Group: 3
Prop: 0.9957</title></rect><rect class="toyplot-Datum" height="99.850014998500143" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.0"><title>Name: p_016s_03
Group: 3
Prop: 0.9986</title></rect><rect class="toyplot-Datum" height="99.950004999500067" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.0"><title>Name: p_016s_04
Group: 3
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.960003999600048" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.0"><title>Name: p_016s_06
Group: 3
Prop: 0.9997</title></rect><rect class="toyplot-Datum" height="99.660033996600333" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.009999000099988"><title>Name: p_016s_15
Group: 3
Prop: 0.9967</title></rect><rect class="toyplot-Datum" height="81.291870812918717" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.009999000099988"><title>Name: p_017s_01
Group: 3
Prop: 0.813</title></rect><rect class="toyplot-Datum" height="95.470452954704527" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.019998000199976"><title>Name: p_017s_03
Group: 3
Prop: 0.9548</title></rect><rect class="toyplot-Datum" height="88.67113288671132" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.0"><title>Name: p_017s_05
Group: 3
Prop: 0.8868</title></rect><rect class="toyplot-Datum" height="0.23997600239976435" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.019998000199976"><title>Name: p_026s_12
Group: 3
Prop: 0.0024</title></rect><rect class="toyplot-Datum" height="0.34996500349964066" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.009999000100002"><title>Name: p_026s_12r
Group: 3
Prop: 0.0035</title></rect><rect class="toyplot-Datum" height="0.1899810018998096" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.009999000099988"><title>Name: p_026s_14
Group: 3
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="17.6082391760824" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.009999000099988"><title>Name: p_027s_01
Group: 3
Prop: 0.1761</title></rect><rect class="toyplot-Datum" height="0.039996000399952436" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.009999000100002"><title>Name: p_027s_03
Group: 3
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.1899810018998167" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.019998000199976"><title>Name: p_027s_06
Group: 3
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="5.2894710528947186" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="50.019998000199976"><title>Name: p_029s_10
Group: 3
Prop: 0.0529</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.0"><title>Name: p_031s_11
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="20.937906209379058" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.009999000099988"><title>Name: p_1027s_12r
Group: 3
Prop: 0.2094</title></rect><rect class="toyplot-Datum" height="0.29997000299970011" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.009999000099988"><title>Name: p_1027s_20
Group: 3
Prop: 0.003</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="tf6c34cb93a964337b0ffc26c1fdcf061" transform="translate(50.0,150.0)translate(0,10.0)"><line style="" x1="0" x2="300.0" y1="0" y2="0"></line><g><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(2.941176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">0</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(61.764705882352935,6)translate(0,7.5)"><tspan style="font-size:10.0px">10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(120.58823529411765,6)translate(0,7.5)"><tspan style="font-size:10.0px">20</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(179.41176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">30</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(238.23529411764704,6)translate(0,7.5)"><tspan style="font-size:10.0px">40</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(297.05882352941177,6)translate(0,7.5)"><tspan style="font-size:10.0px">50</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-3.0" y2="4.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-6"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1", "magnitude2", "magnitude3"], "id": "t7d2b775ff8d545ad9fbdd5536852a520", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0001, 0.0011, 0.0355, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0002, 0.0187, 0.0165, 0.0151, 0.3061, 0.0155, 0.0132, 0.018, 0.0383, 0.017, 0.0176, 0.0119, 0.0116, 0.0161, 0.0355, 0.6042, 0.0003, 0.0002, 0.0001, 0.001, 0.0002, 0.0001, 0.0001, 0.001, 0.0057, 0.004, 0.0132, 0.9513, 0.9488, 0.9524, 0.778, 0.951, 0.9409, 0.0033, 0.001, 0.7464, 0.9455], [0.9999, 0.9988, 0.9364, 1.0, 1.0, 1.0, 0.9999, 1.0, 0.9999, 0.9999, 0.9999, 1.0, 1.0, 0.9999, 0.9998, 0.0001, 0.0, 0.0, 0.0021, 0.0001, 0.001, 0.0001, 0.0018, 0.0001, 0.0001, 0.0009, 0.0061, 0.0016, 0.0091, 0.1389, 0.0001, 0.0001, 0.0001, 0.0002, 0.0011, 0.0001, 0.0, 0.0003, 0.0005, 0.001, 0.0016, 0.0, 0.0, 0.0, 0.0001, 0.004, 0.0009, 0.7379, 0.9989, 0.001, 0.0072], [0.0, 0.0001, 0.027, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0, 0.0001, 0.0001, 0.9779, 0.9828, 0.9838, 0.6822, 0.9797, 0.9843, 0.927, 0.8483, 0.9815, 0.9819, 0.9816, 0.9821, 0.9812, 0.9496, 0.0767, 0.0012, 0.0007, 0.0001, 0.0031, 0.0002, 0.0003, 0.0003, 0.002, 0.1808, 0.0401, 0.0985, 0.0462, 0.0477, 0.0457, 0.0458, 0.0446, 0.0562, 0.2058, 0.0001, 0.0432, 0.0443], [0.0, 0.0001, 0.001, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0034, 0.0007, 0.001, 0.0096, 0.0048, 0.0014, 0.055, 0.1115, 0.0014, 0.0004, 0.0056, 0.0002, 0.0011, 0.0057, 0.1802, 0.9985, 0.9991, 0.9998, 0.9957, 0.9986, 0.9996, 0.9997, 0.9967, 0.813, 0.9548, 0.8868, 0.0024, 0.0035, 0.0019, 0.1761, 0.0004, 0.0019, 0.0529, 0.0001, 0.2094, 0.003]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#ta463b4dfd4374193b9c46ddabe60c6aa .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "ta463b4dfd4374193b9c46ddabe60c6aa";
            var axes = {"tf6c34cb93a964337b0ffc26c1fdcf061": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 300.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>



<div align="center" class="toyplot" id="t2be5d6b9ff50410fbe5652c3b7947498"><svg class="toyplot-canvas-Canvas" height="200.0px" id="t6f6cd53569da45d98c2f3635a88c57da" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 400.0 200.0" width="400.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="tfc51b47a511c410f8880b09e162f7203"><clipPath id="t2550544f3ff34c2db332aef345564e5a"><rect height="120.0" width="320.0" x="40.0" y="40.0"></rect></clipPath><g clip-path="url(#t2550544f3ff34c2db332aef345564e5a)"><g class="toyplot-mark-BarMagnitudes" id="td87fdcffa9ab4d3fae85af129820c3c9" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="150.0"><title>Name: p_001s_02
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.91000899910009"><title>Name: p_001s_03
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="2.3097690230976582" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="147.69023097690234"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0231</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="150.0"><title>Name: p_002s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="150.0"><title>Name: p_002s_09
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="150.0"><title>Name: p_002s_14
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="150.0"><title>Name: p_002s_16
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.11998800119988573" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="149.88001199880011"><title>Name: p_004.5s_05
Group: 0
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="150.0"><title>Name: p_005s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.079992000799933294" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="149.92000799920007"><title>Name: p_005s_11
Group: 0
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="30.536946305369455" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="119.46305369463055"><title>Name: p_006s_01
Group: 0
Prop: 0.3054</title></rect><rect class="toyplot-Datum" height="0.1999800019997906" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="149.80001999800021"><title>Name: p_006s_06
Group: 0
Prop: 0.002</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="149.93000699930008"><title>Name: p_006s_14
Group: 0
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="149.9000099990001"><title>Name: p_009s_13
Group: 0
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="2.8097190280971915" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="147.19028097190281"><title>Name: p_009s_15
Group: 0
Prop: 0.0281</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="150.0"><title>Name: p_011s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="149.99000099989999"><title>Name: p_011s_10
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.079992000799933294" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="149.92000799920007"><title>Name: p_011s_11
Group: 0
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="149.91000899910009"><title>Name: p_012s_06
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.22997700229976203" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="149.77002299770024"><title>Name: p_012s_12
Group: 0
Prop: 0.0023</title></rect><rect class="toyplot-Datum" height="0.33996600339966676" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="149.66003399660033"><title>Name: p_013s_02
Group: 0
Prop: 0.0034</title></rect><rect class="toyplot-Datum" height="40.355964403559639" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="109.64403559644036"><title>Name: p_013s_08
Group: 0
Prop: 0.4036</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="150.0"><title>Name: p_015s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.99000099989999"><title>Name: p_015s_14
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.99000099989999"><title>Name: p_016s_01
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.91000899910009"><title>Name: p_016s_02
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.98000199980001"><title>Name: p_016s_03
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.029997000299999854" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.9700029997"><title>Name: p_016s_04
Group: 0
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.99000099989999"><title>Name: p_016s_06
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.91000899910009"><title>Name: p_016s_15
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="149.9000099990001"><title>Name: p_017s_01
Group: 0
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.039996000399952436" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="149.96000399960005"><title>Name: p_017s_03
Group: 0
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.58994100589941922" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="149.41005899410058"><title>Name: p_017s_05
Group: 0
Prop: 0.0059</title></rect><rect class="toyplot-Datum" height="99.950004999500052" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.049995000499948"><title>Name: p_026s_12
Group: 0
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.860013998600138" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.139986001399862"><title>Name: p_026s_12r
Group: 0
Prop: 0.9987</title></rect><rect class="toyplot-Datum" height="99.870012998700133" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.129987001299874"><title>Name: p_026s_14
Group: 0
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="78.522147785221478" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="71.477852214778522"><title>Name: p_027s_01
Group: 0
Prop: 0.7853</title></rect><rect class="toyplot-Datum" height="83.581641835816413" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="66.418358164183587"><title>Name: p_027s_03
Group: 0
Prop: 0.8359</title></rect><rect class="toyplot-Datum" height="80.451954804519545" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="69.548045195480455"><title>Name: p_027s_06
Group: 0
Prop: 0.8046</title></rect><rect class="toyplot-Datum" height="0.1899810018998096" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="149.81001899810019"><title>Name: p_029s_10
Group: 0
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.91000899910009"><title>Name: p_031s_11
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="69.383061693830626" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="80.616938306169374"><title>Name: p_1027s_12r
Group: 0
Prop: 0.6939</title></rect><rect class="toyplot-Datum" height="86.191380861913814" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="63.808619138086186"><title>Name: p_1027s_20
Group: 0
Prop: 0.862</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="150.0"><title>Name: p_001s_02
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.9000099990001"><title>Name: p_001s_03
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="1.7498250174982957" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="145.94040595940405"><title>Name: p_002.5s_01
Group: 1
Prop: 0.0175</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000099989999"><title>Name: p_002s_03
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000099989999"><title>Name: p_002s_06
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000099989999"><title>Name: p_002s_09
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.99000099989999"><title>Name: p_002s_14
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.99000099989999"><title>Name: p_002s_16
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="11.998800119988005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="137.88121187881211"><title>Name: p_004.5s_05
Group: 1
Prop: 0.12</title></rect><rect class="toyplot-Datum" height="97.630236976302371" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="52.369763023697622"><title>Name: p_005s_06
Group: 1
Prop: 0.9764</title></rect><rect class="toyplot-Datum" height="96.990300969902989" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="52.92970702929707"><title>Name: p_005s_11
Group: 1
Prop: 0.97</title></rect><rect class="toyplot-Datum" height="54.624537546245392" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="64.838516148385153"><title>Name: p_006s_01
Group: 1
Prop: 0.5463</title></rect><rect class="toyplot-Datum" height="32.05679432056796" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="117.74322567743225"><title>Name: p_006s_06
Group: 1
Prop: 0.3206</title></rect><rect class="toyplot-Datum" height="42.695730426957311" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="107.23427657234276"><title>Name: p_006s_14
Group: 1
Prop: 0.427</title></rect><rect class="toyplot-Datum" height="32.216778322167784" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="117.68323167683232"><title>Name: p_009s_13
Group: 1
Prop: 0.3222</title></rect><rect class="toyplot-Datum" height="50.974902509749029" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="96.21537846215378"><title>Name: p_009s_15
Group: 1
Prop: 0.5098</title></rect><rect class="toyplot-Datum" height="97.360263973602656" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="52.639736026397351"><title>Name: p_011s_08
Group: 1
Prop: 0.9737</title></rect><rect class="toyplot-Datum" height="97.580241975802409" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="52.409759024097589"><title>Name: p_011s_10
Group: 1
Prop: 0.9759</title></rect><rect class="toyplot-Datum" height="98.380161983801599" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="51.53984601539846"><title>Name: p_011s_11
Group: 1
Prop: 0.9839</title></rect><rect class="toyplot-Datum" height="78.892110788921116" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="71.01789821017897"><title>Name: p_012s_06
Group: 1
Prop: 0.789</title></rect><rect class="toyplot-Datum" height="76.612338766123401" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="73.157684231576837"><title>Name: p_012s_12
Group: 1
Prop: 0.7662</title></rect><rect class="toyplot-Datum" height="12.898710128987091" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="136.76132386761324"><title>Name: p_013s_02
Group: 1
Prop: 0.129</title></rect><rect class="toyplot-Datum" height="6.3793620637936215" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="103.26467353264674"><title>Name: p_013s_08
Group: 1
Prop: 0.0638</title></rect><rect class="toyplot-Datum" height="0.34996500349964776" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.65003499650035"><title>Name: p_015s_13
Group: 1
Prop: 0.0035</title></rect><rect class="toyplot-Datum" height="0.039996000399924014" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.95000499950007"><title>Name: p_015s_14
Group: 1
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.98000199980001"><title>Name: p_016s_01
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.53994600539945736" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.37006299370063"><title>Name: p_016s_02
Group: 1
Prop: 0.0054</title></rect><rect class="toyplot-Datum" height="0.019998000199962007" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.96000399960005"><title>Name: p_016s_03
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.019998000199933585" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.95000499950007"><title>Name: p_016s_04
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.059994000599914443" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.93000699930008"><title>Name: p_016s_06
Group: 1
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.27997200279969547" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.63003699630039"><title>Name: p_016s_15
Group: 1
Prop: 0.0028</title></rect><rect class="toyplot-Datum" height="4.6895310468952971" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="145.21047895210481"><title>Name: p_017s_01
Group: 1
Prop: 0.0469</title></rect><rect class="toyplot-Datum" height="6.6893310668932884" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="143.27067293270676"><title>Name: p_017s_03
Group: 1
Prop: 0.0669</title></rect><rect class="toyplot-Datum" height="11.758824117588233" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="137.65123487651235"><title>Name: p_017s_05
Group: 1
Prop: 0.1176</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.03999600039996"><title>Name: p_026s_12
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199983323" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.119988001199879"><title>Name: p_026s_12r
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.109989001099883"><title>Name: p_026s_14
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="5.8194180581941879" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="65.658434156584335"><title>Name: p_027s_01
Group: 1
Prop: 0.0582</title></rect><rect class="toyplot-Datum" height="2.0497950204979531" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="64.368563143685634"><title>Name: p_027s_03
Group: 1
Prop: 0.0205</title></rect><rect class="toyplot-Datum" height="3.4496550344965584" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="66.098390160983897"><title>Name: p_027s_06
Group: 1
Prop: 0.0345</title></rect><rect class="toyplot-Datum" height="20.807919208079198" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="129.00209979002099"><title>Name: p_029s_10
Group: 1
Prop: 0.2081</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.91000899910009"><title>Name: p_031s_11
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="2.0397960203979437" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="78.57714228577143"><title>Name: p_1027s_12r
Group: 1
Prop: 0.0204</title></rect><rect class="toyplot-Datum" height="1.9898010198980103" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="61.818818118188176"><title>Name: p_1027s_20
Group: 1
Prop: 0.0199</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.98000199980001"><title>Name: p_001s_02
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.31996800319967633" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.58004199580043"><title>Name: p_001s_03
Group: 2
Prop: 0.0032</title></rect><rect class="toyplot-Datum" height="2.3297670232976486" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="143.6106389361064"><title>Name: p_002.5s_01
Group: 2
Prop: 0.0233</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.98000199980001"><title>Name: p_002s_03
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000099989999"><title>Name: p_002s_06
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000099989999"><title>Name: p_002s_09
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.98000199980001"><title>Name: p_002s_14
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.04999500049993344" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.94000599940006"><title>Name: p_002s_16
Group: 2
Prop: 0.0005</title></rect><rect class="toyplot-Datum" height="87.701229877012281" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.179982001799829"><title>Name: p_004.5s_05
Group: 2
Prop: 0.8771</title></rect><rect class="toyplot-Datum" height="2.2897710228977104" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.079992000799912"><title>Name: p_005s_06
Group: 2
Prop: 0.0229</title></rect><rect class="toyplot-Datum" height="2.8097190280971915" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.119988001199879"><title>Name: p_005s_11
Group: 2
Prop: 0.0281</title></rect><rect class="toyplot-Datum" height="13.378662133786612" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="51.459854014598541"><title>Name: p_006s_01
Group: 2
Prop: 0.1338</title></rect><rect class="toyplot-Datum" height="67.083291670832921" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.659934006599336"><title>Name: p_006s_06
Group: 2
Prop: 0.6709</title></rect><rect class="toyplot-Datum" height="56.994300569943" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.239976002399764"><title>Name: p_006s_14
Group: 2
Prop: 0.57</title></rect><rect class="toyplot-Datum" height="66.243375662433778" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="51.439856014398551"><title>Name: p_009s_13
Group: 2
Prop: 0.6625</title></rect><rect class="toyplot-Datum" height="39.866013398660137" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="56.349365063493643"><title>Name: p_009s_15
Group: 2
Prop: 0.3987</title></rect><rect class="toyplot-Datum" height="2.5197480251974724" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.119988001199879"><title>Name: p_011s_08
Group: 2
Prop: 0.0252</title></rect><rect class="toyplot-Datum" height="2.3797620237976247" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.029997000299964"><title>Name: p_011s_10
Group: 2
Prop: 0.0238</title></rect><rect class="toyplot-Datum" height="1.1198880111988814" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.419958004199579"><title>Name: p_011s_11
Group: 2
Prop: 0.0112</title></rect><rect class="toyplot-Datum" height="20.487951204879508" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.529947005299462"><title>Name: p_012s_06
Group: 2
Prop: 0.2049</title></rect><rect class="toyplot-Datum" height="22.807719228077183" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.349965003499655"><title>Name: p_012s_12
Group: 2
Prop: 0.2281</title></rect><rect class="toyplot-Datum" height="85.841415858414166" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.919908009199077"><title>Name: p_013s_02
Group: 2
Prop: 0.8585</title></rect><rect class="toyplot-Datum" height="41.705829417058311" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="61.558844115588428"><title>Name: p_013s_08
Group: 2
Prop: 0.4171</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.58004199580043"><title>Name: p_015s_13
Group: 2
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.16998300169987601" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.78002199780019"><title>Name: p_015s_14
Group: 2
Prop: 0.0017</title></rect><rect class="toyplot-Datum" height="0.029997000299943011" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.95000499950007"><title>Name: p_016s_01
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.34006599340066"><title>Name: p_016s_02
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.89001099890012"><title>Name: p_016s_03
Group: 2
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.95000499950007"><title>Name: p_016s_04
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.019998000199990429" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.91000899910009"><title>Name: p_016s_06
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.079992000799933294" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.55004499550046"><title>Name: p_016s_15
Group: 2
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="22.117788221177904" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="123.0926907309269"><title>Name: p_017s_01
Group: 2
Prop: 0.2212</title></rect><rect class="toyplot-Datum" height="3.4896510348965251" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="139.78102189781023"><title>Name: p_017s_03
Group: 2
Prop: 0.0349</title></rect><rect class="toyplot-Datum" height="5.7394260573942688" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="131.91180881911808"><title>Name: p_017s_05
Group: 2
Prop: 0.0574</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 2
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.0099990000999952144" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.109989001099883"><title>Name: p_026s_12r
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.079992000799919083" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029997000299964"><title>Name: p_026s_14
Group: 2
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.70992900709929074" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="64.948505149485044"><title>Name: p_027s_01
Group: 2
Prop: 0.0071</title></rect><rect class="toyplot-Datum" height="13.928607139286065" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.439956004399569"><title>Name: p_027s_03
Group: 2
Prop: 0.1393</title></rect><rect class="toyplot-Datum" height="15.898410158984099" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.199980001999798"><title>Name: p_027s_06
Group: 2
Prop: 0.159</title></rect><rect class="toyplot-Datum" height="3.0196980301969916" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="125.982401759824"><title>Name: p_029s_10
Group: 2
Prop: 0.0302</title></rect><rect class="toyplot-Datum" height="0.019998000199962007" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.89001099890012"><title>Name: p_031s_11
Group: 2
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="7.8192180781921792" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="70.757924207579251"><title>Name: p_1027s_12r
Group: 2
Prop: 0.0782</title></rect><rect class="toyplot-Datum" height="10.758924107589237" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="51.059894010598939"><title>Name: p_1027s_20
Group: 2
Prop: 0.1076</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.970002999700029" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 3
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.570042995700447" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.009999000099988"><title>Name: p_001s_03
Group: 3
Prop: 0.9958</title></rect><rect class="toyplot-Datum" height="93.500649935006521" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.109989001099883"><title>Name: p_002.5s_01
Group: 3
Prop: 0.9351</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.009999000099988"><title>Name: p_002.5s_04
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.970002999700029" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009999000099988"><title>Name: p_002s_03
Group: 3
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.009999000099988"><title>Name: p_002s_06
Group: 3
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.009999000099988"><title>Name: p_002s_08
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.009999000099988"><title>Name: p_002s_09
Group: 3
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 3
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.970002999700029" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.009999000099988"><title>Name: p_002s_14
Group: 3
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.940005999400057" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.0"><title>Name: p_002s_16
Group: 3
Prop: 0.9995</title></rect><rect class="toyplot-Datum" height="0.019998000199983323" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.159984001599845"><title>Name: p_004.5s_05
Group: 3
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.079992000799912"><title>Name: p_005s_06
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.119988001199879"><title>Name: p_005s_11
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.26997300269974289" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="51.189881011898798"><title>Name: p_006s_01
Group: 3
Prop: 0.0027</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.649935006499348"><title>Name: p_006s_06
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.14998500149985"><title>Name: p_006s_14
Group: 3
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="51.439856014398551"><title>Name: p_009s_13
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.15998400159984527" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="56.189381061893798"><title>Name: p_009s_15
Group: 3
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.099990000999902"><title>Name: p_011s_08
Group: 3
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.029997000299964"><title>Name: p_011s_10
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.389961003899607"><title>Name: p_011s_11
Group: 3
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.50994900509948593" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.019998000199976"><title>Name: p_012s_06
Group: 3
Prop: 0.0051</title></rect><rect class="toyplot-Datum" height="0.17998200179981438" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.16998300169984"><title>Name: p_012s_12
Group: 3
Prop: 0.0018</title></rect><rect class="toyplot-Datum" height="0.38996100389961441" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.529947005299462"><title>Name: p_013s_02
Group: 3
Prop: 0.0039</title></rect><rect class="toyplot-Datum" height="6.3393660633936619" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="55.219478052194766"><title>Name: p_013s_08
Group: 3
Prop: 0.0634</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.57004299570042"><title>Name: p_015s_13
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.78002199780019"><title>Name: p_015s_14
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.94000599940006"><title>Name: p_016s_01
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019998000199962007" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.32006799320069"><title>Name: p_016s_02
Group: 3
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.10998900109993315" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.78002199780019"><title>Name: p_016s_03
Group: 3
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.0099990001000094253" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.94000599940006"><title>Name: p_016s_04
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.91000899910009"><title>Name: p_016s_06
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.029997000299971432" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.52004799520049"><title>Name: p_016s_15
Group: 3
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.059994000599928654" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="123.03269673032698"><title>Name: p_017s_01
Group: 3
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.099990000999923723" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="139.68103189681031"><title>Name: p_017s_03
Group: 3
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.1099890010998763" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="131.8018198180182"><title>Name: p_017s_05
Group: 3
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999810036" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.099990000999902"><title>Name: p_026s_12r
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029997000299964"><title>Name: p_026s_14
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.019998000199976218" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="64.928507149285068"><title>Name: p_027s_01
Group: 3
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.41995800419957163" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.019998000199998"><title>Name: p_027s_03
Group: 3
Prop: 0.0042</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.129987001299874"><title>Name: p_027s_06
Group: 3
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="73.012698730126971" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="52.969703029697037"><title>Name: p_029s_10
Group: 3
Prop: 0.7302</title></rect><rect class="toyplot-Datum" height="99.870012998700147" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.019998000199976"><title>Name: p_031s_11
Group: 3
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="0.10998900109989052" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="70.64793520647936"><title>Name: p_1027s_12r
Group: 3
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.79992000799919794" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.259974002599741"><title>Name: p_1027s_20
Group: 3
Prop: 0.008</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009999000099988"><title>Name: p_001s_02
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.0"><title>Name: p_001s_03
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999895301" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.009999000099988"><title>Name: p_002.5s_01
Group: 4
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.009999000099988"><title>Name: p_002.5s_04
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.009999000099988"><title>Name: p_002.5s_07
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.009999000099988"><title>Name: p_002s_01
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009999000099988"><title>Name: p_002s_03
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.009999000099988"><title>Name: p_002s_05
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.009999000099988"><title>Name: p_002s_06
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.009999000099988"><title>Name: p_002s_08
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.009999000099988"><title>Name: p_002s_09
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.009999000099988"><title>Name: p_002s_12
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.009999000099988"><title>Name: p_002s_13
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.009999000099988"><title>Name: p_002s_14
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.0"><title>Name: p_002s_16
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.15998400159984527" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.0"><title>Name: p_004.5s_05
Group: 4
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.069993000699923869" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.009999000099988"><title>Name: p_005s_06
Group: 4
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.099990000999902406" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.019998000199976"><title>Name: p_005s_11
Group: 4
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="1.169883011698829" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.019998000199969"><title>Name: p_006s_01
Group: 4
Prop: 0.0117</title></rect><rect class="toyplot-Datum" height="0.64993500649934788" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.0"><title>Name: p_006s_06
Group: 4
Prop: 0.0065</title></rect><rect class="toyplot-Datum" height="0.12998700129985252" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.019998000199998"><title>Name: p_006s_14
Group: 4
Prop: 0.0013</title></rect><rect class="toyplot-Datum" height="1.4298570142985625" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.009999000099988"><title>Name: p_009s_13
Group: 4
Prop: 0.0143</title></rect><rect class="toyplot-Datum" height="6.1793820617938238" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.009999000099974"><title>Name: p_009s_15
Group: 4
Prop: 0.0618</title></rect><rect class="toyplot-Datum" height="0.079992000799926188" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.019998000199976"><title>Name: p_011s_08
Group: 4
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.029997000299964327" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.0"><title>Name: p_011s_10
Group: 4
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.3799620037996192" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.009999000099988"><title>Name: p_011s_11
Group: 4
Prop: 0.0038</title></rect><rect class="toyplot-Datum" height="0.009999000099988109" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.009999000099988"><title>Name: p_012s_06
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.14998500149984295" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.019998000199998"><title>Name: p_012s_12
Group: 4
Prop: 0.0015</title></rect><rect class="toyplot-Datum" height="0.50994900509948593" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.019998000199976"><title>Name: p_013s_02
Group: 4
Prop: 0.0051</title></rect><rect class="toyplot-Datum" height="5.2094790520947782" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="50.009999000099988"><title>Name: p_013s_08
Group: 4
Prop: 0.0521</title></rect><rect class="toyplot-Datum" height="99.550044995500443" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.019998000199976"><title>Name: p_015s_13
Group: 4
Prop: 0.9956</title></rect><rect class="toyplot-Datum" height="99.77002299770021" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009999000099988"><title>Name: p_015s_14
Group: 4
Prop: 0.9978</title></rect><rect class="toyplot-Datum" height="99.940005999400057" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.0"><title>Name: p_016s_01
Group: 4
Prop: 0.9995</title></rect><rect class="toyplot-Datum" height="99.300069993000719" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.019998000199976"><title>Name: p_016s_02
Group: 4
Prop: 0.9931</title></rect><rect class="toyplot-Datum" height="99.77002299770021" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.009999000099988"><title>Name: p_016s_03
Group: 4
Prop: 0.9978</title></rect><rect class="toyplot-Datum" height="99.940005999400057" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.0"><title>Name: p_016s_04
Group: 4
Prop: 0.9995</title></rect><rect class="toyplot-Datum" height="99.900009999000105" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.009999000099988"><title>Name: p_016s_06
Group: 4
Prop: 0.9991</title></rect><rect class="toyplot-Datum" height="99.500049995000509" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.019998000199976"><title>Name: p_016s_15
Group: 4
Prop: 0.9951</title></rect><rect class="toyplot-Datum" height="73.02269773022698" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.009999000099988"><title>Name: p_017s_01
Group: 4
Prop: 0.7303</title></rect><rect class="toyplot-Datum" height="89.671032896710329" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.009999000099988"><title>Name: p_017s_03
Group: 4
Prop: 0.8968</title></rect><rect class="toyplot-Datum" height="81.791820817918222" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.009999000099988"><title>Name: p_017s_05
Group: 4
Prop: 0.818</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.009999000099988"><title>Name: p_026s_12
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.089991000899914297" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.009999000099988"><title>Name: p_026s_12r
Group: 4
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029997000299964"><title>Name: p_026s_14
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="14.928507149285068" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.0"><title>Name: p_027s_01
Group: 4
Prop: 0.1493</title></rect><rect class="toyplot-Datum" height="0.0099990000999952144" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.009999000100002"><title>Name: p_027s_03
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.10998900109989762" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.019998000199976"><title>Name: p_027s_06
Group: 4
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="2.9497050294970393" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="50.019998000199998"><title>Name: p_029s_10
Group: 4
Prop: 0.0295</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.019998000199976"><title>Name: p_031s_11
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="20.637936206379372" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.009999000099988"><title>Name: p_1027s_12r
Group: 4
Prop: 0.2064</title></rect><rect class="toyplot-Datum" height="0.24997500249975246" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.009999000099988"><title>Name: p_1027s_20
Group: 4
Prop: 0.0025</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="t872a592f7af744ee8578cb0de4a6b7c7" transform="translate(50.0,150.0)translate(0,10.0)"><line style="" x1="0" x2="300.0" y1="0" y2="0"></line><g><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(2.941176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">0</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(61.764705882352935,6)translate(0,7.5)"><tspan style="font-size:10.0px">10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(120.58823529411765,6)translate(0,7.5)"><tspan style="font-size:10.0px">20</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(179.41176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">30</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(238.23529411764704,6)translate(0,7.5)"><tspan style="font-size:10.0px">40</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(297.05882352941177,6)translate(0,7.5)"><tspan style="font-size:10.0px">50</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-3.0" y2="4.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-6"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1", "magnitude2", "magnitude3", "magnitude4"], "id": "td87fdcffa9ab4d3fae85af129820c3c9", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0009, 0.0231, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0012, 0.0, 0.0008, 0.3054, 0.002, 0.0007, 0.001, 0.0281, 0.0, 0.0001, 0.0008, 0.0009, 0.0023, 0.0034, 0.4036, 0.0, 0.0001, 0.0001, 0.0009, 0.0002, 0.0003, 0.0001, 0.0009, 0.001, 0.0004, 0.0059, 0.9996, 0.9987, 0.9988, 0.7853, 0.8359, 0.8046, 0.0019, 0.0009, 0.6939, 0.862], [0.0, 0.0001, 0.0175, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0, 0.0001, 0.0001, 0.12, 0.9764, 0.97, 0.5463, 0.3206, 0.427, 0.3222, 0.5098, 0.9737, 0.9759, 0.9839, 0.789, 0.7662, 0.129, 0.0638, 0.0035, 0.0004, 0.0001, 0.0054, 0.0002, 0.0002, 0.0006, 0.0028, 0.0469, 0.0669, 0.1176, 0.0001, 0.0002, 0.0002, 0.0582, 0.0205, 0.0345, 0.2081, 0.0, 0.0204, 0.0199], [0.0002, 0.0032, 0.0233, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0005, 0.8771, 0.0229, 0.0281, 0.1338, 0.6709, 0.57, 0.6625, 0.3987, 0.0252, 0.0238, 0.0112, 0.2049, 0.2281, 0.8585, 0.4171, 0.0007, 0.0017, 0.0003, 0.0003, 0.0007, 0.0, 0.0002, 0.0008, 0.2212, 0.0349, 0.0574, 0.0003, 0.0001, 0.0008, 0.0071, 0.1393, 0.159, 0.0302, 0.0002, 0.0782, 0.1076], [0.9998, 0.9958, 0.9351, 1.0, 1.0, 1.0, 0.9998, 1.0, 0.9999, 1.0, 0.9999, 1.0, 1.0, 0.9998, 0.9995, 0.0002, 0.0, 0.0, 0.0027, 0.0001, 0.0009, 0.0, 0.0016, 0.0002, 0.0, 0.0003, 0.0051, 0.0018, 0.0039, 0.0634, 0.0001, 0.0, 0.0001, 0.0002, 0.0011, 0.0001, 0.0, 0.0003, 0.0006, 0.001, 0.0011, 0.0, 0.0001, 0.0, 0.0002, 0.0042, 0.0007, 0.7302, 0.9988, 0.0011, 0.008], [0.0, 0.0001, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0016, 0.0007, 0.001, 0.0117, 0.0065, 0.0013, 0.0143, 0.0618, 0.0008, 0.0003, 0.0038, 0.0001, 0.0015, 0.0051, 0.0521, 0.9956, 0.9978, 0.9995, 0.9931, 0.9978, 0.9995, 0.9991, 0.9951, 0.7303, 0.8968, 0.818, 0.0, 0.0009, 0.0, 0.1493, 0.0001, 0.0011, 0.0295, 0.0, 0.2064, 0.0025]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#t2be5d6b9ff50410fbe5652c3b7947498 .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "t2be5d6b9ff50410fbe5652c3b7947498";
            var axes = {"t872a592f7af744ee8578cb0de4a6b7c7": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 300.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>



<div align="center" class="toyplot" id="t81bd05ec41654c65847f98e257b38e14"><svg class="toyplot-canvas-Canvas" height="200.0px" id="t5f2c45a1e459428d844e42378f9b49a2" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 400.0 200.0" width="400.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="t1a377e737ead4689b38eada0a5292fd0"><clipPath id="t453de98d026846279681dcb9b4137e91"><rect height="120.0" width="320.0" x="40.0" y="40.0"></rect></clipPath><g clip-path="url(#t453de98d026846279681dcb9b4137e91)"><g class="toyplot-mark-BarMagnitudes" id="t689bd55618f94f4f965cd2ede1951a6e" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000199960008"><title>Name: p_001s_02
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="150.0"><title>Name: p_001s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="2.3495300939811727" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="147.65046990601883"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0235</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="150.0"><title>Name: p_002s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="150.0"><title>Name: p_002s_09
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.99000199960008"><title>Name: p_002s_14
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.99000199960008"><title>Name: p_002s_16
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="72.075584883023396" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="77.924415116976604"><title>Name: p_004.5s_05
Group: 0
Prop: 0.7209</title></rect><rect class="toyplot-Datum" height="1.9296140771845671" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="148.07038592281543"><title>Name: p_005s_06
Group: 0
Prop: 0.0193</title></rect><rect class="toyplot-Datum" height="4.7290541891621558" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="145.27094581083784"><title>Name: p_005s_11
Group: 0
Prop: 0.0473</title></rect><rect class="toyplot-Datum" height="12.727454509098152" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="137.27254549090185"><title>Name: p_006s_01
Group: 0
Prop: 0.1273</title></rect><rect class="toyplot-Datum" height="92.301539692061596" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="57.698460307938404"><title>Name: p_006s_06
Group: 0
Prop: 0.9232</title></rect><rect class="toyplot-Datum" height="78.934213157368532" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="71.065786842631468"><title>Name: p_006s_14
Group: 0
Prop: 0.7895</title></rect><rect class="toyplot-Datum" height="92.26154769046191" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="57.73845230953809"><title>Name: p_009s_13
Group: 0
Prop: 0.9228</title></rect><rect class="toyplot-Datum" height="40.821835632873416" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="109.17816436712658"><title>Name: p_009s_15
Group: 0
Prop: 0.4083</title></rect><rect class="toyplot-Datum" height="1.4797040591881512" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="148.52029594081185"><title>Name: p_011s_08
Group: 0
Prop: 0.0148</title></rect><rect class="toyplot-Datum" height="1.4397120575884799" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="148.56028794241152"><title>Name: p_011s_10
Group: 0
Prop: 0.0144</title></rect><rect class="toyplot-Datum" height="0.99980003999201017" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="149.00019996000799"><title>Name: p_011s_11
Group: 0
Prop: 0.01</title></rect><rect class="toyplot-Datum" height="0.59988002399521179" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="149.40011997600479"><title>Name: p_012s_06
Group: 0
Prop: 0.006</title></rect><rect class="toyplot-Datum" height="0.99980003999201017" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="149.00019996000799"><title>Name: p_012s_12
Group: 0
Prop: 0.01</title></rect><rect class="toyplot-Datum" height="61.127774445110987" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="88.872225554889013"><title>Name: p_013s_02
Group: 0
Prop: 0.6114</title></rect><rect class="toyplot-Datum" height="18.836232753449309" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="131.16376724655069"><title>Name: p_013s_08
Group: 0
Prop: 0.1884</title></rect><rect class="toyplot-Datum" height="0.079984003199371045" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.92001599680063"><title>Name: p_015s_13
Group: 0
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.1699660067986315" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.83003399320137"><title>Name: p_015s_14
Group: 0
Prop: 0.0017</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.99000199960008"><title>Name: p_016s_01
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.98000399920016"><title>Name: p_016s_02
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.98000399920016"><title>Name: p_016s_03
Group: 0
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="150.0"><title>Name: p_016s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.99000199960008"><title>Name: p_016s_06
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089982003599260452" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.91001799640074"><title>Name: p_016s_15
Group: 0
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="17.60647870425916" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="132.39352129574084"><title>Name: p_017s_01
Group: 0
Prop: 0.1761</title></rect><rect class="toyplot-Datum" height="0.96980603879225669" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="149.03019396120774"><title>Name: p_017s_03
Group: 0
Prop: 0.0097</title></rect><rect class="toyplot-Datum" height="0.48990201959608726" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="149.51009798040391"><title>Name: p_017s_05
Group: 0
Prop: 0.0049</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="150.0"><title>Name: p_026s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="150.0"><title>Name: p_026s_12r
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.029994001199753484" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="149.97000599880025"><title>Name: p_026s_14
Group: 0
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="5.2989402119576141" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="144.70105978804239"><title>Name: p_027s_01
Group: 0
Prop: 0.053</title></rect><rect class="toyplot-Datum" height="0.35992801439712707" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="149.64007198560287"><title>Name: p_027s_03
Group: 0
Prop: 0.0036</title></rect><rect class="toyplot-Datum" height="1.6496700659868111" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="148.35032993401319"><title>Name: p_027s_06
Group: 0
Prop: 0.0165</title></rect><rect class="toyplot-Datum" height="7.298540291941606" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="142.70145970805839"><title>Name: p_029s_10
Group: 0
Prop: 0.073</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.99000199960008"><title>Name: p_031s_11
Group: 0
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.65986802639469033" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="149.34013197360531"><title>Name: p_1027s_12r
Group: 0
Prop: 0.0066</title></rect><rect class="toyplot-Datum" height="0.62987402519496527" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="149.37012597480503"><title>Name: p_1027s_20
Group: 0
Prop: 0.0063</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000199960008"><title>Name: p_001s_02
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.069986002799424796" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.93001399720058"><title>Name: p_001s_03
Group: 1
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="1.3297340531894122" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="146.32073585282942"><title>Name: p_002.5s_01
Group: 1
Prop: 0.0133</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="150.0"><title>Name: p_002.5s_04
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="150.0"><title>Name: p_002s_03
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="150.0"><title>Name: p_002s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="150.0"><title>Name: p_002s_09
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.99000199960008"><title>Name: p_002s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.99000199960008"><title>Name: p_002s_16
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.11997600479904236" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="77.804439112177562"><title>Name: p_004.5s_05
Group: 1
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="148.07038592281543"><title>Name: p_005s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.079984003199371045" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="145.19096180763847"><title>Name: p_005s_11
Group: 1
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="21.765646870625915" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="115.50689862027593"><title>Name: p_006s_01
Group: 1
Prop: 0.2177</title></rect><rect class="toyplot-Datum" height="0.15996800639872788" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="57.538492301539677"><title>Name: p_006s_06
Group: 1
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="0.059988002399521179" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="71.005798840231947"><title>Name: p_006s_14
Group: 1
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.059988002399521179" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="57.678464307138569"><title>Name: p_009s_13
Group: 1
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="2.3195360927814477" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="106.85862827434514"><title>Name: p_009s_15
Group: 1
Prop: 0.0232</title></rect><rect class="toyplot-Datum" height="0.0099980003999462497" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="148.5102979404119"><title>Name: p_011s_08
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="148.56028794241152"><title>Name: p_011s_10
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.059988002399506968" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="148.94021195760848"><title>Name: p_011s_11
Group: 1
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.069986002799424796" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="149.33013397320536"><title>Name: p_012s_06
Group: 1
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.18996200759846715" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="148.81023795240952"><title>Name: p_012s_12
Group: 1
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="0.17996400719856354" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="88.69226154769045"><title>Name: p_013s_02
Group: 1
Prop: 0.0018</title></rect><rect class="toyplot-Datum" height="32.203559288142372" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="98.960207958408319"><title>Name: p_013s_08
Group: 1
Prop: 0.3221</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="149.92001599680063"><title>Name: p_015s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="149.83003399320137"><title>Name: p_015s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="149.99000199960008"><title>Name: p_016s_01
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.079984003199371045" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="149.90001999600079"><title>Name: p_016s_02
Group: 1
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="149.98000399920016"><title>Name: p_016s_03
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.029994001199753484" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="149.97000599880025"><title>Name: p_016s_04
Group: 1
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="149.99000199960008"><title>Name: p_016s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.079984003199371045" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="149.83003399320137"><title>Name: p_016s_15
Group: 1
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.059988002399506968" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="132.33353329334133"><title>Name: p_017s_01
Group: 1
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="149.01019796040791"><title>Name: p_017s_03
Group: 1
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.17996400719854933" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="149.33013397320536"><title>Name: p_017s_05
Group: 1
Prop: 0.0018</title></rect><rect class="toyplot-Datum" height="98.640271945610891" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="51.359728054389116"><title>Name: p_026s_12
Group: 1
Prop: 0.9866</title></rect><rect class="toyplot-Datum" height="97.490501899620085" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="52.509498100379922"><title>Name: p_026s_12r
Group: 1
Prop: 0.9751</title></rect><rect class="toyplot-Datum" height="98.250349930014011" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="51.719656068786243"><title>Name: p_026s_14
Group: 1
Prop: 0.9827</title></rect><rect class="toyplot-Datum" height="72.795440911817622" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="71.905618876224764"><title>Name: p_027s_01
Group: 1
Prop: 0.7281</title></rect><rect class="toyplot-Datum" height="63.457308538292338" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="86.182763447310535"><title>Name: p_027s_03
Group: 1
Prop: 0.6347</title></rect><rect class="toyplot-Datum" height="61.787642471505691" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="86.562687462507498"><title>Name: p_027s_06
Group: 1
Prop: 0.618</title></rect><rect class="toyplot-Datum" height="0.1999600079984134" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="142.50149970005998"><title>Name: p_029s_10
Group: 1
Prop: 0.002</title></rect><rect class="toyplot-Datum" height="0.089982003599288873" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.90001999600079"><title>Name: p_031s_11
Group: 1
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="56.618676264747094" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="92.721455708858215"><title>Name: p_1027s_12r
Group: 1
Prop: 0.5663</title></rect><rect class="toyplot-Datum" height="68.556288742251539" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="80.813837232553496"><title>Name: p_1027s_20
Group: 1
Prop: 0.6857</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000199960008"><title>Name: p_001s_02
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099980003999462497" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.92001599680063"><title>Name: p_001s_03
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="1.5996800639871935" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="144.72105578884222"><title>Name: p_002.5s_01
Group: 2
Prop: 0.016</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="149.99000199960008"><title>Name: p_002.5s_04
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000199960008"><title>Name: p_002s_03
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000199960008"><title>Name: p_002s_06
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000199960008"><title>Name: p_002s_09
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 2
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.98000399920016"><title>Name: p_002s_14
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.98000399920016"><title>Name: p_002s_16
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="5.9188162367526473" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="71.885622875424914"><title>Name: p_004.5s_05
Group: 2
Prop: 0.0592</title></rect><rect class="toyplot-Datum" height="95.700859828034396" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="52.369526094781037"><title>Name: p_005s_06
Group: 2
Prop: 0.9572</title></rect><rect class="toyplot-Datum" height="92.651469706058805" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="52.539492101579675"><title>Name: p_005s_11
Group: 2
Prop: 0.9267</title></rect><rect class="toyplot-Datum" height="41.981603679264126" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="73.525294941011808"><title>Name: p_006s_01
Group: 2
Prop: 0.4199</title></rect><rect class="toyplot-Datum" height="6.6586682663467229" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.879824035192954"><title>Name: p_006s_06
Group: 2
Prop: 0.0666</title></rect><rect class="toyplot-Datum" height="17.576484703059378" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="53.429314137172568"><title>Name: p_006s_14
Group: 2
Prop: 0.1758</title></rect><rect class="toyplot-Datum" height="6.578684263147359" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="51.09978004399121"><title>Name: p_009s_13
Group: 2
Prop: 0.0658</title></rect><rect class="toyplot-Datum" height="34.353129374125189" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="72.505498900219948"><title>Name: p_009s_15
Group: 2
Prop: 0.3436</title></rect><rect class="toyplot-Datum" height="96.180763847230537" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="52.329534093181358"><title>Name: p_011s_08
Group: 2
Prop: 0.962</title></rect><rect class="toyplot-Datum" height="96.230753849230155" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="52.329534093181358"><title>Name: p_011s_10
Group: 2
Prop: 0.9625</title></rect><rect class="toyplot-Datum" height="96.75064987002601" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="52.189562087582473"><title>Name: p_011s_11
Group: 2
Prop: 0.9677</title></rect><rect class="toyplot-Datum" height="60.887822435512916" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="88.442311537692447"><title>Name: p_012s_06
Group: 2
Prop: 0.609</title></rect><rect class="toyplot-Datum" height="58.698260347930415" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="90.111977604479108"><title>Name: p_012s_12
Group: 2
Prop: 0.5871</title></rect><rect class="toyplot-Datum" height="6.2387522495500889" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="82.453509298140361"><title>Name: p_013s_02
Group: 2
Prop: 0.0624</title></rect><rect class="toyplot-Datum" height="4.3791241751649466" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="94.581083783243372"><title>Name: p_013s_08
Group: 2
Prop: 0.0438</title></rect><rect class="toyplot-Datum" height="3.2893421315736759" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="146.63067386522695"><title>Name: p_015s_13
Group: 2
Prop: 0.0329</title></rect><rect class="toyplot-Datum" height="2.8394321135772884" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="146.99060187962408"><title>Name: p_015s_14
Group: 2
Prop: 0.0284</title></rect><rect class="toyplot-Datum" height="2.1995600879824053" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="147.79044191161768"><title>Name: p_016s_01
Group: 2
Prop: 0.022</title></rect><rect class="toyplot-Datum" height="4.2991401719656039" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="145.60087982403519"><title>Name: p_016s_02
Group: 2
Prop: 0.043</title></rect><rect class="toyplot-Datum" height="2.2895420915816942" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="147.69046190761847"><title>Name: p_016s_03
Group: 2
Prop: 0.0229</title></rect><rect class="toyplot-Datum" height="2.8194361127774528" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="147.15056988602279"><title>Name: p_016s_04
Group: 2
Prop: 0.0282</title></rect><rect class="toyplot-Datum" height="2.2295540891821588" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="147.76044791041792"><title>Name: p_016s_06
Group: 2
Prop: 0.0223</title></rect><rect class="toyplot-Datum" height="2.9494101179764129" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="146.88062387522496"><title>Name: p_016s_15
Group: 2
Prop: 0.0295</title></rect><rect class="toyplot-Datum" height="5.2089582083583252" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="127.12457508498301"><title>Name: p_017s_01
Group: 2
Prop: 0.0521</title></rect><rect class="toyplot-Datum" height="7.1185762847430283" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="141.89162167566488"><title>Name: p_017s_03
Group: 2
Prop: 0.0712</title></rect><rect class="toyplot-Datum" height="7.5084983003399373" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="141.82163567286543"><title>Name: p_017s_05
Group: 2
Prop: 0.0751</title></rect><rect class="toyplot-Datum" height="1.2697460507898484" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.089982003599268"><title>Name: p_026s_12
Group: 2
Prop: 0.0127</title></rect><rect class="toyplot-Datum" height="2.3895220955808867" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.119976004799035"><title>Name: p_026s_12r
Group: 2
Prop: 0.0239</title></rect><rect class="toyplot-Datum" height="1.5996800639872077" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.119976004799035"><title>Name: p_026s_14
Group: 2
Prop: 0.016</title></rect><rect class="toyplot-Datum" height="5.6388722255549055" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="66.266746650669859"><title>Name: p_027s_01
Group: 2
Prop: 0.0564</title></rect><rect class="toyplot-Datum" height="4.2091581683663293" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="81.973605278944206"><title>Name: p_027s_03
Group: 2
Prop: 0.0421</title></rect><rect class="toyplot-Datum" height="4.449110177964414" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="82.113577284543084"><title>Name: p_027s_06
Group: 2
Prop: 0.0445</title></rect><rect class="toyplot-Datum" height="16.956608678264331" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="125.54489102179565"><title>Name: p_029s_10
Group: 2
Prop: 0.1696</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.89002199560088"><title>Name: p_031s_11
Group: 2
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="4.57908418316336" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="88.142371525694855"><title>Name: p_1027s_12r
Group: 2
Prop: 0.0458</title></rect><rect class="toyplot-Datum" height="4.4391121775644962" style="fill:rgb(55.3%,62.7%,79.6%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="76.374725054989"><title>Name: p_1027s_20
Group: 2
Prop: 0.0444</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.99000199960008"><title>Name: p_001s_02
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="149.92001599680063"><title>Name: p_001s_03
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.059988002399506968" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="144.66106778644271"><title>Name: p_002.5s_01
Group: 3
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="149.99000199960008"><title>Name: p_002.5s_04
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000199960008"><title>Name: p_002s_03
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000199960008"><title>Name: p_002s_06
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="150.0"><title>Name: p_002s_08
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000199960008"><title>Name: p_002s_09
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.98000399920016"><title>Name: p_002s_14
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.98000399920016"><title>Name: p_002s_16
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.12997400519896019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="71.755648870225954"><title>Name: p_004.5s_05
Group: 3
Prop: 0.0013</title></rect><rect class="toyplot-Datum" height="0.07998400319936394" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="52.289542091581673"><title>Name: p_005s_06
Group: 3
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.089982003599281768" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="52.449510097980394"><title>Name: p_005s_11
Group: 3
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.78984203159369315" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="72.735452909418115"><title>Name: p_006s_01
Group: 3
Prop: 0.0079</title></rect><rect class="toyplot-Datum" height="0.20995800839832413" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.669866026794629"><title>Name: p_006s_06
Group: 3
Prop: 0.0021</title></rect><rect class="toyplot-Datum" height="0.10997800439912453" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="53.319336132773444"><title>Name: p_006s_14
Group: 3
Prop: 0.0011</title></rect><rect class="toyplot-Datum" height="0.62987402519496527" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.469906018796244"><title>Name: p_009s_13
Group: 3
Prop: 0.0063</title></rect><rect class="toyplot-Datum" height="4.9290141971605692" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="67.576484703059378"><title>Name: p_009s_15
Group: 3
Prop: 0.0493</title></rect><rect class="toyplot-Datum" height="0.059988002399521179" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="52.269546090781837"><title>Name: p_011s_08
Group: 3
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.0099980003999107225" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="52.319536092781448"><title>Name: p_011s_10
Group: 3
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.31993601279743444" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="51.869626074785039"><title>Name: p_011s_11
Group: 3
Prop: 0.0032</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="88.442311537692447"><title>Name: p_012s_06
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.12997400519896019" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="89.982003599280148"><title>Name: p_012s_12
Group: 3
Prop: 0.0013</title></rect><rect class="toyplot-Datum" height="0.15996800639872788" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="82.293541291741633"><title>Name: p_013s_02
Group: 3
Prop: 0.0016</title></rect><rect class="toyplot-Datum" height="1.8396320735852925" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="92.741451709658079"><title>Name: p_013s_08
Group: 3
Prop: 0.0184</title></rect><rect class="toyplot-Datum" height="96.480703859228157" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.149970005998803"><title>Name: p_015s_13
Group: 3
Prop: 0.965</title></rect><rect class="toyplot-Datum" height="96.910617876424737" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.07998400319935"><title>Name: p_015s_14
Group: 3
Prop: 0.9693</title></rect><rect class="toyplot-Datum" height="97.750449910018006" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.039992001599678"><title>Name: p_016s_01
Group: 3
Prop: 0.9777</title></rect><rect class="toyplot-Datum" height="95.520895820835847" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.07998400319935"><title>Name: p_016s_02
Group: 3
Prop: 0.9554</title></rect><rect class="toyplot-Datum" height="97.500499900020003" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.189962007598474"><title>Name: p_016s_03
Group: 3
Prop: 0.9752</title></rect><rect class="toyplot-Datum" height="97.110577884423122" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.039992001599678"><title>Name: p_016s_04
Group: 3
Prop: 0.9713</title></rect><rect class="toyplot-Datum" height="97.73045390921817" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.029994001199753"><title>Name: p_016s_06
Group: 3
Prop: 0.9775</title></rect><rect class="toyplot-Datum" height="96.740651869626078" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.139972005598878"><title>Name: p_016s_15
Group: 3
Prop: 0.9676</title></rect><rect class="toyplot-Datum" height="64.257148570285949" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="62.867426514697058"><title>Name: p_017s_01
Group: 3
Prop: 0.6427</title></rect><rect class="toyplot-Datum" height="79.464107178564291" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="62.427514497100589"><title>Name: p_017s_03
Group: 3
Prop: 0.7948</title></rect><rect class="toyplot-Datum" height="72.275544891021795" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="69.546090781843631"><title>Name: p_017s_05
Group: 3
Prop: 0.7229</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.089982003599268"><title>Name: p_026s_12
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.089982003599281768" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.029994001199753"><title>Name: p_026s_12r
Group: 3
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.119976004799035"><title>Name: p_026s_14
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="13.737252549490094" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="52.529494101179765"><title>Name: p_027s_01
Group: 3
Prop: 0.1374</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="81.973605278944206"><title>Name: p_027s_03
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.079984003199342624" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="82.033593281343741"><title>Name: p_027s_06
Group: 3
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="2.6994601079784246" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="122.84543091381722"><title>Name: p_029s_10
Group: 3
Prop: 0.027</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.89002199560088"><title>Name: p_031s_11
Group: 3
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="18.046390721855616" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="70.09598080383924"><title>Name: p_1027s_12r
Group: 3
Prop: 0.1805</title></rect><rect class="toyplot-Datum" height="0.18996200759846715" style="fill:rgb(90.6%,54.1%,76.5%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="76.184763047390533"><title>Name: p_1027s_20
Group: 3
Prop: 0.0019</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="149.97000599880025"><title>Name: p_001s_02
Group: 4
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.92981403719255695" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="148.99020195960807"><title>Name: p_001s_03
Group: 4
Prop: 0.0093</title></rect><rect class="toyplot-Datum" height="1.5096980603879331" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="143.15136972605478"><title>Name: p_002.5s_01
Group: 4
Prop: 0.0151</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="149.99000199960008"><title>Name: p_002.5s_04
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="150.0"><title>Name: p_002.5s_07
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="150.0"><title>Name: p_002s_01
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="149.99000199960008"><title>Name: p_002s_03
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="150.0"><title>Name: p_002s_05
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="149.99000199960008"><title>Name: p_002s_06
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="149.98000399920016"><title>Name: p_002s_08
Group: 4
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="149.99000199960008"><title>Name: p_002s_09
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="150.0"><title>Name: p_002s_12
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="150.0"><title>Name: p_002s_13
Group: 4
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="149.97000599880025"><title>Name: p_002s_14
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.059988002399535389" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="149.92001599680063"><title>Name: p_002s_16
Group: 4
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="21.745650869826044" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.009998000399911"><title>Name: p_004.5s_05
Group: 4
Prop: 0.2175</title></rect><rect class="toyplot-Datum" height="2.2695460907818372" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.019996000799836"><title>Name: p_005s_06
Group: 4
Prop: 0.0227</title></rect><rect class="toyplot-Datum" height="2.4195160967806402" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.029994001199753"><title>Name: p_005s_11
Group: 4
Prop: 0.0242</title></rect><rect class="toyplot-Datum" height="22.415516896620666" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.319936012797449"><title>Name: p_006s_01
Group: 4
Prop: 0.2242</title></rect><rect class="toyplot-Datum" height="0.63987202559487599" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.029994001199753"><title>Name: p_006s_06
Group: 4
Prop: 0.0064</title></rect><rect class="toyplot-Datum" height="3.229354129174169" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.089982003599275"><title>Name: p_006s_14
Group: 4
Prop: 0.0323</title></rect><rect class="toyplot-Datum" height="0.44991001799639463" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.01999600079985"><title>Name: p_009s_13
Group: 4
Prop: 0.0045</title></rect><rect class="toyplot-Datum" height="17.376524695060979" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.199960007998399"><title>Name: p_009s_15
Group: 4
Prop: 0.1738</title></rect><rect class="toyplot-Datum" height="2.2595480903819265" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.009998000399911"><title>Name: p_011s_08
Group: 4
Prop: 0.0226</title></rect><rect class="toyplot-Datum" height="2.309538092381537" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009998000399911"><title>Name: p_011s_10
Group: 4
Prop: 0.0231</title></rect><rect class="toyplot-Datum" height="1.7796440711857713" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.089982003599268"><title>Name: p_011s_11
Group: 4
Prop: 0.0178</title></rect><rect class="toyplot-Datum" height="37.992401519696053" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.449910017996395"><title>Name: p_012s_06
Group: 4
Prop: 0.38</title></rect><rect class="toyplot-Datum" height="39.782043591281749" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.199960007998399"><title>Name: p_012s_12
Group: 4
Prop: 0.3979</title></rect><rect class="toyplot-Datum" height="32.083583283343337" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.209958008398296"><title>Name: p_013s_02
Group: 4
Prop: 0.3209</title></rect><rect class="toyplot-Datum" height="38.452309538092379" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="54.2891421715657"><title>Name: p_013s_08
Group: 4
Prop: 0.3846</title></rect><rect class="toyplot-Datum" height="0.11997600479904236" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.029994001199761"><title>Name: p_015s_13
Group: 4
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.069986002799439007" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009998000399911"><title>Name: p_015s_14
Group: 4
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.029994001199767695" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.009998000399911"><title>Name: p_016s_01
Group: 4
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.049990001999596245" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.029994001199753"><title>Name: p_016s_02
Group: 4
Prop: 0.0005</title></rect><rect class="toyplot-Datum" height="0.079984003199356835" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.109978004399117"><title>Name: p_016s_03
Group: 4
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.029994001199767695" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.009998000399911"><title>Name: p_016s_04
Group: 4
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.019996000799836"><title>Name: p_016s_06
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089982003599288873" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.049990001999589"><title>Name: p_016s_15
Group: 4
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="12.797440511897634" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.069986002799425"><title>Name: p_017s_01
Group: 4
Prop: 0.128</title></rect><rect class="toyplot-Datum" height="12.327534493101375" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.099980003999214"><title>Name: p_017s_03
Group: 4
Prop: 0.1233</title></rect><rect class="toyplot-Datum" height="19.426114777044596" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.119976004799035"><title>Name: p_017s_05
Group: 4
Prop: 0.1943</title></rect><rect class="toyplot-Datum" height="0.069986002799431901" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.019996000799836"><title>Name: p_026s_12
Group: 4
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.019996000799836"><title>Name: p_026s_12r
Group: 4
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.089982003599281768" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029994001199753"><title>Name: p_026s_14
Group: 4
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="2.509498100379929" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.019996000799836"><title>Name: p_027s_01
Group: 4
Prop: 0.0251</title></rect><rect class="toyplot-Datum" height="31.683663267346539" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.289942011597667"><title>Name: p_027s_03
Group: 4
Prop: 0.3169</title></rect><rect class="toyplot-Datum" height="31.973605278944213" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.059988002399528"><title>Name: p_027s_06
Group: 4
Prop: 0.3198</title></rect><rect class="toyplot-Datum" height="0.48990201959607305" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="122.35552889422115"><title>Name: p_029s_10
Group: 4
Prop: 0.0049</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="149.87002599480104"><title>Name: p_031s_11
Group: 4
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="19.976004799040204" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.119976004799035"><title>Name: p_1027s_12r
Group: 4
Prop: 0.1998</title></rect><rect class="toyplot-Datum" height="25.374925014997011" style="fill:rgb(65.1%,84.7%,32.9%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.809838032393522"><title>Name: p_1027s_20
Group: 4
Prop: 0.2538</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.960007998400329" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="50.0" y="50.009998000399911"><title>Name: p_001s_02
Group: 5
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="98.980203959208154" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="55.882352941176464" y="50.009998000399911"><title>Name: p_001s_03
Group: 5
Prop: 0.99</title></rect><rect class="toyplot-Datum" height="93.131373725254946" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="61.764705882352942" y="50.019996000799836"><title>Name: p_002.5s_01
Group: 5
Prop: 0.9315</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="67.647058823529406" y="50.009998000399911"><title>Name: p_002.5s_04
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="73.529411764705884" y="50.019996000799836"><title>Name: p_002.5s_07
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="79.411764705882348" y="50.019996000799836"><title>Name: p_002s_01
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764781" x="85.294117647058826" y="50.009998000399911"><title>Name: p_002s_03
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="91.176470588235304" y="50.019996000799836"><title>Name: p_002s_05
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.970005998800247" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="97.058823529411768" y="50.019996000799836"><title>Name: p_002s_06
Group: 5
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.960007998400329" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="102.94117647058823" y="50.019996000799836"><title>Name: p_002s_08
Group: 5
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.970005998800247" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="108.8235294117647" y="50.019996000799836"><title>Name: p_002s_09
Group: 5
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="114.70588235294119" y="50.019996000799836"><title>Name: p_002s_12
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.980003999200164" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="120.58823529411765" y="50.019996000799836"><title>Name: p_002s_13
Group: 5
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.960007998400329" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="126.47058823529412" y="50.009998000399911"><title>Name: p_002s_14
Group: 5
Prop: 0.9998</title></rect><rect class="toyplot-Datum" height="99.900019996000793" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="132.35294117647061" y="50.019996000799836"><title>Name: p_002s_16
Group: 5
Prop: 0.9992</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="138.23529411764704" y="50.009998000399911"><title>Name: p_004.5s_05
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="144.11764705882354" y="50.019996000799836"><title>Name: p_005s_06
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="150.0" y="50.029994001199753"><title>Name: p_005s_11
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.28994201159768807" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="155.88235294117646" y="50.029994001199761"><title>Name: p_006s_01
Group: 5
Prop: 0.0029</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="161.76470588235296" y="50.019996000799836"><title>Name: p_006s_06
Group: 5
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.069986002799439007" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="167.64705882352942" y="50.019996000799836"><title>Name: p_006s_14
Group: 5
Prop: 0.0007</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="173.52941176470588" y="50.01999600079985"><title>Name: p_009s_13
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.16996600679864571" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="179.41176470588238" y="50.029994001199753"><title>Name: p_009s_15
Group: 5
Prop: 0.0017</title></rect><rect class="toyplot-Datum" height="0.0099980003999107225" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="185.29411764705884" y="50.0"><title>Name: p_011s_08
Group: 5
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="191.1764705882353" y="50.009998000399911"><title>Name: p_011s_10
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.059988002399528284" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="197.05882352941177" y="50.029994001199739"><title>Name: p_011s_11
Group: 5
Prop: 0.0006</title></rect><rect class="toyplot-Datum" height="0.42991401719655897" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="202.9411764705882" y="50.019996000799836"><title>Name: p_012s_06
Group: 5
Prop: 0.0043</title></rect><rect class="toyplot-Datum" height="0.16996600679864571" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="208.8235294117647" y="50.029994001199753"><title>Name: p_012s_12
Group: 5
Prop: 0.0017</title></rect><rect class="toyplot-Datum" height="0.18996200759848136" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="214.70588235294122" y="50.019996000799814"><title>Name: p_013s_02
Group: 5
Prop: 0.0019</title></rect><rect class="toyplot-Datum" height="4.2591481703659255" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="220.58823529411762" y="50.029994001199775"><title>Name: p_013s_08
Group: 5
Prop: 0.0426</title></rect><rect class="toyplot-Datum" height="0.019996000799849867" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="226.47058823529412" y="50.009998000399911"><title>Name: p_015s_13
Group: 5
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="232.35294117647055" y="50.009998000399911"><title>Name: p_015s_14
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="238.23529411764704" y="50.009998000399911"><title>Name: p_016s_01
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.009998000399917828" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764355" x="244.11764705882354" y="50.019996000799836"><title>Name: p_016s_02
Group: 5
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099980003999206701" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764923" x="249.99999999999997" y="50.009998000399911"><title>Name: p_016s_03
Group: 5
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="255.88235294117646" y="50.009998000399911"><title>Name: p_016s_04
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="261.76470588235293" y="50.019996000799836"><title>Name: p_016s_06
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.019996000799835656" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="267.64705882352939" y="50.029994001199753"><title>Name: p_016s_15
Group: 5
Prop: 0.0002</title></rect><rect class="toyplot-Datum" height="0.04999000199958914" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="273.52941176470591" y="50.019996000799836"><title>Name: p_017s_01
Group: 5
Prop: 0.0005</title></rect><rect class="toyplot-Datum" height="0.07998400319936394" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="279.41176470588232" y="50.01999600079985"><title>Name: p_017s_03
Group: 5
Prop: 0.0008</title></rect><rect class="toyplot-Datum" height="0.099980003999199596" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="285.29411764705884" y="50.019996000799836"><title>Name: p_017s_05
Group: 5
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="291.1764705882353" y="50.019996000799836"><title>Name: p_026s_12
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="297.05882352941177" y="50.019996000799836"><title>Name: p_026s_12r
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="302.94117647058823" y="50.029994001199753"><title>Name: p_026s_14
Group: 5
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099980003999249334" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764071" x="308.82352941176475" y="50.009998000399911"><title>Name: p_027s_01
Group: 5
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.27994401119775603" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="314.70588235294116" y="50.009998000399911"><title>Name: p_027s_03
Group: 5
Prop: 0.0028</title></rect><rect class="toyplot-Datum" height="0.029994001199767695" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="320.58823529411762" y="50.029994001199761"><title>Name: p_027s_06
Group: 5
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="72.345530893821234" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411765207" x="326.47058823529409" y="50.009998000399911"><title>Name: p_029s_10
Group: 5
Prop: 0.7236</title></rect><rect class="toyplot-Datum" height="99.860027994401122" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="332.35294117647061" y="50.009998000399911"><title>Name: p_031s_11
Group: 5
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="0.099980003999199596" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="338.23529411764707" y="50.019996000799836"><title>Name: p_1027s_12r
Group: 5
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.77984403119376822" style="fill:rgb(100%,85.1%,18.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" width="5.8823529411764639" x="344.11764705882354" y="50.029994001199753"><title>Name: p_1027s_20
Group: 5
Prop: 0.0078</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="t140142f5c7bc44e39da518abe9fb3fb7" transform="translate(50.0,150.0)translate(0,10.0)"><line style="" x1="0" x2="300.0" y1="0" y2="0"></line><g><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(2.941176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">0</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(61.764705882352935,6)translate(0,7.5)"><tspan style="font-size:10.0px">10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(120.58823529411765,6)translate(0,7.5)"><tspan style="font-size:10.0px">20</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(179.41176470588235,6)translate(0,7.5)"><tspan style="font-size:10.0px">30</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(238.23529411764704,6)translate(0,7.5)"><tspan style="font-size:10.0px">40</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:middle" transform="translate(297.05882352941177,6)translate(0,7.5)"><tspan style="font-size:10.0px">50</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-3.0" y2="4.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-6"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1", "magnitude2", "magnitude3", "magnitude4", "magnitude5"], "id": "t689bd55618f94f4f965cd2ede1951a6e", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0001, 0.0, 0.0235, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0001, 0.7209, 0.0193, 0.0473, 0.1273, 0.9232, 0.7895, 0.9228, 0.4083, 0.0148, 0.0144, 0.01, 0.006, 0.01, 0.6114, 0.1884, 0.0008, 0.0017, 0.0001, 0.0002, 0.0002, 0.0, 0.0001, 0.0009, 0.1761, 0.0097, 0.0049, 0.0, 0.0, 0.0003, 0.053, 0.0036, 0.0165, 0.073, 0.0001, 0.0066, 0.0063], [0.0, 0.0007, 0.0133, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0012, 0.0, 0.0008, 0.2177, 0.0016, 0.0006, 0.0006, 0.0232, 0.0001, 0.0, 0.0006, 0.0007, 0.0019, 0.0018, 0.3221, 0.0, 0.0, 0.0, 0.0008, 0.0, 0.0003, 0.0, 0.0008, 0.0006, 0.0002, 0.0018, 0.9866, 0.9751, 0.9827, 0.7281, 0.6347, 0.618, 0.002, 0.0009, 0.5663, 0.6857], [0.0, 0.0001, 0.016, 0.0001, 0.0, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0001, 0.0, 0.0, 0.0001, 0.0001, 0.0592, 0.9572, 0.9267, 0.4199, 0.0666, 0.1758, 0.0658, 0.3436, 0.962, 0.9625, 0.9677, 0.609, 0.5871, 0.0624, 0.0438, 0.0329, 0.0284, 0.022, 0.043, 0.0229, 0.0282, 0.0223, 0.0295, 0.0521, 0.0712, 0.0751, 0.0127, 0.0239, 0.016, 0.0564, 0.0421, 0.0445, 0.1696, 0.0001, 0.0458, 0.0444], [0.0, 0.0, 0.0006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0013, 0.0008, 0.0009, 0.0079, 0.0021, 0.0011, 0.0063, 0.0493, 0.0006, 0.0001, 0.0032, 0.0, 0.0013, 0.0016, 0.0184, 0.965, 0.9693, 0.9777, 0.9554, 0.9752, 0.9713, 0.9775, 0.9676, 0.6427, 0.7948, 0.7229, 0.0, 0.0009, 0.0, 0.1374, 0.0, 0.0008, 0.027, 0.0, 0.1805, 0.0019], [0.0002, 0.0093, 0.0151, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0002, 0.0, 0.0, 0.0, 0.0001, 0.0006, 0.2175, 0.0227, 0.0242, 0.2242, 0.0064, 0.0323, 0.0045, 0.1738, 0.0226, 0.0231, 0.0178, 0.38, 0.3979, 0.3209, 0.3846, 0.0012, 0.0007, 0.0003, 0.0005, 0.0008, 0.0003, 0.0001, 0.0009, 0.128, 0.1233, 0.1943, 0.0007, 0.0001, 0.0009, 0.0251, 0.3169, 0.3198, 0.0049, 0.0002, 0.1998, 0.2538], [0.9998, 0.99, 0.9315, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9999, 0.9998, 0.9999, 1.0, 1.0, 0.9998, 0.9992, 0.0, 0.0, 0.0, 0.0029, 0.0001, 0.0007, 0.0, 0.0017, 0.0001, 0.0, 0.0006, 0.0043, 0.0017, 0.0019, 0.0426, 0.0002, 0.0, 0.0, 0.0001, 0.001, 0.0, 0.0, 0.0002, 0.0005, 0.0008, 0.001, 0.0, 0.0, 0.0, 0.0001, 0.0028, 0.0003, 0.7236, 0.9988, 0.001, 0.0078]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#t81bd05ec41654c65847f98e257b38e14 .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "t81bd05ec41654c65847f98e257b38e14";
            var axes = {"t140142f5c7bc44e39da518abe9fb3fb7": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 300.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>


### Make a slightly fancier plot and save to file


```python
## save plots for your favorite value of K
table = struct.get_clumpp_table(kpop=2)
#table = table.ix[myorder]
```

    mean scores across 20 replicates.



```python
## further styling of plot with css 
style = {"stroke":toyplot.color.near_black, 
         "stroke-width": 2}

## build barplot
canvas = toyplot.Canvas(width=600, height=250)
axes = canvas.cartesian(bounds=("5%", "95%", "5%", "45%"))
axes.bars(table, title=hover(table), style=style)

## add names to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "12px"}
axes.x.spine.style = style
axes.y.show = False
    
## options: uncomment to save plots. Only html retains hover.
import toyplot.svg
import toyplot.pdf
import toyplot.html
toyplot.svg.render(canvas, "struct.svg")
toyplot.pdf.render(canvas, "struct.pdf")
toyplot.html.render(canvas, "struct.html")

## show in notebook
canvas
```




<div align="center" class="toyplot" id="t44b7937ceec34cb0a4cb790ae0c90e1c"><svg class="toyplot-canvas-Canvas" height="250.0px" id="t93105a17dd334b268b1deab459bbf524" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" viewBox="0 0 600.0 250.0" width="600.0px" xmlns="http://www.w3.org/2000/svg" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink"><g class="toyplot-coordinates-Cartesian" id="ta87ae9e5527940ca904ffd01b99ee21a"><clipPath id="tc5854255c9f54280be5eeabd49e6756f"><rect height="120.0" width="560.0" x="20.0" y="2.5"></rect></clipPath><g clip-path="url(#tc5854255c9f54280be5eeabd49e6756f)"><g class="toyplot-mark-BarMagnitudes" id="t02127190be944094922ea043a8235794" style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2"><g class="toyplot-Series"><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117645" x="30.0" y="112.5"><title>Name: p_001s_02
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117645" x="40.588235294117645" y="112.5"><title>Name: p_001s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="5.9194080591940832" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="51.17647058823529" y="106.58059194080592"><title>Name: p_002.5s_01
Group: 0
Prop: 0.0592</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="61.764705882352942" y="112.5"><title>Name: p_002.5s_04
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="72.35294117647058" y="112.5"><title>Name: p_002.5s_07
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="82.941176470588232" y="112.5"><title>Name: p_002s_01
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="93.529411764705884" y="112.5"><title>Name: p_002s_03
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="104.11764705882354" y="112.5"><title>Name: p_002s_05
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="114.70588235294117" y="112.5"><title>Name: p_002s_06
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="125.29411764705883" y="112.5"><title>Name: p_002s_08
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="135.88235294117646" y="112.5"><title>Name: p_002s_09
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="146.47058823529412" y="112.5"><title>Name: p_002s_12
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="157.05882352941177" y="112.5"><title>Name: p_002s_13
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="167.64705882352939" y="112.5"><title>Name: p_002s_14
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="178.23529411764707" y="112.5"><title>Name: p_002s_16
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="188.82352941176472" y="12.509999000099992"><title>Name: p_004.5s_05
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="199.41176470588235" y="12.509999000099992"><title>Name: p_005s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="210.0" y="12.509999000099992"><title>Name: p_005s_11
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.890010998900109" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="220.58823529411765" y="12.609989001099885"><title>Name: p_006s_01
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="231.17647058823528" y="12.509999000099992"><title>Name: p_006s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.7900209979002" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="241.76470588235296" y="12.709979002099793"><title>Name: p_006s_14
Group: 0
Prop: 0.998</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="252.35294117647058" y="12.509999000099992"><title>Name: p_009s_13
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.870012998700133" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="262.94117647058823" y="12.62998700129987"><title>Name: p_009s_15
Group: 0
Prop: 0.9988</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="273.52941176470591" y="12.519998000199983"><title>Name: p_011s_08
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="284.11764705882354" y="12.509999000099992"><title>Name: p_011s_10
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.890010998900109" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="294.70588235294116" y="12.609989001099885"><title>Name: p_011s_11
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="95.100489951004903" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="305.29411764705878" y="17.399510048995104"><title>Name: p_012s_06
Group: 0
Prop: 0.9511</title></rect><rect class="toyplot-Datum" height="99.890010998900109" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="315.88235294117646" y="12.609989001099885"><title>Name: p_012s_12
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="94.570542945705427" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="326.47058823529414" y="17.929457054294573"><title>Name: p_013s_02
Group: 0
Prop: 0.9458</title></rect><rect class="toyplot-Datum" height="94.04059594040595" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="337.05882352941177" y="18.459404059594043"><title>Name: p_013s_08
Group: 0
Prop: 0.9405</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="347.64705882352945" y="12.509999000099992"><title>Name: p_015s_13
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="358.23529411764702" y="12.509999000099992"><title>Name: p_015s_14
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="368.8235294117647" y="12.509999000099992"><title>Name: p_016s_01
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="379.41176470588238" y="12.519998000199983"><title>Name: p_016s_02
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.890010998900109" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="390.0" y="12.609989001099885"><title>Name: p_016s_03
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="400.58823529411768" y="12.509999000099992"><title>Name: p_016s_04
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="411.1764705882353" y="12.509999000099992"><title>Name: p_016s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.98000199980001" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="421.76470588235293" y="12.519998000199983"><title>Name: p_016s_15
Group: 0
Prop: 0.9999</title></rect><rect class="toyplot-Datum" height="99.950004999500052" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="432.35294117647055" y="12.549995000499948"><title>Name: p_017s_01
Group: 0
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.90000999900009" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="442.94117647058823" y="12.599990000999904"><title>Name: p_017s_03
Group: 0
Prop: 0.9991</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="453.52941176470591" y="12.509999000099992"><title>Name: p_017s_05
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="464.11764705882354" y="12.509999000099992"><title>Name: p_026s_12
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="474.70588235294116" y="12.509999000099992"><title>Name: p_026s_12r
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117737" x="485.29411764705878" y="12.509999000099992"><title>Name: p_026s_14
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.950004999500052" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="495.88235294117652" y="12.549995000499948"><title>Name: p_027s_01
Group: 0
Prop: 0.9996</title></rect><rect class="toyplot-Datum" height="99.770022997700238" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117737" x="506.47058823529409" y="12.729977002299764"><title>Name: p_027s_03
Group: 0
Prop: 0.9978</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="517.05882352941182" y="12.509999000099992"><title>Name: p_027s_06
Group: 0
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="26.407359264073591" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="527.64705882352939" y="86.092640735926409"><title>Name: p_029s_10
Group: 0
Prop: 0.2641</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="538.23529411764707" y="112.5"><title>Name: p_031s_11
Group: 0
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="99.890010998900109" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="548.82352941176464" y="12.609989001099885"><title>Name: p_1027s_12r
Group: 0
Prop: 0.999</title></rect><rect class="toyplot-Datum" height="97.900209979002099" style="fill:rgb(40%,76.1%,64.7%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="559.41176470588232" y="14.599790020997904"><title>Name: p_1027s_20
Group: 0
Prop: 0.9791</title></rect></g><g class="toyplot-Series"><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117645" x="30.0" y="12.509999000099992"><title>Name: p_001s_02
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117645" x="40.588235294117645" y="12.509999000099992"><title>Name: p_001s_03
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="94.070592940705922" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="51.17647058823529" y="12.509999000099992"><title>Name: p_002.5s_01
Group: 1
Prop: 0.9408</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="61.764705882352942" y="12.509999000099992"><title>Name: p_002.5s_04
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="72.35294117647058" y="12.509999000099992"><title>Name: p_002.5s_07
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="82.941176470588232" y="12.509999000099992"><title>Name: p_002s_01
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="93.529411764705884" y="12.509999000099992"><title>Name: p_002s_03
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="104.11764705882354" y="12.509999000099992"><title>Name: p_002s_05
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="114.70588235294117" y="12.509999000099992"><title>Name: p_002s_06
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117638" x="125.29411764705883" y="12.509999000099992"><title>Name: p_002s_08
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="135.88235294117646" y="12.509999000099992"><title>Name: p_002s_09
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="146.47058823529412" y="12.509999000099992"><title>Name: p_002s_12
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="157.05882352941177" y="12.509999000099992"><title>Name: p_002s_13
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="167.64705882352939" y="12.509999000099992"><title>Name: p_002s_14
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="178.23529411764707" y="12.509999000099992"><title>Name: p_002s_16
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="188.82352941176472" y="12.509999000099992"><title>Name: p_004.5s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999916617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="199.41176470588235" y="12.5"><title>Name: p_005s_06
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="210.0" y="12.509999000099992"><title>Name: p_005s_11
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.099990000999893525" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="220.58823529411765" y="12.509999000099992"><title>Name: p_006s_01
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="231.17647058823528" y="12.509999000099992"><title>Name: p_006s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.19998000199980126" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="241.76470588235296" y="12.509999000099992"><title>Name: p_006s_14
Group: 1
Prop: 0.002</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117652" x="252.35294117647058" y="12.509999000099992"><title>Name: p_009s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.11998800119987862" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="262.94117647058823" y="12.509999000099992"><title>Name: p_009s_15
Group: 1
Prop: 0.0012</title></rect><rect class="toyplot-Datum" height="0.0099990000999916617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="273.52941176470591" y="12.509999000099992"><title>Name: p_011s_08
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="284.11764705882354" y="12.509999000099992"><title>Name: p_011s_10
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.099990000999893525" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="294.70588235294116" y="12.509999000099992"><title>Name: p_011s_11
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="4.8895110488951126" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="305.29411764705878" y="12.509999000099992"><title>Name: p_012s_06
Group: 1
Prop: 0.0489</title></rect><rect class="toyplot-Datum" height="0.099990000999893525" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="315.88235294117646" y="12.509999000099992"><title>Name: p_012s_12
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="5.4194580541945818" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="326.47058823529414" y="12.509999000099992"><title>Name: p_013s_02
Group: 1
Prop: 0.0542</title></rect><rect class="toyplot-Datum" height="5.9494050594940511" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="337.05882352941177" y="12.509999000099992"><title>Name: p_013s_08
Group: 1
Prop: 0.0595</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="347.64705882352945" y="12.509999000099992"><title>Name: p_015s_13
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="358.23529411764702" y="12.509999000099992"><title>Name: p_015s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="368.8235294117647" y="12.509999000099992"><title>Name: p_016s_01
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999916617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="379.41176470588238" y="12.509999000099992"><title>Name: p_016s_02
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.099990000999893525" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="390.0" y="12.509999000099992"><title>Name: p_016s_03
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="400.58823529411768" y="12.509999000099992"><title>Name: p_016s_04
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="411.1764705882353" y="12.509999000099992"><title>Name: p_016s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999916617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="421.76470588235293" y="12.509999000099992"><title>Name: p_016s_15
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.029997000299964327" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="432.35294117647055" y="12.519998000199983"><title>Name: p_017s_01
Group: 1
Prop: 0.0003</title></rect><rect class="toyplot-Datum" height="0.089991000899912521" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="442.94117647058823" y="12.509999000099992"><title>Name: p_017s_03
Group: 1
Prop: 0.0009</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="453.52941176470591" y="12.509999000099992"><title>Name: p_017s_05
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="464.11764705882354" y="12.509999000099992"><title>Name: p_026s_12
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.0099990000999916617" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117624" x="474.70588235294116" y="12.5"><title>Name: p_026s_12r
Group: 1
Prop: 0.0001</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117737" x="485.29411764705878" y="12.509999000099992"><title>Name: p_026s_14
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="0.039996000399955989" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="495.88235294117652" y="12.509999000099992"><title>Name: p_027s_01
Group: 1
Prop: 0.0004</title></rect><rect class="toyplot-Datum" height="0.21997800219977215" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117737" x="506.47058823529409" y="12.509999000099992"><title>Name: p_027s_03
Group: 1
Prop: 0.0022</title></rect><rect class="toyplot-Datum" height="0.0" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="517.05882352941182" y="12.509999000099992"><title>Name: p_027s_06
Group: 1
Prop: 0.0</title></rect><rect class="toyplot-Datum" height="73.582641735826414" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="527.64705882352939" y="12.509999000099992"><title>Name: p_029s_10
Group: 1
Prop: 0.7359</title></rect><rect class="toyplot-Datum" height="99.990000999900005" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.588235294117567" x="538.23529411764707" y="12.509999000099992"><title>Name: p_031s_11
Group: 1
Prop: 1.0</title></rect><rect class="toyplot-Datum" height="0.099990000999893525" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="548.82352941176464" y="12.509999000099992"><title>Name: p_1027s_12r
Group: 1
Prop: 0.001</title></rect><rect class="toyplot-Datum" height="2.0897910208979127" style="fill:rgb(98.8%,55.3%,38.4%);fill-opacity:1.0;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" width="10.58823529411768" x="559.41176470588232" y="12.509999000099992"><title>Name: p_1027s_20
Group: 1
Prop: 0.0209</title></rect></g></g></g><g class="toyplot-coordinates-Axis" id="td49238a6a1f6418ab3cd52961cc4189f" transform="translate(30.0,112.5)translate(0,10.0)"><line style="stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:2" x1="0" x2="540.0" y1="0" y2="0"></line><g><line style="" x1="5.294117647058823" x2="5.294117647058823" y1="0" y2="-5"></line><line style="" x1="15.882352941176471" x2="15.882352941176471" y1="0" y2="-5"></line><line style="" x1="26.470588235294116" x2="26.470588235294116" y1="0" y2="-5"></line><line style="" x1="37.05882352941177" x2="37.05882352941177" y1="0" y2="-5"></line><line style="" x1="47.64705882352941" x2="47.64705882352941" y1="0" y2="-5"></line><line style="" x1="58.235294117647065" x2="58.235294117647065" y1="0" y2="-5"></line><line style="" x1="68.8235294117647" x2="68.8235294117647" y1="0" y2="-5"></line><line style="" x1="79.41176470588236" x2="79.41176470588236" y1="0" y2="-5"></line><line style="" x1="90.0" x2="90.0" y1="0" y2="-5"></line><line style="" x1="100.58823529411765" x2="100.58823529411765" y1="0" y2="-5"></line><line style="" x1="111.17647058823529" x2="111.17647058823529" y1="0" y2="-5"></line><line style="" x1="121.76470588235294" x2="121.76470588235294" y1="0" y2="-5"></line><line style="" x1="132.35294117647058" x2="132.35294117647058" y1="0" y2="-5"></line><line style="" x1="142.94117647058823" x2="142.94117647058823" y1="0" y2="-5"></line><line style="" x1="153.52941176470588" x2="153.52941176470588" y1="0" y2="-5"></line><line style="" x1="164.1176470588235" x2="164.1176470588235" y1="0" y2="-5"></line><line style="" x1="174.7058823529412" x2="174.7058823529412" y1="0" y2="-5"></line><line style="" x1="185.2941176470588" x2="185.2941176470588" y1="0" y2="-5"></line><line style="" x1="195.88235294117646" x2="195.88235294117646" y1="0" y2="-5"></line><line style="" x1="206.47058823529412" x2="206.47058823529412" y1="0" y2="-5"></line><line style="" x1="217.05882352941177" x2="217.05882352941177" y1="0" y2="-5"></line><line style="" x1="227.64705882352942" x2="227.64705882352942" y1="0" y2="-5"></line><line style="" x1="238.23529411764704" x2="238.23529411764704" y1="0" y2="-5"></line><line style="" x1="248.8235294117647" x2="248.8235294117647" y1="0" y2="-5"></line><line style="" x1="259.4117647058824" x2="259.4117647058824" y1="0" y2="-5"></line><line style="" x1="270.0" x2="270.0" y1="0" y2="-5"></line><line style="" x1="280.5882352941177" x2="280.5882352941177" y1="0" y2="-5"></line><line style="" x1="291.1764705882353" x2="291.1764705882353" y1="0" y2="-5"></line><line style="" x1="301.7647058823529" x2="301.7647058823529" y1="0" y2="-5"></line><line style="" x1="312.3529411764706" x2="312.3529411764706" y1="0" y2="-5"></line><line style="" x1="322.94117647058823" x2="322.94117647058823" y1="0" y2="-5"></line><line style="" x1="333.5294117647059" x2="333.5294117647059" y1="0" y2="-5"></line><line style="" x1="344.1176470588235" x2="344.1176470588235" y1="0" y2="-5"></line><line style="" x1="354.70588235294116" x2="354.70588235294116" y1="0" y2="-5"></line><line style="" x1="365.29411764705884" x2="365.29411764705884" y1="0" y2="-5"></line><line style="" x1="375.88235294117646" x2="375.88235294117646" y1="0" y2="-5"></line><line style="" x1="386.47058823529414" x2="386.47058823529414" y1="0" y2="-5"></line><line style="" x1="397.05882352941177" x2="397.05882352941177" y1="0" y2="-5"></line><line style="" x1="407.6470588235294" x2="407.6470588235294" y1="0" y2="-5"></line><line style="" x1="418.2352941176471" x2="418.2352941176471" y1="0" y2="-5"></line><line style="" x1="428.8235294117647" x2="428.8235294117647" y1="0" y2="-5"></line><line style="" x1="439.4117647058824" x2="439.4117647058824" y1="0" y2="-5"></line><line style="" x1="450.0" x2="450.0" y1="0" y2="-5"></line><line style="" x1="460.5882352941176" x2="460.5882352941176" y1="0" y2="-5"></line><line style="" x1="471.1764705882353" x2="471.1764705882353" y1="0" y2="-5"></line><line style="" x1="481.764705882353" x2="481.764705882353" y1="0" y2="-5"></line><line style="" x1="492.35294117647055" x2="492.35294117647055" y1="0" y2="-5"></line><line style="" x1="502.94117647058823" x2="502.94117647058823" y1="0" y2="-5"></line><line style="" x1="513.5294117647059" x2="513.5294117647059" y1="0" y2="-5"></line><line style="" x1="524.1176470588235" x2="524.1176470588235" y1="0" y2="-5"></line><line style="" x1="534.7058823529412" x2="534.7058823529412" y1="0" y2="-5"></line></g><g><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(5.294117647058823,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_001s_02</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(15.882352941176471,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_001s_03</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(26.470588235294116,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002.5s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(37.05882352941177,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002.5s_04</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(47.64705882352941,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002.5s_07</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(58.235294117647065,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(68.8235294117647,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_03</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(79.41176470588236,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_05</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(90.0,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(100.58823529411765,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_08</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(111.17647058823529,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_09</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(121.76470588235294,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_12</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(132.35294117647058,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_13</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(142.94117647058823,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_14</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(153.52941176470588,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_002s_16</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(164.1176470588235,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_004.5s_05</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(174.7058823529412,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_005s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(185.2941176470588,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_005s_11</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(195.88235294117646,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_006s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(206.47058823529412,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_006s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(217.05882352941177,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_006s_14</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(227.64705882352942,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_009s_13</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(238.23529411764704,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_009s_15</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(248.8235294117647,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_011s_08</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(259.4117647058824,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_011s_10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(270.0,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_011s_11</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(280.5882352941177,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_012s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(291.1764705882353,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_012s_12</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(301.7647058823529,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_013s_02</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(312.3529411764706,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_013s_08</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(322.94117647058823,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_015s_13</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(333.5294117647059,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_015s_14</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(344.1176470588235,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(354.70588235294116,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_02</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(365.29411764705884,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_03</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(375.88235294117646,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_04</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(386.47058823529414,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(397.05882352941177,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_016s_15</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(407.6470588235294,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_017s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(418.2352941176471,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_017s_03</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(428.8235294117647,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_017s_05</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(439.4117647058824,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_026s_12</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(450.0,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_026s_12r</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(460.5882352941176,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_026s_14</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(471.1764705882353,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_027s_01</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(481.764705882353,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_027s_03</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(492.35294117647055,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_027s_06</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(502.94117647058823,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_029s_10</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(513.5294117647059,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_031s_11</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(524.1176470588235,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_1027s_12r</tspan></text><text style="font-weight:normal;stroke:none;text-anchor:start" transform="translate(534.7058823529412,10.0)rotate(60)translate(0,4.5)"><tspan style="font-size:12.0px">p_1027s_20</tspan></text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0" x1="0" x2="0" y1="-5.0" y2="7.5"></line><text style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle" x="0" y="-10.0"></text></g></g></g></svg><div class="toyplot-interactive"><ul class="toyplot-mark-popup" onmouseleave="this.style.visibility='hidden'" style="background:rgba(0%,0%,0%,0.75);border:0;border-radius:6px;color:white;cursor:default;list-style:none;margin:0;padding:5px;position:fixed;visibility:hidden">
            <li class="toyplot-mark-popup-title" style="color:lightgray;cursor:default;padding:5px;list-style:none;margin:0"></li>
            <li class="toyplot-mark-popup-save-csv" onmouseout="this.style.color='white';this.style.background='steelblue'" onmouseover="this.style.color='steelblue';this.style.background='white'" style="border-radius:3px;padding:5px;list-style:none;margin:0">
                Save as .csv
            </li>
        </ul><script>
        (function()
        {
          var data_tables = [{"title": "Bar Data", "names": ["left", "right", "baseline", "magnitude0", "magnitude1"], "id": "t02127190be944094922ea043a8235794", "columns": [[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5], [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0592, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.999, 1.0, 0.998, 1.0, 0.9988, 0.9999, 1.0, 0.999, 0.9511, 0.999, 0.9458, 0.9405, 1.0, 1.0, 1.0, 0.9999, 0.999, 1.0, 1.0, 0.9999, 0.9996, 0.9991, 1.0, 1.0, 1.0, 1.0, 0.9996, 0.9978, 1.0, 0.2641, 0.0, 0.999, 0.9791], [1.0, 1.0, 0.9408, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0001, 0.0, 0.001, 0.0, 0.002, 0.0, 0.0012, 0.0001, 0.0, 0.001, 0.0489, 0.001, 0.0542, 0.0595, 0.0, 0.0, 0.0, 0.0001, 0.001, 0.0, 0.0, 0.0001, 0.0003, 0.0009, 0.0, 0.0, 0.0001, 0.0, 0.0004, 0.0022, 0.0, 0.7359, 1.0, 0.001, 0.0209]], "filename": "toyplot"}];

          function save_csv(data_table)
          {
            var uri = "data:text/csv;charset=utf-8,";
            uri += data_table.names.join(",") + "\n";
            for(var i = 0; i != data_table.columns[0].length; ++i)
            {
              for(var j = 0; j != data_table.columns.length; ++j)
              {
                if(j)
                  uri += ",";
                uri += data_table.columns[j][i];
              }
              uri += "\n";
            }
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = data_table.filename + ".csv";

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
          }

          function open_popup(data_table)
          {
            return function(e)
            {
              var popup = document.querySelector("#t44b7937ceec34cb0a4cb790ae0c90e1c .toyplot-mark-popup");
              popup.querySelector(".toyplot-mark-popup-title").innerHTML = data_table.title;
              popup.querySelector(".toyplot-mark-popup-save-csv").onclick = function() { popup.style.visibility = "hidden"; save_csv(data_table); }
              popup.style.left = (e.clientX - 50) + "px";
              popup.style.top = (e.clientY - 20) + "px";
              popup.style.visibility = "visible";
              e.stopPropagation();
              e.preventDefault();
            }

          }

          for(var i = 0; i != data_tables.length; ++i)
          {
            var data_table = data_tables[i];
            var event_target = document.querySelector("#" + data_table.id);
            event_target.oncontextmenu = open_popup(data_table);
          }
        })();
        </script><script>
        (function()
        {
            function _sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function _mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function _log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function _in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(_in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return _mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return _sign(segment.domain.min) * Math.pow(base, _mix(_log(segment.domain.min, base), _log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            function display_coordinates(e)
            {
                var current = svg.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            var root_id = "t44b7937ceec34cb0a4cb790ae0c90e1c";
            var axes = {"td49238a6a1f6418ab3cd52961cc4189f": [{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 50.5, "min": -0.5}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 540.0, "min": 0.0}, "scale": "linear"}]};

            var svg = document.querySelector("#" + root_id + " svg");
            svg.addEventListener("click", display_coordinates);
        })();
        </script></div></div>



### Calculating the best K 
I haven't gotten around to writing the code for this yet (contributors are welcome!). For now, I like using the site http://taylor0.biology.ucla.edu/structureHarvester/. It's great. Super easy. Zip up all the files in our structure directory, submit them to the site, and you're done. 


```python
%%bash -s "$STRUCTDIR"

## creates zip dir of all files ending with _f
zip structure-files-K2and3.zip *_f
```

### Copying this notebook to your computer/cluster
You can easily copy this notebook and then just replace my file names with your filenames to run your analysis. Just click on the [Download Notebook] link at the top of this page. Then run `jupyter-notebook` from a terminal and open this notebook from the dashboard.
