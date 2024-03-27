# pyINETA

pyINETA is a Python package for analyzing INADEQUATE NMR spectra.

pyINETA can perform basic tasks such as reading and referencing (using basic shifting) the INADEQUATE spectra and peak picking. It is designed to filter these picked peaks to identify networks of peaks (ideally) coming from the same compound which it then matches to a simulated INADEQUATE database of metabolites to identify metabolites present in the query INADEQUATE spectra.


## INSTALLATION INSTRUCTIONS:

1. Clone this repository locally:

`git clone https://github.com/rhilt13/PyINETA.git`

2. Install the required packages.
Preferably, if you use anaconda as a package manager, you can create a new environment with the required packages using the provided environment.yml file.

`conda env create -f <path_to_pyineta_repo>/environment.yml`

3. Activate the created environment before running the scripts.

`conda activate pyineta`

4. Run the run_pyineta.py script with the required options.

`python <path_to_pyineta_repo>/run_pyineta.py -c config.ini <other options>`

For help and options, please run:

```
python <path_to_pyineta_repo>/run_pyineta.py -h

usage: run_pyineta.py [-h] -c CONFIGFILE [-o OUTDIR] [-s STEPS] [-n NET]
                      [-d DBNAME] [-f FIGURE]

Script to run the INETA pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIGFILE, --configfile CONFIGFILE
                        Required: A config file with all the options and
                        parameters required for the INETA run.
  -o OUTDIR, --outdir OUTDIR
                        Optional: Full path to the output folder.(Default:
                        Current folder)
  -s STEPS, --steps STEPS
                        Optional: Specify which steps you want to run. Can be
                        one of {all,load,pick,cluster,find,match,plot,summary,
                        singleplot,load+,pick+,cluster+,find+,match+}. Adding
                        a + to the end of option runs all steps after the
                        specified step.
  -n NET, --net NET     Required with -s singlePlot: Specify which Network you
                        want to plot.
  -d DBNAME, --dbname DBNAME
                        Required with -s singlePlot: Specify which database
                        metabolite you want to plot.
  -f FIGURE, --figure FIGURE
                        Optional: Generate figures- yes or no. Default: Yes
```

## EXAMPLE RUN:

Inside the example/ folder, 2 examples are provided.

### `example1/`

To run this example,

`cd <path_to_pyineta_repo>/example/example1/`
`python ../../run_pyineta.py -c config.ini -o example1_output`

All outputs will be written to a folder called `example1_output`.
You can cross-check the results with the provided output in folder `example1_default_output`.

### `example2/`

Since the input and temporary files are larger than 100 Mb, please contact the Edison lab for getting the raw files to try out example2, if needed.
The output files and a config file is provided for reference.



