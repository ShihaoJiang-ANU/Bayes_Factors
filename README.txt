
Prerequisites:

Numpy/Scipy/Matplotlib (generally installed)

Dynesty (Installable from pip: >pip install dynesty)

The basic steps for doing the analysis of a single data set are as follows:

1. Extracting/grouping data from Excel into a single file

> python3 getdata.py -c SiO2 -o sio2_basin.txt

There are a couple of options to do with the source of the data:
 -G|--global-plate
 -P|--petro
 -L|--seismic

Many of your columns have Unicode characters in the names, eg for
lambda 1 etc. You will need to copy and paste from the Excel column
heading, eg open excel, select the column heading you want, <copy>, go
to command line and do python3 get_data.py -c <paste> ... This worked
fine for me.
 
Check the first few lines of the __init__ method in platedata.py to see what these mean in
terms of the files that are loaded, in short:
 - Seismic means that the x coordinate is loaded from the x_mean column of OIB_compiled_location_x_seis.xlsx
 - Global Plate mean that the x coordinate is loaded from global_mean_plate column instead of x_mean_plate
 - Petrolog means that data is read from OIB_geochemistry_petrolog_final.xlsx

2. Run an analysis

> python3 evaluate_single_hierarchical.py -i sio2_basin.txt -y 35 -Y 55 -m 0

This will run a quick inversion using a constant model. The y range of the data needs to be configured
for the posterior histogram that is used to gauge the degree of uncertainty (the blue shading in the
plots).

The -m parameter is for the type of model: 0 = constant, 1 = linear, 2 = bilinear

By default the results will simply be plotted on the screen. To save them you specify an
output prefix with the -o <prefix> option. This will output a number of npy files with the
given prefix recording all the information relevant to the analysis.

Also by default the Nested sampler only runs a small number of slices to be quick (50). For a
proper analysis this should be increased to 500 .. 1000. The option to do this -N 500 for example.

See the Inversion Targets section in the Makefile for some examples of this step.

In general you will need to set for each dataset:
 -y <ymin> -Y <ymax> 
 -m <model type: 0, 1, 2>
 -c <change point mean> -C <change point std> these default to 60 and 5
 -N <sampling amount >= 500 for final analysis, 50 for quick run>

The difference from the previous results is that I have changed the default standard deviation
for the prior on the location of the change point to 5 (previously 20). This removes the issue with
change points at around 120km deep as shown yesterday.

3. Plotting

Assuming you have run the above with an "-o test" option added, ie
> python3 evaluate_single_hierarchical.py -i sio2_basin.txt -y 35 -Y 55 -m 0 -o test

You can then plot with 
> python3 plot.py -i test -d sio2_basin.txt -l "SiO\$_2\$"

Note that you can use latex commands in the yaxis label, but the $'s have to be escaped (with \ directly
in the shell or with $$ in the Makefile)

If you use the bilinear model

> python3 evaluate_single_hierarchical.py -i sio2_basin.txt -y 35 -Y 55 -m 2 -o test_bilinear

You can then plot with the -b switch to also plot the histogram of the location of the
change or break point.

> python3 plot.py -i test_bilinear -d sio2_basin.txt -l "SiO\$_2\$" -b

Lastly, the Makefile contains all the steps to generate the Figures for the
report I sent yesterday. If you just do:

> make

And then

> pdflatex analysis_notes.tex
> bibtex analysis_notes
> pdflatex analysis_notes.tex

You should get the same report I sent previously (with slightly different results for
SiO2 due to the change of prior).

