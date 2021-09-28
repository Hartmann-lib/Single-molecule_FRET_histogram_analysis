# Single-molecule_FRET_histogram_analysis

One widely used method of single-molecule FRET spectroscopy is the acquisition of intensity time traces from immobilized molecules using TIRF microscopy. This widefield technique allows a quick recording of hundreds to thousands of time traces in parallel with a time resolution of a few milliseconds. In such experiments, it is common to collect FRET efficiency histograms using all recorded time points of selected regions. This can lead to a relatively fast accumulation of smooth distributions. However, in the case where the time average of the system has not quite reached the ensemble average, molecules with long time traces are prone to over representation. This in turn can lead to wrong conclusions. In the attached graphical user interface (FitXGauss), fractions of FRET populations can be quantified by up to five Gaussian distributions. To this end, histogram contributions of individual molecules are weighted with respect to their time trace length. Boot strapping is performed to evaluate the confidence interval of each derived FRET population. This approach presents a more robust way of FRET efficiency histogram analysis.


![Figure_gui](https://user-images.githubusercontent.com/58071484/135050815-3ce29f00-b619-403c-bea0-17c839b0c554.jpg)


![Figure_output](https://user-images.githubusercontent.com/58071484/135051420-a5df01b6-8cde-41d0-b93b-55f2df0093c7.jpg)
