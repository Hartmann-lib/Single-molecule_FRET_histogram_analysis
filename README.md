# Single-molecule_FRET_histogram_analysis

One widely used method of single-molecule FRET spectroscopy is the acquisition of intensity time traces from immobilized molecules using TIRF microscopy. This widefield technique allows a quick recording of hundreds to thousands of time traces in parallel with a time resolution of a few milliseconds. In such experiments, it is common to collect FRET efficiency histograms using all recorded time points of selected regions. This can lead to a relatively fast accumulation of smooth distributions. However, in the case where the time average of the system has not quite reached the ensemble average, molecules with long time traces are prone to over representation. This in turn can lead to wrong conclusions. 

### Effect of molecule weighting on the FRET efficiency histograms:

![Figure_1_FitXGauss](https://user-images.githubusercontent.com/58071484/135234365-af626e6a-85c0-4bc9-979d-4bf981e35b55.png)

In the attached graphical user interface (FitXGauss), fractions of FRET populations can be quantified by up to five Gaussian distributions. To this end, histogram contributions of individual molecules are weighted with respect to their time trace length. Boot strapping is performed to evaluate the confidence interval of each derived FRET population. This approach presents a more robust way of FRET efficiency histogram analysis.

### Graphical user interface:

![Figure_2_FitXGauss](https://user-images.githubusercontent.com/58071484/135234611-f09c2dde-92e2-4bd2-9f6e-5ab6ff4729da.png)

### Provided representations:

![Figure_3_FitXGauss](https://user-images.githubusercontent.com/58071484/135234646-ef5d0198-15ac-408a-862b-432869a4b12c.png)

### Data import

Data are imported from multiple ascii files, where each file contains the FRET efficiency histogram of a single molecule. The format within the files should contain two columns: the edges of the FRET efficiency histogram and the counts per bin, both separated by a delimiter. For details see 'example' folder.

### Data export

The main figure containing the fit result can be exported via File -> Export Graph (Ctrl+E). Under the menu point 'View' detailed information about the individual histograms, the overall confidence interval and the population fractions can be exported in separated figures. Beside that, fitting results can be saved as mat-file using the third button (feather).
