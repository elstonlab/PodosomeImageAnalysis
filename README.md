# Podosome Image Analysis 

Accompanies Herron JC\*, Hu S\*, Watanabe T\*, Nogueira AT, Aaron J, Taylor A, Pablo M, Chew TL, Elston TC, Hahn KM. <i>Nano-architecture of phagocytic podosomes</i>. (\*Equal contribution)


__Includes code to:__
- Locate podosomes and phagocytosis sites using a persistent homology-based pipeline ('Persistent_Homology_Example.ipynb', 'helper_functions.py').
	- ExampleData contains the files needed to demo.
- Validate this pipeline ('PH_Pipeline_Validation_Analysis.ipynb').
- Analyze podosomes from iPALM data and quantify relevant features ('iPALM_Podosome_Analysis.ipynb').
	- Demo included ('iPALM_Podosome_Analysis_Demo.ipynb') using Z-stack data located in ExampleData.
- Analyze podosomes and relevant species from TIRF-SIM and 3D-SIM data ('perform_pod_species_analysis.py', 'podosome_analysis_other_species.py', 'TIRFSIM_OtherSpeciesAnalysis.ipynb', '3DSIM_Analysis_Pax.ipynb', 'helper_functions.py')
	- Note that myosin and a-actinin were analyzed similar to '3DSIM_Analysis_Pax.ipynb', but this was not included due to redundancy.
- Perform a 3D contour rendering of the average podosome ('PreprocessPodForMayavi.ipynb', 'rendering_figure.py')


__Installation:__
Note: This code was created using Python v 3.7.10 and Jupyter Notebook 6.1.4
1. Requires Python (https://www.python.org/downloads/) and Jupyter Notebook (https://jupyter.org/install).
2. Requires installation of the 'unusual' Python packages Dionysus2 (https://mrzv.org/software/dionysus2/, v2.0.8, __required__), and statannot (https://pypi.org/project/statannot/, used for plotting statistical annotations). 

All analyses conducted on a Macintosh computer (2.3 GHz Quad-Core Intel Core i7) and should only require a maximum of a few minutes per step.




