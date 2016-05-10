# Phospho-Enrichment-Correction
This small Java applet can be used to clean the output of Thermo Proteome Discoverer™ Mass Spectrometry software for PhosphoRS 3.0 that produces incorrect peptide sequences given their Phospho-enrichment site probabilities.
## How To Use
1. Export the peptide sequences with their probabilities of Phospho-enrichment (including any other data that will be used thereafter) as comma separated values series.
2. Import into the Phospho-Enrichment-Correction program and select the necessary parameters;
  - Remove duplicate scans remove redundant data from mass spec scans
  - Remove any unused sequences (definded by the data sheet produced by the mass spec software)
  - Delete duplicate scans to remove duplicate peptide sequences with the same probabilities
3. The applet then normalises, cleans and outputs the data to two separate files of corrected peptide sequences (one under and one over the specified probability threshold) in the CSV format that it was originally entered

![Screenshot](http://i.imgur.com/lsHnZ3M.png)
