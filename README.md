# `eht-mcp` 

*Developer:* Joseph Farah, SAO (joseph.farah@cfa.harvard.edu)

**What is it?** `eht-mcp` is a complete **model comparison pipeline** with a focus on ring-like structures.

  - **Modular**: don't like the included ring metrics? Add your own!
  - **Fast**: it only takes a second to parameterize a dozen image attributes!
  - **Powerful**: built on `PontiFEX`, `Hyperion` and `eht-imaging`, this pipeline comes equipped with a robust and versatile backend custom-made for EHT images and model fits!


## Quick Links
 * [`eht-mcp` at `master` branch](https://github.com/sao-eht/model-comparison-pipeline)
 * [`eht-imaging`](https://github.com/achael/eht-imaging), which form the file IO backbone of the pipeline
 * [`PontiFEX` and `Hyperion`](https://github.com/sao-eht/pontifex), which form the feature extraction backbone of the pipeline
 * [Challenge ID database](https://docs.google.com/spreadsheets/d/11oCD7T6okr3iRfJjvlFesowY3O8QLhv87HbOpP79v3c/edit?usp=sharing), used for naming your files!
 * [Frequently updated documentation of metrics used.](https://github.com/sao-eht/eht-mcp/blob/master/docs/eht-mcp_metric_descriptions.pdf)
## How do I submit an image to the pipeline?
Each challenge will have its own ID, folder and output within the pipeline. To submit your image for comparison and analysis...

  1) Ensure that your image file types meet the **submission guidelines** (see below).
  
  2) Properly **name your file(s)** (see below).

  3) In the MCFE WG Dropbox, navigate to `Results/eht-mcp/submissions/<challenge_name><challenge-ID>/`. 

  4) Copy and paste your ***images only*** into the folder.

### Submission guidelines

The pipeline analyzes images only in order to provide a uniform understanding of how different models compare. Therefore, we request that only images be submitted to the submission folders. The images must comply with the following requirements:

  1) The **file type** must be `.fits`. 

  2) The FITS header must include the **size of a pixel in the x direction** in microarcseconds, stored in `CDELT1` keyword, and the **size of a pixel in the y direction** in microarcseconds, stored in `CDELT2` keyword. 

  3) It would be very helpful if some **padding of the field of view** is included, to facilitate the estimation of some ring metrics (especially sharpness, orientation, and flux distribution).

### Submission name guidelines

All files must use the below name convention. The general format is:

```
./LastnameFirstnameExtraInfo_YYMMDD_ID.fits
```

 * `Lastname`, `Firstname` (self explanatory) ex: `JohnsonMichael`
 * `ExtraInfo` -- submitting multiple model fits? distinguish them here. ex: `LL` or `HOPS`
 * `YYMMDD`, date of SUBMISSION (e.g. November 8th, 2018 becomes `181108`)
 * `ID`, dataset ID unique to each challenge. For the first challenge, this tag is: `WCROMH`. See: https://docs.google.com/spreadsheets/d/11oCD7T6okr3iRfJjvlFesowY3O8QLhv87HbOpP79v3c/edit?usp=sharing

Example: Michael Johnson submits two model fits for the sample challenge, one using the `HOPS` pipeline for data and one using the `CASA` pipeline, on Nov 8th, 2018. His resulting files will look like:

```
./JohnsonMichaelHOPS_181108_WCROMH.fits
./JohnsonMichaelCASA_181108_WCROMH.fits
``` 

## How do I see the results of the pipeline?

Work in progress. You can always check the `Jupyter` notebook that is the frontend for the latest comparison graphs for any particular challenge. **Coming soon**: output PDF containing all the graphs and comparisons, as well as diagnostic info for the fits files!


## How can I contribute to the pipeline?

The pipeline has a very specific structure designed for front-end ease of use, backend versatility, and long-term sustainability. The pipeline supports a huge smorgasbord of metrics efficiently, and is designed to remember how to do everything we need while allowing us to be selective about what we use in our final comparisons. If you would like to contribute a feature extraction metric to the pipeline, please follow these steps:

  1) **Fork** the pipeline locally. 

  2) Open a **new branch** with a descriptive name. 

  3) Write your metric--**outside** of the pipeline. The metirc should meet the following requirements:

  * The metric should be written as a single main function, with **small** helper functions implemented as needed. 
  * The metric should take in only `self` as an input, and use `self.image` from the `PontiFEX` `ImageAnalysisObject` for the actual analysis. 
  * The function should return a single object (preferably), easily compared and ready to be graphed. 
  * A time complexity of less than N^2 is preferred. 

4) Once the metric has been functionalized as described in 3), it's time to add it to the pipeline! The files you will need to change are `Pontifex.py` and `eht-mcp.ipynb`. Make the following changes to `Pontifex.py`:

* Add any imports to the top of the file.
* Find the two lines: 

```
### function attributes specific to the model comparison pipeline ###
### end function attributes specific to the model comparison pipeline ###
```
* Place your metric function between those two lines. Make it the last function in the list!

Next, make the following changes to `eht-mcp.ipynb`:

* Add the name of your metric output as a key to the `self.analysisResults` dictionary within the `Submission` object. Use camelCase for the key name, and make the default value something unphysical, or zero. 
* Define a function below using the **same name as the dictionary key**, and pass it **only the argument `self`**.
* Call the function you just added to `PontiFEX` from the function you just added to the `Submission` object, using `ImageAnalysisObject`. Assign the output to the key you created in the `self.analysisResults` dictionary. 

5) Test the pipeline locally using `sample_data/` and make sure it works, at least on the base case. 

6) Open a pull request on Github from your fork, and request a branch merge with master.


