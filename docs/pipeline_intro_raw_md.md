<h1><a id="ehtmcp_0"></a><code>eht-mcp</code></h1>
<p><em>Developer:</em> Joseph Farah, SAO (<a href="mailto:joseph.farah@cfa.harvard.edu">joseph.farah@cfa.harvard.edu</a>)</p>
<p><strong>What is it?</strong> <code>eht-mcp</code> is a complete <strong>model comparison pipeline</strong> with a focus on ring-like structures.</p>
<ul>
<li><strong>Modular</strong>: don’t like the included ring metrics? Add your own!</li>
<li><strong>Fast</strong>: it only takes a second to parameterize a dozen image attributes!</li>
<li><strong>Powerful</strong>: built on <code>PontiFEX</code>, <code>Hyperion</code> and <code>eht-imaging</code>, this pipeline comes equipped with a robust and versatile backend custom-made for EHT images and model fits!</li>
</ul>
<h2><a id="Quick_Links_11"></a>Quick Links</h2>
<ul>
<li><a href="https://github.com/sao-eht/model-comparison-pipeline"><code>eht-mcp</code> at <code>master</code> branch</a></li>
<li><a href="https://github.com/achael/eht-imaging"><code>eht-imaging</code></a>, which form the file IO backbone of the pipeline</li>
<li><a href="https://github.com/sao-eht/pontifex"><code>PontiFEX</code> and <code>Hyperion</code></a>, which form the feature extraction backbone of the pipeline</li>
<li><a href="https://docs.google.com/spreadsheets/d/11oCD7T6okr3iRfJjvlFesowY3O8QLhv87HbOpP79v3c/edit?usp=sharing">Challenge ID database</a>, used for naming your files!</li>
</ul>
<h2><a id="How_do_I_submit_an_image_to_the_pipeline_17"></a>How do I submit an image to the pipeline?</h2>
<p>Each challenge will have its own ID, folder and output within the pipeline. To submit your image for comparison and analysis…</p>
<ol>
<li>
<p>Ensure that your image file types meet the <strong>submission guidelines</strong> (see below).</p>
</li>
<li>
<p>Properly <strong>name your file(s)</strong> (see below).</p>
</li>
<li>
<p>In the MCFE WG Dropbox, navigate to <code>Results/eht-mcp/submissions/&lt;challenge_name&gt;&lt;challenge-ID&gt;/</code>.</p>
</li>
<li>
<p>Copy and paste your <strong><em>images only</em></strong> into the folder.</p>
</li>
</ol>
<h3><a id="Submission_guidelines_28"></a>Submission guidelines</h3>
<p>The pipeline analyzes images only in order to provide a uniform understanding of how different models compare. Therefore, we request that only images be submitted to the submission folders. The images must comply with the following requirements:</p>
<ol>
<li>
<p>The <strong>file type</strong> must be <code>.fits</code>.</p>
</li>
<li>
<p>The FITS header must include the <strong>size of a pixel in the x direction</strong> in microarcseconds, stored in <code>CDELT1</code> keyword, and the <strong>size of a pixel in the y direction</strong> in microarcseconds, stored in <code>CDELT2</code> keyword.</p>
</li>
<li>
<p>It would be very helpful if some <strong>padding of the field of view</strong> is included, to facilitate the estimation of some ring metrics (especially sharpness, orientation, and flux distribution).</p>
</li>
</ol>
<h3><a id="Submission_name_guidelines_38"></a>Submission name guidelines</h3>
<p>All files must use the below name convention. The general format is:</p>
<pre><code>./LastnameFirstnameExtraInfo_YYMMDD_ID.fits
</code></pre>
<ul>
<li><code>Lastname</code>, <code>Firstname</code> (self explanatory) ex: <code>JohnsonMichael</code></li>
<li><code>ExtraInfo</code> – submitting multiple model fits? distinguish them here. ex: <code>LL</code> or <code>HOPS</code></li>
<li><code>YYMMDD</code>, date of SUBMISSION (e.g. November 8th, 2018 becomes <code>181108</code>)</li>
<li><code>ID</code>, dataset ID unique to each challenge. For the first challenge, this tag is: <code>WCROMH</code>. See: <a href="https://docs.google.com/spreadsheets/d/11oCD7T6okr3iRfJjvlFesowY3O8QLhv87HbOpP79v3c/edit?usp=sharing">https://docs.google.com/spreadsheets/d/11oCD7T6okr3iRfJjvlFesowY3O8QLhv87HbOpP79v3c/edit?usp=sharing</a></li>
</ul>
<p>Example: Michael Johnson submits two model fits for the sample challenge, one using the <code>HOPS</code> pipeline for data and one using the <code>CASA</code> pipeline, on Nov 8th, 2018. His resulting files will look like:</p>
<pre><code>./JohnsonMichaelHOPS_181108_WCROMH.fits
./JohnsonMichaelCASA_181108_WCROMH.fits
</code></pre>
<h2><a id="How_do_I_see_the_results_of_the_pipeline_58"></a>How do I see the results of the pipeline?</h2>
<p>Work in progress. You can always check the <code>Jupyter</code> notebook that is the frontend for the latest comparison graphs for any particular challenge. <strong>Coming soon</strong>: output PDF containing all the graphs and comparisons, as well as diagnostic info for the fits files!</p>
<h2><a id="How_can_I_contribute_to_the_pipeline_63"></a>How can I contribute to the pipeline?</h2>
<p>The pipeline has a very specific structure designed for front-end ease of use, backend versatility, and long-term sustainability. The pipeline supports a huge smorgasbord of metrics efficiently, and is designed to remember how to do everything we need while allowing us to be selective about what we use in our final comparisons. If you would like to contribute a feature extraction metric to the pipeline, please follow these steps:</p>
<ol>
<li>
<p><strong>Fork</strong> the pipeline locally.</p>
</li>
<li>
<p>Open a <strong>new branch</strong> with a descriptive name.</p>
</li>
<li>
<p>Write your metric–<strong>outside</strong> of the pipeline. The metirc should meet the following requirements:</p>
</li>
</ol>
<ul>
<li>The metric should be written as a single main function, with <strong>small</strong> helper functions implemented as needed.</li>
<li>The metric should take in only <code>self</code> as an input, and use <code>self.image</code> from the <code>PontiFEX</code> <code>ImageAnalysisObject</code> for the actual analysis.</li>
<li>The function should return a single object (preferably), easily compared and ready to be graphed.</li>
<li>A time complexity of less than N^2 is preferred.</li>
</ul>
<ol start="4">
<li>Once the metric has been functionalized as described in 3), it’s time to add it to the pipeline! The files you will need to change are <code>Pontifex.py</code> and <code>eht-mcp.ipynb</code>. Make the following changes to <code>Pontifex.py</code>:</li>
</ol>
<ul>
<li>Add any imports to the top of the file.</li>
<li>Find the two lines:</li>
</ul>
<pre><code>### function attributes specific to the model comparison pipeline ###
### end function attributes specific to the model comparison pipeline ###
</code></pre>
<ul>
<li>Place your metric function between those two lines. Make it the last function in the list!</li>
</ul>
<p>Next, make the following changes to <code>eht-mcp.ipynb</code>:</p>
<ul>
<li>Add the name of your metric output as a key to the <code>self.analysisResults</code> dictionary within the <code>Submission</code> object. Use camelCase for the key name, and make the default value something unphysical, or zero.</li>
<li>Define a function below using the <strong>same name as the dictionary key</strong>, and pass it <strong>only the argument <code>self</code></strong>.</li>
<li>Call the function you just added to <code>PontiFEX</code> from the function you just added to the <code>Submission</code> object, using <code>ImageAnalysisObject</code>. Assign the output to the key you created in the <code>self.analysisResults</code> dictionary.</li>
</ul>
<ol start="5">
<li>
<p>Test the pipeline locally using <code>sample_data/</code> and make sure it works, at least on the base case.</p>
</li>
<li>
<p>Open a pull request on Github from your fork, and request a branch merge with master.</p>
</li>
</ol>