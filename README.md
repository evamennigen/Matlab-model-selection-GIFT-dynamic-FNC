# Matlab-model-selection-GIFT-dynamic-FNC
Two Matlab script to a) perform backward model selection and univariate tests on the selected models and b) to create result summaries.
These scripts need the GIFT toolbox from http://mialab.mrn.org/software/gift/index.html
GIFT offers the same steps for static FNC but not for dynamic FNC states.
These scripts perform model selection and univariate tests on each dynamic FNC state separately; subjects mean frame-wise displacement is added to the univariate models to further control for motion-related artifacts (there are other motion correction processes in earlier steps of the pipeline).
