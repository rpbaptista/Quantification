Author: Renata Porciuncula Baptista
e-mail: renata.porciunculabaptista@cea.fr / re.porci.rb@gmail.com
Date: 24/07/2019 - 14h30

------------------- BASICS QUANTIFICATION SCRIPTS ----------------

Order scripts: 
- computing_applying_tsc_model.py
	TO be modified: FA, TR, input file
	Parameters to be measured verified: factor_correction_liquid, T1 tubes
	output: model

- generate CSF, WM, GM mask
	Tested:
		- direct from sodium : no need to worry about PSF, nor order of axis images, low quality
	To be verified:
		generate from anatomy : need to use mask lr(low resolution), higher quality
		Consequential changes, use the lr mask, well chose resolution in voxels 5mm/(resolution proton)

- correctingTSC.py
	receives the model and the masks, and compute TSC, factor to be better 
	output: input_TSC.nii, sodium aparent
	

