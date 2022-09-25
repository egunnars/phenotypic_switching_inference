# Statistical inference of the rates of cell proliferation and phenotypic switching in cancer

This repository holds MATLAB codes which implement the estimation framework from the paper "Statistical inference of the rates of cell proliferation and phenotypic switching in cancer" by Gunnarsson, Leder and Foo.

The repository includes background codes and six example script files which demonstrate how the framework can be used. To run the script files, they must be saved in the same folder as the background codes.

The script files are the following:

1.	"Script_number" implements estimation for cell number data under the statistical model from Section “Estimation for cell number data” of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term.

2.	"Script_number_reducible" implements estimation for cell number data under the statistical model from Section “Estimation for cell number data” of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term, for a model with reducible switching dynamics. The model structure is displayed in Figure 10 of Appendix “Estimation for reducible switching dynamics”.

3.	"Script_fractional_Model_I" implements estimation for cell fraction data under Model I from Section “Application: Transition between stem and non-stem cell states in SW620 colon cancer” of the accompanying paper.

4.	"Script_fractional_Model_II" implements estimation for cell fraction data under Model II from Section “Application: Transition between stem and non-stem cell states in SW620 colon cancer”  of the accompanying paper.

5.	"Script_fractional_Model_Ia" implements estimation for cell fraction data under Model Ia from Section “Application: Transition between stem and non-stem cell states in SW620 colon cancer” of the accompanying paper.

6.	"Script_fractional_Model_IIa" implements estimation for cell fraction data under Model IIa from Section “Application: Transition between stem and non-stem cell states in SW620 colon cancer”  of the accompanying paper.

7.	"Script_number_generation_estimation" generates parameter regimes and artifical cell number data as described in Appendix “Generation of artificial data” of the accompanying paper. It also estimates parameters from the artificial data under the statistical model from section “Estimation for cell number data” of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term.

8.	"Script_fractional_generation_estimation" generates parameter regimes and artifical cell fraction data as described in Appendix “Generation of artificial data”  of the accompanying paper. It also estimates parameters from the artificial data under the statistical model from Section “Estimation for cell fraction data” of the accompanying paper with E_(i,l)^(frac) = 0, i.e. no measurement error term.
