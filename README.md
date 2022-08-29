# Statistical inference of the rates of cell proliferation and phenotypic switching in cancer

This repository holds MATLAB codes which implement the estimation framework from the paper "Statistical inference of the rates of cell proliferation and phenotypic switching in cancer" by Gunnarsson, Leder and Foo.

The repository includes background codes and six example script files which demonstrate how the framework can be used. To run the script files, they must be saved in the same folder as the background codes.

The script files are the following:

1. "Script_number" implements estimation for cell number data under the statistical model (5) from Section 4.2 of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term.

2. "Script_number_reducible" implements estimation for cell number data under the statistical model (5) from Section 4.2 of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term, for a model with reducible switching dynamics. The model structure is displayed in Figure 10 of Appendix B of the accompanying paper.

3. "Script_fractional_Model_I" implements estimation for cell fraction data under Model I from Section 8 of the accompanying paper.

4. "Script_fractional_Model_II" implements estimation for cell fraction data under Model II from Section 8 of the accompanying paper.

5. "Script_number_generation_estimation" generates parameter regimes and artifical cell number data as described in Appendix D of the accompanying paper. It also estimates parameters from the artificial data under the statistical model (5) from Section 4.2 of the accompanying paper with E_(i,l)^(num) = 0, i.e. no measurement error term.

6. "Script_fractional_generation_estimation" generates parameter regimes and artifical cell fraction data as described in Appendix D of the accompanying paper. It also estimates parameters from the artificial data under the statistical model (11) from Section 5.2 of the accompanying paper with E_(i,l)^(frac) = 0, i.e. no measurement error term.
