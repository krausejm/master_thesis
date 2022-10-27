# Determination of the beam asymmetry $\Sigma$ in $\eta$ and $\eta'$ photoproduction using Bayesian statistics
## Master thesis at Helmholtz-Institute (HISKP), University of Bonn
This repository includes analysis scripts as well as the TeX-sources for the thesis and talks. A brief table of contents is given below. Complementary analysis scripts are only available internally for HISKP members. The main goal of the thesis was to test different regression approaches to determine observables in particular hadronic photoproduction reactions. An introduction into the physical concepts is given in the [thesis](/TeX/master_thesis_krause_jakob_3134010.pdf)

* `DPG2022` - Subscription and Talk held at the 2022 DPG conference
* `RooFit` - Scripts and plots used to perform an unbinned fit on selected data in $\eta'$ photoproduction, using Root
* `TeX` - All TeX sources for talks and the final thesis. The template used to create the thesis is [ubonnthesis](https://www.pi.uni-bonn.de/lehre/uni-bonn-thesis).
* `bayes` - Summarizes all Bayesian fits that were performed during the thesis using [Stan](https://mc-stan.org)
    - `etap_event_based_fit` - Event based fit using data from $\eta'$ photoproduction
    - `event_based_fit` - Event based fit using data from $\eta$ photoproduction (provided by [Farah Afzal](https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8064))
    - `realdeal` - Binned fit using data from $\eta$ photoproduction (provided by [Farah Afzal](https://bonndoc.ulb.uni-bonn.de/xmlui/handle/20.500.11811/8064))
    - `toyMC` - Toy data sets with known parameters used to test the employed methods
  
* `demonstration` - Short snippets of code used for demonstration in the thesis
* `etap_PWA` - PWA predictions for the polarization observable $\Sigma$ in $\eta'$ photoproduction
* `figs` - Plots and graphs generated during the analysis of data acquired at CBELSA/TAPS for $\eta'$ photoproduction
* `prev_results` - Results for the polarization observable $\Sigma$ in $\eta'$ photoproduction from previous measurements at GRAAL and CLAS
* `ressources` - Literature that was used during the course of the thesis

