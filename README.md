[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Integrated Spatial Diagnostics for Clinically Deployable, Quantitative Tissue Pathology
<br>
Jia-Yun Chen, Yu-An Chen, Jia-Ren Lin, Yilin Xu1, Claire T. Avillach, Clemens B. Hug, Crystal Chiu1, Sabrina Chan, Terri Woo, Helmut G. Rennke, Joseph V. Bonventre, Astrid Weins, Peter K. Sorger, Sandro Santagata<br>
<br>
<br>
<img src="./docs/KidneyDiagnosticManuscript-landingPage.png)" style="max-width:500px;width:100%">

## SUMMARY
Advances in spatial multi-omic tissue atlases have revealed the molecular and architectural complexity of human disease, yet these insights remain largely inaccessible in routine diagnostics. Here, we introduce a clinically compatible spatial pathology framework that unifies one-shot multiplex immunofluorescence, multimodal computational pathology, and interactive digital interpretation within the operational constraints of standard histopathology. From a single biopsy section, this approach generates quantitative, compartment-resolved molecular maps and multimodal predictions that remain aligned with established diagnostic conventions and turnaround times. Applied to 199 consecutively collected kidney biopsies – a testbed rich in mechanistically informative biomarkers and diagnostic complexity – the framework identified disease-defining molecular signatures; quantified key processes such as complement activation, paraprotein imbalance, tubular injury, mesangial expansion, and fibrosis; and uncovered continuous spatial trajectories of progressive podocyte loss and glomerular remodeling. Integrating H&E-derived morphological embeddings with spatial molecular features markedly improved diagnostic accuracy and mechanistic interpretation compared with either modality alone. These findings demonstrate a clinically deployable architecture that incorporates rich molecular and spatial information into routine diagnostic workflows while preserving the interpretive logic of pathology. Although evaluated in kidney disease, the modular design supports rapid extension to other tissues, providing a generalizable route toward quantitative, predictive, and mechanism-informed tissue pathology.


## ACCESS THE DATA
Example images can be accessed through Harvard tissue Atlas website (https://s3.amazonaws.com/www.cycif.org/152-kidney-imaging/LSP20571/index.html). Cell count tables have been released via Github (https://github.com/labsyspharm/Chen-Santagata-Kidney-2026/).
Additional imagesa and data will be provide upon request.
<br>

## Codes for imaging data processing
High-plex whole-slide images were acquired using tissue cyclic immunofluorescence (t-CyCIF; Lin et al, 2018) and then stitched and registered using ASHLAR (https://github.com/labsyspharm/ashlar). 
