# estructural-equation-models

The first part of this script estimate the scale of effect (SoF) of three landscape metrics: % of
foret cover, number of patches, and edge density for 169 primate communities in the neotropical region, using
the multifit function (Huais, 2018). The landscape metrics was extracted from 11 different buffers (1000m to 5000m
with 400 m interval from each distance), so the SoF for all three metrics was assessed based on these 11 buffers.

The second part of this script evaluate the potential spatial autocorrelation presented in the data, by testing five
different spatial correlation structures (exponential, gaussian, linear, rational quadratic, and spheric).

The third and last part of this script (depicted in the image below) I implement the landscape metrics in their respectve scales of effect in structural
equation models to investigate the direct and indirect effects of habitat loss and fragmentation on the species
richness and functional diversity of primate communities in the neotropics.

![alt text](https://github.com/Lucas-a-pereira/estructural-equation-models/blob/main/SEM_diagram_neotropics.png?raw=true)
