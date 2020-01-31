# Economics of Connectomics

Code for the paper, [*Efficient Coding in the Economics of Human Brain Connectomics*](https://www.biorxiv.org/content/10.1101/2020.01.14.906842v2).

## Setup

```
platform       x86_64-pc-linux-gnu         
arch           x86_64                      
os             linux-gnu                   
system         x86_64, linux-gnu           

language       R                           
version        R version 3.5.1 (2018-07-02)

MATLAB:        8.4.0.150421 (R2014b)
```

### Functions

Download the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/)

Download [polyfix](https://www.mathworks.com/matlabcentral/fileexchange/54207-polyfix-x-y-n-xfix-yfix-xder-dydx)

Download [spin test](https://github.com/frantisekvasa/rotate_parcellation)

Download [GenLouvain](http://netwiki.amath.unc.edu/GenLouvain)

Move to /scripts/

## Author

Dale Zhou (dalezhou [at] pennmedicine.upenn.edu)

## Project Organization

```

    ├── data                                    <- Data goes here.
    ├── figures                                 <- Figures for manuscript.
    │   ├── fig1                   
    │   ├── fig2
    │   ├── fig3
    │   ├── fig4
    │   ├── fig5
    │   ├── fig6
    │   ├── fig7
    │   ├── fig8
    │   ├── supplement
    │ 
    ├── scripts                                 <- Downloaded functions go here
    │
    ├── fig1.R                                  <- CBF vs CMRGlu.
    ├── fig2.R                                  <- CBF of paths.
    ├── fig3.R                                  <- CBF/network structure trade-offs.
    ├── fig4_randomVsBrainRD.m                  <- Rate-distortion function.
    ├── fig4.R                                  <- Compression eff vs age*sex.
    ├── fig5_resourceEfficiencyConcatenate.R    <- Send/rcv compression eff.
    ├── fig5_fitRD.m                            <- Individual compression eff.
    ├── fig5.R                                  <- Random walk vs chemotaxis.
    ├── fig6.R                                  <- High & low fidelity regimes.
    ├── fig7.R                                  <- Send/rcv vs scaling, myelin.
    ├── fig7_atlasCorrelations.m                <- Spin test for spatial corr.
    ├── fig8_preprocessRegional.m               <- Regional send/recv compression eff.
    ├── fig8.R                                  <- Rich-club hubs & speed-accuracy.
    ├── resource_efficiency_wei.m               <- Resource efficiency for WU networks.
    ├── supplement_richClubCBF.R                <- CBF of rich-club vs non-rich-club.
    ├── supplement_subnetworks.R                <- Yeo 7 subnetwork levels.
    │
    ├── README.md

```

## Order of scripts

1. fig1.R
2. fig2.R
3. fig3.R
4. fig4.R
5. fig4_randomVsBrainRD.m
6. fig5_resourceEfficiencyConcatenate.R
7. fig5_fitRD.m
8. fig5.R
9. fig6.R
10. fig7.R
11. fig7_atlasCorrelations.m
12. fig8_preprocessRegional.m
13. fig8.R
14. supplement_richClubCBF.R 
15. supplement_subnetworks.R

### Notes

Will need cluster for extra memory for some of these steps, such as Fig 2.

qlogin -q qlogin.himem.q -l h_vmem=15G,s_vmem=15G