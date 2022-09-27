Supporting data and code for the publication:
# Unraveling the Palindromic and Non-Palindromic Motifs of Retroviral Integration Sites by Statistical Mixture Models

The study-supporting data are separated into three main parts:
- **Data**
- **Code**
- **Results**

The GitHub ***IS_Motifs***  repository contains only the '**Code**' directory. The '**Data**' and the '**Results**' directories are beeing saved elsewhere. The three parts, however, make a complete study data set. Therefore the whole data set is described in this file.

The '**Data**' directory includes some of the coordinates of the integration sites that were not previously publicly available. Sources for additional publicky available data sets are named in the '**IS**' directory as well as in the publication.
  - **IS** - the directory where information describing inttegration sites is stored. The directory also contains other sample-specific directories where sequences of ranges around the integration sites are stored. When running the code, coordinates of the integration sites should be stored in sample-specific directories.
  - **Annotations** - other-than-IS data needed for results reproduction. The data stored here relates to the Alu sequences and "HIV hotspot" part of the study.

The '**Code**' directory contains code used for processing the data and production of the figures presented in the publication. Further information can be found in the readme file of the '**Code**' directory. Sub-directories include:
  - **01_PreMix**
  - **02_Mix**
  - **03_PostMix**
  - **04_HotSpot**
  - **05_Figures**
  
The directories are ordered as they were run during the study. Firs three folders (**01_PreMix**, **02_Mix**, **03_PostMix**) contain the code related to the analysis of the mixture models. Directory **04_HotSpot** contains the code related to the analysis of the HIV-1 integration hotspot. The **05_Figures** directory contains the code that produces the figures present in the publication. Furter detailed information about the content of the directories can be found in the directory README_code.md file.

The '**Results**' directory contains the data that were produced in the study. These data can be produced by the code in th '**Code**' directory and are needed for the production of the figures. This directory does NOT contain all the data produced by the code - however, the intermediate data can be produced by the the code in the '**Code**' directory. 
  - **Figures**
  - **HotSpot**
  - **Mixture**

The '**Mixtures**' directory contains data produced by the EM algorithm. The data produced during the analysis of the HIV-1 intra-Alu integration hotspot are stored in the '**HotSpot**' directory. The '**Figures**' folder contains only the figures that cannot be reproduced by any of the code present in the '**Code**' directory. The data in the 'Results' folder should be sufficient to reproduce the Figures present in the publication.
