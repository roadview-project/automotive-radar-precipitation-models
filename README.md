# automotive-radar-precipitation-models
Models and analysis for understanding and estimating the distribution of radar backscatter of precipitation over the radar cube/detection space typical for automotive radars
# file description
The project contains two Matlab files that should be possible to run as is.
1. precipitationRadarBackscatterDistribution.m
Contains analysis of unit volume RCS, effect of range on detection bin volume, and effect of velocity on bin volume/probability of finding rain drop in each velocity bin. It also conatains a simulation part where drops are genereated according to a uniform random distribution in space and randomly assigned a drop diameter according to the Marshall-Palmer drop distribution. The resulting radar backscatter from the simulation is then compared to the analytically derived values.
2. points_on_sphere.m
Contains analysis on the effect of angular position on detection bin volume/probability of finding a drop in some specific angular bin/region. It also contains analysis showing the pdf of one single drop over velocity bins.

# Publications 
These models were used in the work to be presented as "Rain Reflectivity Distribution and Detectability in Automotive Radar" at the 2025 IEEE 28th International Conference on Intelligent Transportation Systems (ITSC 2025)

# Acknowledgment
Co-Funded by the European Union (grant no. 101069576). Views and opinions expressed are however those
of the author(s) only and do not necessarily reflect those of the European Union or the European Climate,
Infrastructure and Environment Executive Agency (CINEA). Neither the European Union nor the granting
authority can be held responsible for them. UK and Swiss participants in this project are supported by Innovate
UK (contract no. 10045139) and the Swiss State Secretariat for Education, Research and Innovation (contract
no. 22.00123) respectively
