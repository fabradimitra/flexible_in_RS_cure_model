# flexible_in_RS_cure_model
We relaxed the parametric assumption of the relation between cure fraction and covariates in the context of relative survival Weibull cure model, by using Generalized Additive Models (GAM) and one-hidden layer Neural Networks. This repository contains three simulation study used to assess the performance of the proposed models against the standard model which uses a Generalized Linear Model (GLM) instead.

#### (i) GLM model: Weibull mixture cure model with GLM in the cure fraction.
#### (ii) GAM model: Weibull mixture cure model with GAM in the cure fraction.
#### (iii) Nnet_x model: Weibull mixture cure model with Nnet in the cure fraction with x-sized hidden layer.

In scenario 1 we specify the cure fraction as main effects and interaction of three covariates. In scenario 2 we extend scenario 1 using seven covariates. In scenario 3 we consider a only main effects of seven covariates. 

Each model subfolder contains a folder "data", where the results of the simulation are saved when the simulation is run; a folder "functions", containing the building blocks used to run the simulation, as well as the script "em_<model>.R", which runs the estimation algorithm; and a folder "log", in which log files are created and updated during the simulation. In each model's folder, there is a shell script file to run the simulation on a computing cluster and the "main_<model>.R" file, which orchestrates the simulation.



