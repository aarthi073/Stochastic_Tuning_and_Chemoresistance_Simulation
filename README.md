## Stochasting Tuning Simulation in Response to Cytotoxicity
A MATLAB simulation using the Gillespie algorithm to model tumor resistance and random gene expression through fitness-driven epigenetic changes that optimize overall cell health through cell plasticity.

# Project Breakdown
Initiaize gene, RNA, and protein states with an array of zeros and a tracked single RNA molecule and 
protein with an empty array. 

Continuous Monte Carlo simulation to model coupled reactions: drug-induced gene activation and inactivation,
transcription and RNA degradation, translation and degradation, cell division, and drug-induced death. Reactions are tracked at the cell adn population level, and the sampled waiting time for the next reaction is distributed exponentially. 
The number of RNA, protein, and cells depends on the weight of the reaction for each computational run. 

1. Defines various dosages: 1 micrometer, 10 micrometers, 50 micrometers, 250 micrometers, 500 micrometers, and 1000 micrometers.
2. Simultaneously tracks single-cell trajectories and population-level behavior to analyze stochasticity at teh gene level while analyzing population patterns to determine tumor resistance. Gene activation and cell death were the only reactions that were affected by the drug concentration.
 



Drug-protein interactions follow a sigmoidal function illustrated by a hill curve that models the inhibition rate
where the resistance protein acts as a competitive inhibitor to drug-induced apoptosis. 

Parameters and equations are explain in "Simulation_Parameters".

# Repository Structure
|   README.md
|   Simulation Parameters.pdf
|
+---docs
|       Random Walks Project Paper.pdf
|
+---src
|       Stochastic_Tuning_in_Response_to_Cytotoxicity.m
|
+---Stochastic_Tuning_and_Chemoresistance_Simulation
|       README.md
|
\---Tables
        Simulation Parameters.csv


#Outputs and Interpretation
To optimize a general analysis of stochasticity in RNA production, protein production, and final populatio, the 
inactivation rate was set to 100, protein level to overcome the effects of the drug (K) were set to 0.001, 
transcription rate was set to 100, translation rate was set to 70, division rate was set to 0.02, mRNA
decay rate was set to 0.1, protein decay rate of 0.5, and drug_delta is set to 0.1 (Refer to Simulation_Parameters table). 

Overall, the RNA production pattern fluctuates and transends to the protein expression pattern, generation
substantial heterogeneity as time persists. The "bursting" pattern in transcription and translation are random
with each run for the given parameters, but this suggests that the overall trend generates noise in variable
RNA and protein numbers. Tumor population survival response against different drug conentration varies for different hill coefficients. Therefore, population dynamics are sensitive to transcription regulation and stochasticity in response to varying drug concentrations.

There is susbtantial overdispersion relative to the Poisson distribution in all trajectories of 
RNA and protein distributions at the population level while the fano factor (ratio of variance to the mean)
is much greater than 1. 

The results are consistent with the hypothesis that stochastic tuning contributes to phenotypic heterogeneity and transient resistance as a stress response to cytotoxicity. 
Non-genetic variability can be a targetted process in oncology treatments instead of focusing on any one stage of the disease, as resistance prevails against evolving drugs.


#Limitations
Consideration of cancer stages, mesenchymal-epithelial shift dynamics, limited number of RNA, and limited number of proteins. Additionally,
treatment cycles, drug administration, and rest periods were ignored in this simulation.

#Author
Aarthi Bharathan B.S. Computational & Applied Mathematics and Statistics in Mathematical Biology, College of William & Mary

