# LandscapeVirusExpansion
"The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" 
Donnelly and Gilligan 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Instructions for navigating files uploaded to GitHub repository ‘LandscapeVirusExpansion’
The code in the repository generates the results for Fig. 1 Fig. 2 Fig. 3 Fig. 4 Table 1 Fig. S1.1 Fig. S2.1 from the above manuscript as described below.
All relevant code and data are in the repository folders ‘/MatlabFiles’ and ‘/rFiles’
Where the principle .r, .stan, and .m files are cited below we follow the file name with parentheses in bold summarising the main purpose of the code. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fig. 1
‘/MatlabFiles’
‘theorySchematic_Fig1.m’ (graphs)
Notes… above generates (purely illustratively by combining logistic curves) figure with 4 subplots corresponding to Fig. 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fig. 2, Fig. S1.1, Fig. S2.1
‘/MatlabFiles’
‘theorySimulation_Figs2_S1_S2.m’ (simulates and graphs)
(calls ‘aggreg2_wCuttings.m’)
(calls ‘midScapeEvents.m’)
Notes… above generates 3x2 figure… this figure corresponds to Fig. 2, Fig. S1.1, Fig. S2.1 in turn by setting the appropriate scenario:
i.e., 
propCleanSeed=1 for Fig. 2
propCleanSeed=0 and is_discrim=1 for Fig. S2.1 (cutFromWithin=1 for case 1 cutFromWithin=0  for case 2)
propCleanSeed=0 and is_discrim=0 for Fig. S2.2 (cutFromWithin=1 for case 3 cutFromWithin=0  for case 4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fig. 3, Fig. 4, Table 1
‘/rFiles’
‘BayesLandscapeEstimation_Table1.r’ (fits)
(‘waveProfileMEmodel.stan’)
(requires dataset ‘rFiles/DigitiseFiles_Experiment/ExperimentOrientS.csv’ and ‘rFiles/DigitiseFiles_Survey/SurveyCentralOrientS.csv’ and
and ‘rFiles/DigitiseFiles_Survey/SurveyEasternOrientS.csv’)
Notes… above generates Table 1 and in addition ‘MatlabFiles /empiricalLexp_Fig3.m’ (graphs) and ‘MatlabFiles /empiricalLsurv_Fig4.m’ (graphs) produce Fig. 3, Fig. 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

See also Appendix S5 for description of digitising referring to contents of folders:
‘/rFiles/DigitiseFiles_Experiment’
and
‘/rFiles/DigitiseFiles_Survey’

