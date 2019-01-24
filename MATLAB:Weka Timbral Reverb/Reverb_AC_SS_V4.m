% This script adds the necessary paths for the Two!Ears Auditory Front-End and RT60 estimation. extraction framework and set up the needed pathes. Make sure to clear
% yourself the Matlab workspace, if that is necessary. This script also generate the ".arff" file and will load the model and test set for evaluation. The output will 
%of this script is the prediction for each test file. 
clear all;
cd('wav files directory'); %location which contains all the wav files
[filenames, pathname, filterindex] = uigetfile( '*.wav', 'WAV-files (*.wav)', 'Pick a file', 'MultiSelect', 'on'); 
cd('output directory');
fid = fopen('outputfile.arff','w'); %name of the outfile file which will contain all the extracted features
fprintf(fid, '@relation RTlandRTRandBLandsREVandFLandpClarandLlowandITDfandpASWandITDbandpLEVandsLEV\n@attribute RTL numeric\n@attribute RTR numeric\n@attribute BL numeric\n@attribute sREV numeric\n@attribute FL numeric\n@attribute pClar numeric\n@attribute Llow numeric\n@attribute ITDf numeric\n@attribute pASW numeric\n@attribute ITDb numeric\n@attribute pLEV numeric\n@attribute sLEV numeric\n@attribute reverb {L, H}\n\n@data\n');
for K = 1 : length(filenames) %loop over the test files
  thisfullname = fullfile(pathname, filenames{K});
  [y,fs]=audioread(thisfullname);  % Reading Audio files
         parConf.mu_psi = -0.034;            % Defining parameters
         parConf.mu_psi_dip = -160.39*1e-3; % for estimation
         parConf.T_min = 63.1*1e-3;
         parConf.a_P_REV = 0.207; 
         parConf.b_P_REV = 21.4;
         parConf.a_P_CLA = 2.76; 
         parConf.b_P_CLA = 2.84;
         parConf.mu_ASW = 2e-2;      
         parConf.nu_ASW = 5.63e2;    
         parConf.a_P_ASW = 2.41;
         parConf.b_P_ASW = 2.64;
         parConf.mu_LEV = 2.76e-2;   
         parConf.nu_LEV = 6.80e2;    
         parConf.a_P_LEV = 1.17; 
         parConf.b_P_LEV = 2.20;
  [par, psi] = afeRAA(y, fs, parConf); % obtaining room acoustic aspects and store it in the par as a structure.
  fprintf('\n RT estimation over time \n\n')
tic
[ rt_est_L, par2 ] = RT_estimation_my( y(:,1), fs ); %calculating the RT60 for
[ rt_est_R, par2 ] = RT_estimation_my( y(:,2), fs ); %left and right channels
toc
rtL_M=mean(rt_est_L); 
rtR_M=mean(rt_est_R);
  fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,?\n', rtL_M, rtR_M, par.BL, par.sREV, par.FL, par.pClar, par.Llow, par.ITDf, par.pASW, par.ITDb, par.pLEV, par.sLEV); %storing features in the .arff file
end
fclose(fid);
cd('C:\Program Files\Weka-3-8');
dosCommand=('java -classpath weka.jar weka.classifiers.functions.Logistic -p 0 -l Logistic_classifier_Weka_D5.6.model -T C:\Users\Fred\Desktop\temp\outputfile.arff');
[status,cmdout] = dos(dosCommand)