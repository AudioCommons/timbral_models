cd('C:\Users\Fred\Documents\AudioComm\Andy\Deliverable 5.3 archive\Deliverable 5.3 archive\Listening test data\WithinSource\Reverb\Stimuli');
[filenames, pathname, filterindex] = uigetfile( '*.wav', 'WAV-files (*.wav)', 'Pick a file', 'MultiSelect', 'on');
%filenames = cellstr(filename);   %in case only one selected
for K = 1 : length(filenames)
  thisfullname = fullfile(pathname, filenames{K});
  [speechIn6,FS6]=audioread(thisfullname);
  [rt,drr,cte,cfs,edt] = irStats(thisfullname);
  cd('C:\Users\Fred\Documents\MATLAB\AudioComm');
  csvwrite([filenames{K}(1:end-4),'.csv'],[rt;drr;cte;edt]);
  %savefig([filenames{K}(1:end-4),'.fig']);
  %speechIn6 = myVAD(speechIn6);
  %fMatrix6(1,o) = {mfccf(ncoeff,speechIn6,FS6)};
end