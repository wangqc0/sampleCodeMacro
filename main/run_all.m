% This is the highest-level script of the replicable code. Run this script
% to generate all output in the sample project.

% initialization
cd(fileparts(mfilename('fullpath')))

% clean data
job_clean_data = batch('data');
wait(job_clean_data)
delete(job_clean_data);

% analysis
job_analysis = batch('analysis');
wait(job_analysis)
delete(job_analysis);
