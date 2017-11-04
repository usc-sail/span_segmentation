% You need to run this script once on the hpc cluster
% Before running, create in your cluster home folder a MATLAB subfolder 
%   and replace the full path here:
matlab_folder = '/auto/rcf-proj2/sn/toutios/MATLAB';
% And replace your usc e-mail address here
email = 'toutios@usc.edu';

MDCSprofile = parallel.cluster.Torque;
MDCSprofile.JobStorageLocation = matlab_folder;
MDCSprofile.ResourceTemplate = sprintf('-l procs=^N^ -l software=MDCS+^N^ -m abe -M %s', email);
MDCSprofile.RshCommand = 'ssh';
MDCSprofile.RcpCommand = 'scp';
MDCSprofile.SubmitArguments='-l walltime=23:59:59 -A lc_sn';

MDCSprofile.saveAsProfile('24h_MultiNodeProfile');

