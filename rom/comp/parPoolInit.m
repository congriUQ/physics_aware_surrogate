function [ppool] = parPoolInit(N_Threads)
%% Initializes parallel pool
current_pool = gcp('nocreate');
if(~numel(current_pool))
    %Create with N_Threads workers
    delete(gcp('nocreate'));
    
    %Create cluster object
    pc = parcluster('local');
    %Create temporary folder for par job data
    timestr = datestr(now, 'mmddHHMMSSFFF');
    if(strcmp((java.net.InetAddress.getLocalHost.getHostName),...
            'workstation1-room0436') || ...
            strcmp((java.net.InetAddress.getLocalHost.getHostName),...
            'constantin-ThinkPad-T430s') || ...
            strcmp((java.net.InetAddress.getLocalHost.getHostName),...
            'master'))
        foldername = strcat('/home/constantin/.par_job_files/', timestr);
    else
        %on subnodes via ethernet
        foldername = strcat('/home_eth/constantin/.par_job_files/', timestr);
    end
%     foldername = strcat('/home/constantin/.par_job_files/', timestr);
    mkdir(foldername);
    pc.JobStorageLocation = foldername;
    
    if(nargin == 0 || N_Threads > pc.NumWorkers)
        N_Threads = pc.NumWorkers;
    end
    ppool = parpool(pc, N_Threads);
else
    ppool = gcp;
end
ppool.IdleTimeout = Inf;

end

