addpath ./mfile/

ccc
%% raw data path
files = dir('./RawData/*.mat');

%% change the bellowing for running a different dataset
para.dir.raw_file   = fullfile(files.folder, files.name);
para.dir.output_dir = './ReconData';
if ~isfolder(para.dir.output_dir)
    mkdir(para.dir.output_dir)
end

%% set parameters
para.setting.ifplot = 1;        % plot convergence during reconstruction 
para.setting.ifGPU  = 0;        % set to 1 when you want to use GPU
para.Recon.time_frames = 101:200;   %set to 'all' for reconstructe all time frames (could take a while)

%% set parameters
para.weight_tTV = 0.08;                         % temporal TV regularizaiton parameter (normalized by F^T d)
para.weight_sTV = 0;                            % spatial TV regularizaiton parameter (normalized by F^T d)

para.Recon.narm         = 2;                    % number of arms per time frame
para.Recon.FOV          = 1.25;                 % reconstruction FOV
para.Recon.epsilon      = eps('single');        % small vale to avoid singularity in TV constraint
para.Recon.step_size    = 2;                    % initial step size
para.Recon.noi          = 150;                  % number of CG iterations
para.Recon.type         = '2D Spiral server';   % stack of spiral
para.Recon.break        = 1;                    % stop iteration if creteria met. Otherwise will run to noi
para.Recon.matrix_size  = [84, 84];             % acquisition image matrix size

%% do the recon

methods = {'ES', 'LF', 'no'};

for i = 1:3
    method = methods{i};
    para.Recon.method = method;
    para.dir.save_recon = fullfile(para.dir.output_dir, sprintf('%s_narm_%g_t_%.5f_s_%.5f_iter_%g_%s.mat', files.name(1:end-4), para.Recon.narm, para.weight_tTV, para.weight_sTV, para.Recon.noi, method));
    reconstruction(para);
end

