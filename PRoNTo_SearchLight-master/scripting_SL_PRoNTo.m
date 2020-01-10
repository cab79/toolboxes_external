% Test variation of accuracy with region size.
% Usee radii of increasing size: 5, 10, 15, 20, 25

% RECOMMEND TO USE WITH NO MORE THAN 2-FOLD CROSS-VALIDATION FOR SPEED

%opt.searchtype = 'spacetime';Rlist = 5;%5:5:25;
opt.searchtype = 'time';Rlist = [inf inf 10];%5:5:25;
opt.parallel = 1;
opt.i_model = 1; % Model index to use
opt.loadF = 1; % load all features
opt.savImg = 0;
opt.plot = 0;
opt.permStat = 0; % permutations
Pin = 'C:\Data\Catastrophising study\SPMstats\pronto\t-3000_-2_b-3000_-2500_m_-2200_-1800_Grp_Exp_Subject_orig_cleaned_SPNall_prt_Exp_earlycue2_gpc_ROI\PRT.mat';
for ii = 1:numel(Rlist)
    opt.R = Rlist(ii,:);
    [SLres{ii},Pout{ii}] = crc_parSLdev(Pin,opt); %#ok<*SAGROW>
end

