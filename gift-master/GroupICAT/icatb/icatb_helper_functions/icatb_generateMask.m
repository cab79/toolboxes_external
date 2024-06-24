function icatb_generateMask(files, varargin)
%% Generate default mask. Default mask is generated by averaging the individual subject masks. Individual subject masks includes voxels above or equalling the mean.
% After the mask is generated, batch file is written out using the current mask.
%
% Inputs:
%
% 1. files - Character array of file names
% 2. varargin - Variable number of arguments
%   a. multiplier - Multiplier applied to the mean. Default is 1.
%   b. threshold - Average mask threshold.
%   c. prefix - Output prefix
%   d. outputDir - Output directory
%   e. corr_threshold - Threshold to exclude outliers based on correlation
%

icatb_defaults;
global FONT_COLOR;
global DEFAULT_MASK_SBM_MULTIPLIER;

modalityType = icatb_get_modality;

if (strcmpi(modalityType, 'eeg'))
    error('Mask feature is not supported for EEG modality');
end

%% Initialze vars
% Default mask multiplier (data >= mult*mean)
mult = 1;
% Mask Threshold
threshold = 0.7;
% Prefix
prefix = 'ica_analysis';
% Output directory
outputDir = '';

if (strcmpi(modalityType, 'smri'))
    mult = DEFAULT_MASK_SBM_MULTIPLIER;
end


for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'multiplier'))
        mult = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'prefix'))
        prefix = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'outputdir'))
        outputDir = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'corr_threshold'))
        corr_threshold = varargin{n + 1};
    end
end

showGUI = 0;
if (isempty(outputDir))
    outputDir = icatb_selectEntry('title', 'Select analysis output directory', 'typeEntity', 'directory');
end

if (isempty(outputDir))
    outputDir = pwd;
end

if (~exist('files', 'var') || isempty(files))
    showGUI = 1;
    files = icatb_selectEntry('title', 'Select Nifti files of all subjects ...', 'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.img;*.nii');
end

drawnow;

files = cellstr(files);

if (showGUI)
    
    prompt= {'Enter mask multiplier (voxels >= multiplier*mean(voxels)', 'Enter average mask threshold', 'Enter output prefix'};
    defaultanswer = {num2str(mult), num2str(threshold), prefix};
    
    answers = icatb_inputdlg2(prompt, 'Mask options', 1, defaultanswer);
    
    if (isempty(answers))
        error('Input gui is closed');
    end
    
    mult = str2num(answers{1});
    threshold = str2num(answers{2});
    prefix = answers{3};
    
end

drawnow;

disp('Generating mask ...');

%time_points = zeros(1, length(files));

masks = cell(1, length(files));

filesN = files;

parfor nF = 1:length(files)
    currentFileN = deblank(files{nF});
    if (size(currentFileN, 1) == 1)
        [pathstr, fN, extn] = fileparts(currentFileN);
        filesP = icatb_listFiles_inDir(pathstr, [fN, extn]);
        if (isempty(filesP))
            error(['Files doesn''t exist. Please check the file pattern ', currentFileN]);
        end
        filesP = icatb_fullFile('directory', pathstr, 'files', filesP);
    else
        filesP = currentFileN;
    end
    filesP = icatb_rename_4d_file(filesP);
    %time_points(nF) = size(filesP, 1);
    filesP = deblank(filesP(1,:));
    disp(['Loading ', filesP, ' ...']);
    [tmp_dat, HInfo] = icatb_loadData(filesP);
    filesN{nF} = filesP;
    tmp_mask = double(tmp_dat >= mult*mean(tmp_dat(:)));
    
    masks{nF} = tmp_mask(:);
    
    %     if (nF == 1)
    %         masks = zeros(numel(tmp_mask), length(files));
    %     end
    
    %masks(:, nF) = tmp_mask(:);
    
end


[Vw, HInfo] = icatb_returnHInfo(deblank(filesN{1}(1, :)));


masks = cat(2, masks{:});

avgMaskT = mean(masks, 2);

correlations = abs(icatb_corr(avgMaskT, masks));

avgMaskT = double(avgMaskT >= threshold);

if (isempty(find(avgMaskT == 1)))
    error('No voxels found in brain. Change the threshold and mask multiplier settings');
end

avgMaskT = reshape(avgMaskT, HInfo.DIM(1:3));

[L, NUM] = icatb_spm_bwlabel(avgMaskT, 26);

comps = zeros(1, NUM);

for ii = 1:NUM
    comps(ii) = length(find(L == ii));
end

[numInds, max_inds] = max(comps);

avgMask = double(L == max_inds);


%% Write mask
writeV = HInfo.V(1);
writeV.n(1) = 1;
maskFile = fullfile(outputDir, [prefix, 'Mask.nii']);
%writeV.fname = maskFile;
%icatb_write_vol(writeV, avgMask);
icatb_write_nifti_data(maskFile, writeV, avgMask);

axesTitle = 'Individual mask correlations w.r.t average mask';
fH = icatb_getGraphics(axesTitle, 'timecourse', 'mask_corr', 'on');
axH = axes('units', 'normalized', 'position', [0.15, 0.15, 0.7, 0.7]);
set(fH, 'resize', 'on');
plot((1:length(files)), correlations, 'm', 'linewidth', 1.5, 'parent', axH);
title(axesTitle, 'parent', axH);
set(axH, 'Xcolor', FONT_COLOR);
set(axH, 'Ycolor', FONT_COLOR);
xlabel('Subjects', 'parent', axH);
ylabel('Correlations', 'parent', axH);
axis(axH, 'tight');


if (~exist('corr_threshold', 'var'))
    prompt={'Enter correlation threshold to remove subject outliers from the analysis.'};
    name = 'Correlation threshold';
    numlines = 1;
    defaultanswer = {num2str(0.8)};
    answer = icatb_inputdlg2(prompt, name, numlines, defaultanswer);
    
    if (~isempty(answer))
        corr_threshold = str2num(answer{1});
    end
end

filesIn = find(correlations >= corr_threshold);


masks = reshape((masks == 1), [HInfo.DIM(1:3), size(masks, 2)]);
save(fullfile(outputDir, [prefix, '_mask_info.mat']), 'masks', 'avgMask', 'mult', 'threshold', 'corr_threshold', 'correlations');

if (~isempty(filesIn))
    
    files = files(filesIn);
    
    
    [eig_summary, num_eigs] = get_eig_summary(files, maskFile);
    
    save(fullfile(outputDir, [prefix, '_mask_info.mat']), 'eig_summary', '-append');
    
    chk_eigs = find(num_eigs >= 2);
    
    if (isempty(chk_eigs))
        disp('Eigen values summary: ');
        disp(eig_summary);
        error('Files don''t have enough degrees of freedom to do group ICA.');
    end
    
    
    files =  files (chk_eigs);
    
    min_df = min(num_eigs);
    
    numComp = min([min_df, 20]);
    if (~strcmpi(modalityType, 'smri'))
        PC1 = min([min_df, ceil(1.5*numComp)]);
    else
        PC1 = numComp;
    end
    
    
    %     if (~strcmpi(modalityType, 'smri'))
    %         time_points = time_points(filesIn);
    %     else
    %         time_points = length(files);
    %     end
    %     numComp = min([min(time_points), 20]);
    %     if (~strcmpi(modalityType, 'smri'))
    %         PC1 = min([min(time_points), ceil(1.5*numComp)]);
    %     else
    %         PC1 = numComp;
    %     end
    
    %% Write input data to a text file
    txtFileName = fullfile(outputDir, [prefix, '_input_files.txt']);
    dlmwrite(txtFileName, char(files), '');
    
    %% Write batch file with necessary info
    batchInfo.modalityType = modalityType;
    comments{1} = '%% Modality Type';
    batchInfo.outputDir = outputDir;
    comments{end + 1} = '%% Output Directory';
    batchInfo.prefix = prefix;
    comments{end + 1} = '%% All the output files will be preprended with the specified prefix.';
    batchInfo.perfType = 1;
    comments{end + 1} = '%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m.';
    batchInfo.txtFileName = txtFileName;
    comments{end + 1} = '%% Input data file names written to a text file. ';
    batchInfo.dataSelectionMethod = 4;
    comments{end + 1} = '%% Data selection option. If option 4 is specified, file names must be entered in input_data_file_patternss';
    batchInfo.input_data_file_patterns = 'textread(txtFileName, ''%s'', ''delimiter'', ''\n'');';
    comments{end + 1} = '%% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of data-sets.';
    if (~strcmpi(modalityType, 'smri'))
        batchInfo.input_design_matrices = {};
        comments{end + 1} = '%% Design matrix/matrices.';
        batchInfo.dummy_scans = 0;
        comments{end + 1} = '%% Number of dummy scans.';
    end
    batchInfo.maskFile = maskFile;
    comments{end + 1} = '%% Full file path of the mask file.';
    if (~strcmpi(modalityType, 'smri'))
        batchInfo.backReconType = 'gica';
        comments{end + 1} = '%% Back-reconstruction type.';
        batchInfo.preproc_type = 1;
        comments{end + 1} = '%% Data pre-processing option. By default, Remove Mean Per Timepoint is used';
        batchInfo.numReductionSteps = 2;
        comments{end + 1} = '%% Number of data reduction steps used.';
    end
    batchInfo.scaleType = 2;
    comments{end + 1} = '%% Scale components. By default, components are converted to z-scores.';
    batchInfo.algoType = 'infomax';
    comments{end + 1} = '%% ICA algorithm to use. Infomax is the default choice.';
    batchInfo.numOfPC1 = PC1;
    comments{end + 1} = '%% Number of principal components used in the first PCA step.';
    if (~strcmpi(modalityType, 'smri'))
        batchInfo.numOfPC2 = numComp;
        comments{end + 1} = '%% Number of principal components used in the second PCA step. Also the number of independent components extracted from the data.';
    else
        comments{end} = [comments{end}, ' Also the number of independent components extracted from the data.'];
    end
    
    if (~strcmpi(modalityType, 'smri'))
        batchFileName = fullfile(outputDir, [prefix, '_gift_batch_file.m']);
    else
        batchFileName = fullfile(outputDir, [prefix, '_sbm_batch_file.m']);
    end
    
    fnames = fieldnames(batchInfo);
    
    %defs = repmat({''}, length(fnames), 1);
    defs = {};
    for nF = 1:length(fnames)
        tmp = batchInfo.(fnames{nF});
        if (isempty(tmp))
            tmp = [fnames{nF}, ' = {};'];
        else
            if (~strcmpi(fnames{nF}, 'input_data_file_patterns'))
                if (isnumeric(tmp))
                    tmp = sprintf('%s = %s;', fnames{nF}, mat2str(tmp, 4));
                else
                    tmp = sprintf('%s = ''%s'';', fnames{nF}, tmp);
                end
            else
                tmp = ['input_data_file_patterns = ', tmp];
            end
        end
        
        defs{end + 1} = comments{nF};
        defs{end + 1} = tmp;
        defs{end + 1} = '';
        
    end
    
    dlmwrite(batchFileName, char(defs), '');
    disp(char(['Input information is saved in batch file ', batchFileName, '.'], 'Edit the batch file according to your needs and ', ...
        ['use icatb_batch_file_run(''', batchFileName, '''); to run the batch analysis']));
    fprintf('\n\n');
    
else
    
    warning('!!!None of the subjects surpassed this threshold');
    
end


function [eig_summary, num_eigs] = get_eig_summary(files, maskFile)
%% Get eigen values summary
%

num_eigs = zeros(length(files), 1);
eig_summary = cell(length(num_eigs), 1);

tol = eps('double');

parfor nF = 1:length(files)
    
    disp(['Reading subject ', num2str(nF), ' ...']);
    drawnow;
    current_files = files{nF};
    if (size(current_files, 1) == 1)
        [pathstr, fp, extn] = fileparts(current_files);
        current_files = icatb_listFiles_inDir(pathstr, [fp, extn]);
        current_files = icatb_fullFile('directory', pathstr, 'files', current_files);
    end
    data = icatb_remove_mean(icatb_read_data(current_files, [], maskFile));
    D = eig(icatb_cov(data));
    
    D = abs(D);
    
    %     if (nF == 1)
    %         tol = eps(class(D));
    %     end
    
    nums = length(find(D > tol));
    
    
    num_eigs(nF) = nums;
    
    [p, fn, extn] = fileparts(deblank(current_files(1,:)));
    extn = icatb_parseExtn(extn);
    
    eig_summary{nF} = [fullfile(p, [fn, extn]), ', ', num2str(nums)];
    
end

