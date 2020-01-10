
%% Save classification results into images
function Pout = save_classif_images(Vmsk,pth,PRTw,SLres,DIM,nVx,i_model,R,lVx)
% 1/ prepre files: Vacc, Vbacc, Vcacc(1/2/...)
Vacc = Vmsk;
Vacc.fname = fullfile(pth, ...
    sprintf('SL_R%s_acc_%s.nii',turn_num2char(R,2),PRTw.model(i_model).model_name));
Vacc.dt(1) = 16; % save as float
Vacc.descrip = 'PRoNTo Search Light accuracy';
Vacc = spm_create_vol(Vacc);
Vbacc = Vacc;
Vbacc.fname = fullfile(pth, ...
    sprintf('SL_R%s_bacc_%s.nii',turn_num2char(R,2),PRTw.model(i_model).model_name));
Vbacc.descrip = 'PRoNTo Search Light balanced accuracy';
Vbacc = spm_create_vol(Vbacc);
nClasses = numel(SLres(end).c_acc);
Vcacc(nClasses) = Vacc;
for ii=1:nClasses
    Vcacc(ii) = Vacc;
    Vcacc(ii).fname = fullfile(pth, ...
        sprintf('SL_R%s_cacc%d_%s.nii',turn_num2char(R,2),ii,PRTw.model(i_model).model_name));
    Vcacc(ii).descrip = ['PRoNTo Search Light class ',num2str(ii),' accuracy'];
    Vcacc(ii) = spm_create_vol(Vcacc(ii));
end
% 2/ prepare & collect values
val_acc = zeros(prod(DIM),1)+NaN;
val_bacc = zeros(prod(DIM),1)+NaN;
val_cacc = zeros(prod(DIM),nClasses)+NaN;
for ivx=1:nVx
    if ~isempty(SLres(ivx).acc)
        i_vx = lVx(ivx);
        val_acc(i_vx) = SLres(ivx).acc;
        val_bacc(i_vx) = SLres(ivx).b_acc;
        val_cacc(i_vx,:) = SLres(ivx).c_acc;
    end
end
val_acc = reshape(val_acc,DIM);
val_bacc = reshape(val_bacc,DIM);
val_cacc = reshape(val_cacc,[DIM nClasses]);
% 3/ save values into images
Vacc = spm_write_vol(Vacc,val_acc);
Vbacc = spm_write_vol(Vbacc,val_bacc);
for ii=1:nClasses
    Vcacc(ii) = spm_write_vol(Vcacc(ii),squeeze(val_cacc(:,:,:,ii)));
end
% 4/ return filenames
Pout{2} = Vacc.fname;
Pout{3} = Vbacc.fname;
for ii=1:nClasses
    Pout{3+ii} = Vcacc(ii).fname; %#ok<*AGROW>
end

end


%% turn index into chars for file numbering
function i_char = turn_num2char(i_num,nMxZ)
% Function turning an index into a char, with zero padding
if nargin<2
    % Upper limit on number to convert, by default:
    % 1e10 - 1 = 999999999, i.e. 9 digits
    nMxZ = 9;
    Zpad = '000000000';
else
    nMxZ = round(nMxZ);
    Zpad = repmat('0',1,nMxZ);
end

i_num = round(i_num);
i_char = num2str(i_num);
i_char = [Zpad(1:(nMxZ-length(i_char))) , i_char];

end