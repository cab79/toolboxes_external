%% Save regression results into images
function Pout = save_regression_images(Vmsk,pth,PRTw,SLres,DIM,nVx,i_model,R,lVx)
% 1/ prepre files: Vmse, Vcorr, Vr2
Vmse = Vmsk;
Vmse.fname = fullfile(pth, ...
    sprintf('SL_R%s_mse_%s.nii',turn_num2char(R,2),PRTw.model(i_model).model_name));

Vmse.dt(1) = 16; % save as float
Vmse.descrip = 'PRoNTo Search Light MSE';
Vmse = spm_create_vol(Vmse);
Vcorr = Vmse;
Vcorr.fname = fullfile(pth, ...
    sprintf('SL_R%s_corr_%s.nii',turn_num2char(R,2),PRTw.model(i_model).model_name));
Vcorr.descrip = 'PRoNTo Search Light Correlation';
Vcorr = spm_create_vol(Vcorr);
Vr2 = Vmse;
Vr2.fname = fullfile(pth, ...
    sprintf('SL_R%s_r2_%s.nii',turn_num2char(R,2),PRTw.model(i_model).model_name));
Vr2.descrip = 'PRoNTo Search Light R2';
Vr2 = spm_create_vol(Vr2);

% 2/ prepare & collect values
val_mse = zeros(prod(DIM),1)+NaN;
val_corr = zeros(prod(DIM),1)+NaN;
val_r2 = zeros(prod(DIM),1)+NaN;
for ivx=1:nVx
    if ~isempty(SLres(ivx).mse)
        i_vx = lVx(ivx);
        val_mse(i_vx) = SLres(ivx).mse;
        val_corr(i_vx) = SLres(ivx).corr;
        val_r2(i_vx) = SLres(ivx).r2;
    end
end
val_mse = reshape(val_mse,DIM);
val_corr = reshape(val_corr,DIM);
val_r2 = reshape(val_r2,DIM);
% 3/ save values
Vmse = spm_write_vol(Vmse,val_mse);
Vcorr = spm_write_vol(Vcorr,val_corr);
Vr2 = spm_write_vol(Vr2,val_r2);
Pout{2} = Vmse.fname;
Pout{3} = Vcorr.fname;
Pout{4} = Vr2.fname;

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