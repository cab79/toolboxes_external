% Denoises ERPs and plots the original and denoised average ERPs, coefficients,
% bands, single trials, and contour plot 
function EEG = EP_den_EEG(EEG,sc,elec,den_type,remove,keep,plotchan)
%clc
close all

handles.par.sr = EEG.srate;                              % sampling rate
handles.par.stim = -EEG.xmin*EEG.srate+1;                         % stim
handles.par.samples = EEG.pnts;               % number of samples
handles.par.scales = sc;                           % number of scales
handles.par.plot_type='coeff';           % 'coeff', 'bands', 'single', 'contour'
handles.par.den_type= den_type;     %' do_den' or 'load_den_coeff'  or 'create_den_coeff'
handles.par.auto_den_type='NZT';  % 'Neigh' or 'NZT'

samples=handles.par.samples;
stim=handles.par.stim;
sr=handles.par.sr;
sc=handles.par.scales;
plot_type= handles.par.plot_type;
den_type=handles.par.den_type;
auto_den_type=handles.par.auto_den_type;
max_trials=100;
max_contour=100;

if ~isempty(remove)
    removemat = ones(sc+1,samples);
    rem_samples = dsearchn(EEG.times',1000*remove.times')';
    removemat(remove.scales,rem_samples(1):rem_samples(2))=0;
end
if ~isempty(keep)
    keepmat = ones(sc+1,samples);
    keep_samples = dsearchn(EEG.times',1000*keep.times')';
    keepmat(keep.scales,keep_samples(1):keep_samples(2))=0;
end

den_coeffs_all = [];
coeffs_all = [];
mean_elec=0;
apply_elec = elec;
if isempty(elec)
    elec = 1:EEG.nbchan;
elseif elec==0
    elec=1;
    mean_elec=1;
    apply_elec = 1:EEG.nbchan;
end
for e = elec
    e
    %% Denoising
    %path='C:\CORE\CORE\CORE\Supporting_functions\EP_den_auto\ERPs\';
    %filename='S4_O2T.asc';
    %save_path='C:\CORE\CORE\CORE\Supporting_functions\EP_den_auto\ERPs\den_results\'
    %if ~exist(save_path,'dir');mkdir(save_path);end
    %save_filename='S4_O2T_den';
    if mean_elec==1
        x=reshape(squeeze(mean(EEG.data(:,:,:),1)),EEG.trials*EEG.pnts,1);
    else
        x=reshape(squeeze(mean(EEG.data(e,:,:),1)),EEG.trials*EEG.pnts,1);
    end
    
    x=x(:);
    sweeps =length(x)/samples;
    xx=reshape(x,samples,sweeps)';
    av=mean(xx,1);
    switch den_type
        case 'do_den'
            switch auto_den_type
                case 'Neigh'
                  [coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);
                case 'NZT'
                  [coeff,denav,den_coeff,y]= Run_NZT(av,handles);
            end
        case 'do_den_all' 
            switch auto_den_type
                case 'Neigh'
                  [coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);
                case 'NZT'
                  [coeff,denav,den_coeff,y]= Run_NZT(av,handles);
            end
        case 'load_den_coeff'        
             [coeff,denav,den_coeff,y,yo]= Run_NZT(av,handles);
             [filename, pathname] = uigetfile('*.mat','Select file');
             matfile=load([pathname filename]);
             den_coeff=matfile.den_coeff;
             [denav,y,den_coeff]=st_den(av,den_coeff,handles);
        case 'create_den_coeff'        
             [coeff,denav,~,yo]= Run_NZT(av,handles);
             %[filename, pathname] = uigetfile('*.mat','Select file');
             den_coeff=coeff;
             if ~isempty(remove)
            %    removemat=load(remove);
                den_coeff=den_coeff.*removemat;
             end
             if ~isempty(keep)
             %   keepmat=load(keep);
                den_coeff=den_coeff.*keepmat;
             end
             [denav,y,den_coeff]=st_den(av,den_coeff,handles);
    end
    if strcmp(den_type,'do_den_all') ||  mean_elec==1
        den_coeffs_all(e,:,:) = den_coeff;
        coeffs_all(e,:,:) = coeff;
    else
        YDEN=st_den(x,den_coeff,handles);
        EEG.data(e,:,:) = YDEN';
    end
    if e == plotchan || e == 1
        plot_av = av;
        plot_denav = denav;
        plot_coeff = coeff;
        plot_den_coeff = den_coeff;
    end
    
        
end
if strcmp(den_type,'do_den_all') ||  mean_elec==1
    den_coeff = squeeze(mean(den_coeffs_all,1));
    coeff = squeeze(mean(coeffs_all,1));
    %den_coeff = reshape(max(reshape(coeffs_all,size(coeffs_all,1),size(coeffs_all,2)*size(coeffs_all,3))),size(coeffs_all,2),size(coeffs_all,3));
    for e = apply_elec
        e
        x=reshape(squeeze(mean(EEG.data(e,:,:),1)),EEG.trials*EEG.pnts,1);
        x=x(:);
        YDEN=st_den(x,den_coeff,handles);
        EEG.data(e,:,:) = YDEN';
        if e == plotchan
            plot_denav = mean(YDEN,1);
        end
    end
    plot_coeff = coeff;
    plot_den_coeff = den_coeff;
end
      
%% Plotting
 %plot average
 set(0,'DefaultFigureColor','w')
 figure('Position',[500 500 700 700])
 subplot(7,1,1:2)       
 plot(((1:samples)-stim+1)/sr,plot_av, 'color',[0.6 0.6 0.6])
 hold on
 plot(((1:samples)-stim+1)/sr,plot_denav, 'color','r')
 xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
 title('Average ERP','fontsize',14)
 xlabel('Time (sec)','fontsize',10)

 switch plot_type
    case 'coeff'
       subplot(7,1,3:7)  
       step = 1/(sc+2):1/(sc+2):1;
        for i=1:sc+1
            scaling_factor = 1.5 * max(abs(plot_coeff(i,:))) * (sc+1);
            aux1= plot_coeff(i,:)/ scaling_factor;
            aux2=plot_den_coeff(i,:)/scaling_factor;
            plot(((1:samples)-stim+1)/sr,aux1+step(sc+2-i),'color', [0.6 0.6 0.6])
            hold on
            plot(((1:samples)-stim+1)/sr,aux2+step(sc+2-i),'r')
        end
        for i=1:sc
            texto =['D' num2str(i)];
            text(-1.1,step(sc+2-i)+0.01,texto);
        end
        texto =['A' num2str(sc)];
        text(-1.1,step(1)+0.01,texto);
        axis off
         text(-0.1,0.05,'Time (sec)');
          text(-1.1,0.94,'(b) Wavelet Coefficients');
        
case 'bands'
    subplot(7,1,3:7)  
    step = 1/(sc+2):1/(sc+2):1;
    scaling_factor = 1.5 * max(max(abs(yo))) * (sc+1);
    aux = y/ scaling_factor;
    aux_all = yo/ scaling_factor;
    for i=1:sc+1
        plot(((1:samples)-stim+1)/sr,aux_all(i,:)+step(sc+2-i),'color', [0.6 0.6 0.6])
        hold on
        plot(((1:samples)-stim+1)/sr,aux(i,:)+step(sc+2-i),'r')
    end   
    for i=1:sc
        texto =['D' num2str(i)];
        text(-1.1,step(sc+2-i)+0.01,texto);
    end
    texto =['A' num2str(sc)];
    text(-1.1,step(1)+0.01,texto);
    axis off
        
case 'single'
    subplot(7,1,3:7)  
    nr_sweeps = min(sweeps,max_trials);
    scaling_factor = 1.7 * max(max(abs(xx))) * nr_sweeps;
    step = 1/(nr_sweeps+1):1/(nr_sweeps+1):1;
    aux_all = xx / scaling_factor;
    aux = YDEN / scaling_factor;                  
    for i=1:nr_sweeps;
        plot(((1:samples)-stim+1)/sr,aux_all(1+nr_sweeps-i,:)+step(i),'color', [0.6 0.6 0.6])
        hold on
        plot(((1:samples)-stim+1)/sr,aux(1+nr_sweeps-i,:)+step(i) ,'r')
    end
    for i=1:nr_sweeps
        texto =['#' num2str(1+nr_sweeps-i)];
        text(-1.1,step(i)+0.01,texto);
    end
    axis off
        
case 'contour'
    subplot(7,1,3:7)
    axis on            
    nr_sweeps = min(sweeps,max_contour);
    [c,h]=contourf(((1:samples)-stim+1)/sr,1:nr_sweeps,YDEN(1:nr_sweeps,:),10);      
    set(h,'Edgecolor','none')
 end      


%% Saving
%save([save_path,save_filename],'YDEN','coeff','den_coeff','y','yo','denav','av','xx')


