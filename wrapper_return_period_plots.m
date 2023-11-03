%% wrapper script to plot return interval comparisons 

% yoram terleth 
% 21/08/2023

%% plot speedup frequency of the modelled output runs in parameter space 
% this script is to visualise the resulting surge envelopes of 4 model configuration options:
%     - no seasonality, no runoff to the bed at zero velocity
%     - no seasonality, 10% runoff to the bed at zero velocity
%     - seasonality, no runoff to the bed at zero vel
%     - seasonality, 10% runoff to the bed at zero vel

% it collects the four configurations from respective folders, that should be generated each in individual with MAIN_envelope_seasonality.m


% coding of filenames: the 10 or zero is for how much melt is let through
% at velocity 0 ! So 10 means beta2, melt is let through
%% 
close all
clear all
clc 
addpath([pwd '\matlab_helpers\'])

%% input paths to model run outputs here: %%%%%%%%%%%%%%%
pathlist = {'D:/SEP/no_seasonality_B1_0/', 'D:/SEP/no_seasonality_B1_10/',...
    'D:/OCT/seasonality_B1_0/','D:/OCT/seasonality_B1_10/'} 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% titles for subplots:
title_list = {'no seasonality, $\beta_{1}$', 'no seasonality, $\beta_{2}$', ...
    'seasonality, $\beta_{1}$', 'seasonality, $\beta_{2}$'} ; 

fig = figure ; 

for plt = 1:4 

% path to draw from %
save_path_cur = string(pathlist(plt)); %%


%% define vectors ! this should be similar to those defined in MAIN_envelope_seasonality.m !
a_AA = 0:0.1:1 ; 
Ta_AA = -15:1:0 ;

% make matrix to fill in with return periods
return_intervals = ones(length(Ta_AA),length(a_AA)) ; 

% initialise legend
leg=[] ; 

% loop over grids
for ii = 1:length(a_AA)
    for jj = 1:length(Ta_AA)

        
        %% select file 
        filename_cur = strcat(save_path_cur, 'params_a_', num2str(a_AA(ii)), '_Ta_', num2str(Ta_AA(jj)), '.mat') ;

        disp(filename_cur)
        load(filename_cur) ; 
        

        % build legend entry 
        leg = [leg ; {['a= ' num2str(a_AA(ii)) ', Ta= ' num2str(Ta_AA(jj))]}] ; 

        %% test if output is real 
        if ~isreal(y)
            disp('Skipping because non real output.')
            return_intervals(jj,ii) = nan ; 
            continue
        end 

         %% recalculate additional values 
         Eplus = max(y(:,2)*C.E0,0)/C.E0 ;
         N = min(y(:,1)/C.chi, 1./Eplus) ; 
        
        u = P.slope^(1/C.p) * y(:,1).^(1+(1/C.p)) .* N.^(-C.q/C.p) ;

        % find velocity peaks: they need to be at least 10% of the maximal
        % velocities.
        [pks, idx] = findpeaks(u); 
        idx(pks/max(u)<0.1)= nan ; 
        pks(pks/max(u)<0.1)= nan ; 

        % remove peaks that were deemed too small
        idx = rmmissing(idx) ; 

        % quantify the median return interval between remaining velocity peaks
        t_idx = t(idx).*C.t_0 ; 
        re_t = median(diff(t_idx)) ; 
        
        % optional add on, comment out if you don't mind shrinking glaciers...
        % gets rid of the result if the long term ice thickness is < .6 of its initial value at the end of the simulation 
        % i.e., the glacier is in a non sustainable climatic regime (upper left corner of the envelope)
        if y(end,1)/max(y(:,1)) < 0.6 
            return_intervals(jj,ii) = nan ; 
        elseif isnan(y(end,1)/max(y(:,1))) 
            return_intervals(jj,ii) = nan ;    
        else 
            return_intervals(jj,ii) =  re_t ; % (max(t)*C.t_0)/length(rmmissing(pks)) ; % re_t ; % 
        end

    end 
end 

%% viualisation 
subplot(2,2,plt) ; 

pcolor(a_AA,Ta_AA/C.T_0,return_intervals), shading flat , hold on 
cmp = flipud(magma(128)) ;% cbrewer('div','RdBu',128,'linear') ; 
colormap(cmp)
xlabel('$a/a_{0}$',Interpreter='latex')
ylabel('$Ta/Ta_{0}$',Interpreter='latex')
% cb = colorbar ; 
% ylabel(cb, 'velocity peak return period (years)',Interpreter='latex')
set(gca, 'colorscale','log')
clim([1,10^3])

contour(a_AA,Ta_AA/C.T_0,return_intervals,[8 8],Color="#77AC30",Linewidth=2,ShowText='on',LabelFormat='%1.0f years')
contour(a_AA,Ta_AA/C.T_0,return_intervals,[70 70],Color="#77AC30",Linewidth=2,ShowText='on',LabelFormat='%1.0f years')

set(gca,'Color',rgb('light grey'))
grid on 
set(gca, 'Layer','top')
title(string(title_list(plt)), Interpreter='latex') 

if plt ==4 
    cb = colorbar ; 
    cb.Position = [.92,.32,0.01,.4] ; 
    ylabel(cb, 'velocity peak return period (years)',Interpreter='latex')


end 
fontsize(fig, 14,"points")


