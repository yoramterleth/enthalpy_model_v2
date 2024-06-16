%% plot_return_periods

% Terleth, May 2024. 

% warning, the fft method is not perfect for finding the best return
% interval. It is a good way to quickly glance at what the data does but it
% is a good idea to manually check the model output and assess what the
% main periodicity really is. 

close all
clearvars
clc

addpath('C:\Users\Yoram\OneDrive - University of Idaho\Desktop\matlab_helpers\')

%% give arrays to plot.

% values of beta_min
array_p = sort([0:0.025:0.4,0.5,0.75,1,0.01,0.04,0.06,0.09]) ; 

% values of [times] k 
array_k =  [0.001,0.0025,0.003,0.004, 0.005,0.006,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1,5, 10] ; % [0.0025,0.0075,0.02] ;  % [0.001,0.005,0.01,0.05,0.1,0.5,1] ; 

% pre assign grid 
envelope = zeros([length(array_k),length(array_p)]) ; 

%% now run through all the .mat fiels to pull out the return period.

for jj = 1:length(array_p)

    for ii = 1:length(array_k)

        % load the correct mat file 
        load(['D:/APRIL_envelope/params_gp_' num2str(array_p(jj)) '_kp_' num2str(array_k(ii)) '.mat'])

        % snip only the second 100 years from everything 
        y = y(C.t_0*t>=100, :) ; 
        t = t(C.t_0*t>=100) ; 

        % recalculate values
        h = y(:,1) ; 
        s = y(:,3) ; 
        Eplus = max(y(:,2)*C_mod.E0,0)/C_mod.E0 ;
        N = min(y(:,1)/C_mod.chi, 1./Eplus) ; 
        u = P.slope^(1/C_mod.p) * y(:,1).^(1+(1/C_mod.p)) .* N.^(-C_mod.q/C_mod.p) ;

        
        if ~isempty(u)

            % Perform FFT
            N = length(u);  % Length of the signal
            Fs = 1 / (t(2)*C.t_0 - t(1)*C.t_0);  % Sampling frequency in samples/year 
            f = Fs*(0:(N/2))/N;  % Frequency vector
            Y = fft(u);  % Compute FFT
            P = abs(Y/N);  % Compute power spectrum
            P = P(1:round(N/2+1));  % Take only the positive frequencies
            P(2:end-1) = 2*P(2:end-1);  % Double the power (except DC and Nyquist)

            % use find peaks to find the main periodicity
            [~, idx] = findpeaks(P, 'SortStr','descend','NPeaks',3)  ; 

            main_periodicities = 1 ./ f(idx);  % Convert frequencies to periods
            disp('Main Periodicities:');
            disp(main_periodicities);

            % assign the longest of th 3 as the main return period
            return_period = max(main_periodicities(~isinf(main_periodicities))) ; 


            % assign that return period to the grid that will be
            % plotted.
            envelope(ii,jj) = return_period ; 

        else 
            envelope(ii,jj) = nan ; 
        end 

    end 
end 

%% figure     
    
fig = figure ; 


sp = pcolor(array_p,array_k,envelope); shading flat ; hold on  ; 
sp.FaceColor='interp' ; 
sp.FaceAlpha = .8 ; 
contour(array_p,array_k,envelope,[2,2],'k',LineWidth=2)
contour(array_p,array_k,envelope,[10,10],'k--',LineWidth=2)
contour(array_p,array_k,envelope,[100,100],'k:',LineWidth=2)
colormap(plasma(128)) ; 
c = colorbar ; 
ax = gca ; 

ax.YAxis.Scale = 'log' ; 
xlabel('$\beta_{min}$ ', Interpreter='latex') 
ylabel('$K$ ($\times K_{0}$)',Interpreter='latex')
ylabel(c,{'highest amplitude', 'speedup return period (years)'},Interpreter='latex')
ax.XLim = [0,.25] ; 
ax.YLim = [0.005,10] ; 
ax.XAxisLocation = 'top' ; 
ax.ColorScale ='log';
yticks([0,0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1,2,5,10]) ; 
xticks(0:0.03:.3) ; 
grid on 
set(gca, 'Layer','Top')
ax.Color = [rgb('grey'),0.3] ; 


legend('','2 years', '10 years', '100 years')
