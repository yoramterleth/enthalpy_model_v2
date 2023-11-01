%% make surging envelopes of Ta and a. 


%% 
close all
clear all
clc 
addpath('C:\Users\Yoram\OneDrive - University of Idaho\Desktop\matlab_helpers\')

%% PATH TO SAVE TO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_path = 'D:\OCT\' ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load constants 

C = constants() ; 

%% define vectors
a_array = 0:0.1:1 ; 
Ta_array = -15:1:0 ; 

%% define other parameters
l = 10e3 ; 
slope = 1 ; 

seasonal_amp = 15 ; 

P.u1 = 10 ;              % m a-1 
P.u2 = 100 ;        % m a-1


%% normalise other parameters 
P.l = l / C.l_0 ; 
P.slope = slope / 1 ; 

%% ODE solver parameters  

% intial conditions (H/Ho and E/Eo and S/So)
init = [1,1,1] ; 

% timespan (t1/to tend/t0) 
time = 0:1e-3:15 ; % used to be 1e-4 timestep 

% solver options 
options = odeset('RelTol',1e-6,'Stats','off','OutputFcn',[]) ;  %   @odeplot) ; % @odephas2); % @odepphas2 


%% initialise output grid 
dHdend = zeros(length(Ta_array),length(a_array)) ; 
dEdend = zeros(length(Ta_array),length(a_array)) ; 

envelope = nan(length(Ta_array),length(a_array)) ; 

%% iterate over options 

for i = 1:length(a_array)
    for j = 1:length(Ta_array) 


        %% define changing parameters 
        % changing 
        a = a_array(i) ; 
        Ta = Ta_array(j); 

        %% Apply seasonality to Ta  
        [TaS,TaConstant] = make_seasonality(time*C.t_0,Ta,seasonal_amp) ; 
        timeS = time ; 

        %% normalise parameters 
        P.a = a / C.a_0 ; 
        TaSnorm = TaS / C.T_0 ; 
        P.l = l / C.l_0 ; 
        P.slope = slope / 1 ; 
        
        % compute melt 
        m = C.DDF * (Ta * C.T_0 - C.Toffset) / C.a_0 ; 
        P.m = max(m,0) ;
        disp('Normalised parameters:')
        disp(P)
        
        %% compute the adjusted DDF
        T1 = max(TaConstant - C.Toffset,0) ; 
        T2 = max(TaSnorm .* C.T_0 - C.Tm,0) ;  
        DDF2 = C.DDF .* sum(T1)/sum(T2) ; 
        
        C.DDF2 = DDF2 ; 

        %% solve ODE 
        [t,y] = ode23s(@(t,y)ode_CASE3_seasonality(t,y,C,P,timeS,TaSnorm),time,init, options) ; 
        
        %% save the associated workspace
        disp('Saving...')
        save([save_path 'params_a_' num2str(a_array(i)) '_Ta_' num2str(Ta_array(j)) '.mat'])
   
    end 
end 

%save('D:/saved_runs_mega_43/FINAL_workspace.mat')
