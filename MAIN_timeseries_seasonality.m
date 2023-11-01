%% mAIN script to model timeseries with seasonality introduced.  

% this script can be used to evaluate all three cases

%% 
close all
clear all
clc 

%% load constants 

C = constants() ; 

% set case nb 
case_nb = 3 ; 

%% define parameters 
a = 0.7 ; 
Ta = -8; 
l = 10e3 ; 
slope = 1 ; 
seasonal_amp = 15 ; 

P.u1 = 0 ;              % m a-1 
P.u2 =  100 ; %100 ;        % m a-1
    
%% ODE solver parameters  

% timespan (t1/to tend/t0) 
time = 0:1e-4:5 ; 

% intial conditions (H/Ho and E/Eo)
if case_nb == 3 
    init = [1,1,1] ; 
else
    init = [1.3,0] ;
end 

% solver options 
options = odeset('RelTol',1e-6,'Stats','on','OutputFcn',@odeplot) ; % @odephas2); % @odepphas2 \

%% Apply seasonality to Ta

% convert time vector to acutal times in years  
 
[TaS,TaConstant] = make_seasonality(time*C.t_0,Ta,seasonal_amp) ; 
timeS = time ; 

% plot the seasonality 
figure
plot(timeS*C.t_0, TaS)
xlim([0,5])

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
%C.a_0 = 10 ; 

%% solve ODE 

if case_nb == 1 
    [t,y] = ode45(@(t,y)ode_CASE1_seasonality(t,y,C,P,timeS,TaSnorm),time,init, options) ; 

elseif case_nb == 2 
    [t,y] = ode45(@(t,y)ode_CASE2_seasonality(t,y,C,P,timeS,TaSnorm),time,init, options) ; 

elseif case_nb == 3
    [t,y] = ode23s(@(t,y)ode_CASE3_seasonality(t,y,C,P,timeS,TaSnorm),time,init, options) ;
end 


% figure 
% scatter(y(:,1),y(:,2),40,t)
% xlim([0,3]),ylim([-1,2])

%% visualise 
if case_nb == 3 
    [fig] = plot_timeseries_CASE3(C,P,t,y)  ;
else 
    [fig] = plot_timeseries(C,P,t,y) ; 
end 


[fig2] = plot_seasonality(C,P,TaSnorm,timeS) ; 


%% save workspace 
 % save([pwd '/saved_runs/case3_seasonality_05_07_2023.mat'])
