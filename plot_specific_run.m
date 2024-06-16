%% plot_specific_run

% Terleth, May 2024. 
clearvars
close all
clc


%% select which model run to choose: 

% proportions of water let through to the base.
array_p = 0.1 ; 

% values of [times] k 
array_k =  0.01 ; 


% load the correct mat file 
filename = ['D:/APRIL_envelope/params_gp_' num2str(array_p) '_kp_' num2str(array_k) '.mat'] ; 
        
disp(filename)
load(filename)
         
        %% recalculate additional values 
        C = constants() ;  

         % snip only the second 100 years from everything 
         y = y(C.t_0*t>=100, :) ; 
         t = t(C.t_0*t>=100) ; 

        % recalculate constants 
        K  = C.K * array_k   ;          
        [C,P] = adjust_constants(C,P,K) ; 
   
        % recalculate variables
        Eplus = max(y(:,2)*C.E0,0)/C.E0 ;
        N = min(y(:,1)/C.chi, 1./Eplus) ; 
        u = P.slope^(1/C.p) * y(:,1).^(1+(1/C.p)) .* N.^(-C.q/C.p) ;
        Phi = min(1,(Eplus./(y(:,1)/C.chi))) ; 
        Q = (1/P.l)*(P.slope .* Eplus.^(C.alpha) + (Phi .* P.slope^(1/2) .* y(:,3).^(4/3))) ; 

        % seasonality 
        try 
           Ta_plot = interp1(time, TaS,t,'linear') ; 
        catch 
            Ta_plot = ones(length(t),1) * Ta_A ;
        end 


        % melt C.DDF2 implements the concentration of all the melt to
        % only when air temperature is >0 deg C! 
        try 
            m_plot = max(0, C.DDF2 * (Ta_plot - C.Tm) / C.a_0) ; 
        catch 
            m_plot = max(0, C.DDF * (Ta_plot - C.Toffset) / C.a_0) ;
        end 
        
        % compute how much melt actually reached the bed 
        Beta = min(max(array_p,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 
        m_bed_plot = m_plot .* Beta ; 

        
        %% visualise  

        fig2 = figure  ; 

        % axis limits 
        xlims = [100,110] ;
        ylims = [-.05,1.1] ; 
        
        % melt 
        plot(t*C.t_0,m_plot/max(m_plot),'--'), hold on 
        area(t*C.t_0,m_bed_plot/max(m_bed_plot),FaceColor=[0 0.4470 0.7410],FaceAlpha=.9,EdgeColor='none')
        %plot(t*C.t_0 , Beta,'r')

        % E 
        E = y(:,2) ; 
        area(t*C.t_0,E/max(E),FaceColor=[0.4940 0.1840 0.5560],FaceAlpha=0.2,EdgeColor='none')

        % velocity 
        plot(t*C.t_0,u/max(u),'LineWidth',2,Color=[0.4660 0.6740 0.1880])

        % N 
        plot(t*C.t_0,N/max(N),'LineWidth',2,Color=[0.9290 0.6940 0.1250])

        % S 
        S = y(:,3) ; 
        plot(t*C.t_0,S/max(S),'LineWidth',2,Color=[0.6350 0.0780 0.1840])

        % Q 
        plot(t*C.t_0,Q/max(Q),'LineWidth',2,Color=[0.8500 0.3250 0.0980])

        % H 
        H = y(:,1) ; 
        plot(t*C.t_0, H/max(H),'LineWidth',1.5,Color='k')

        xlim(xlims)
        ylim(ylims)
        ylabel('$(y/y_{max}$)',Interpreter='latex')
        xlabel('time (years)',Interpreter='latex')  
        ax = gca ; 
        box(ax,"on")
        ax.YAxis.Color = 'k' ; 
        ax.XAxis.Color = 'k' ; 
      
        grid on ; 
        lg = legend('w surf','w bed', 'E', 'u', 'N', 'S', 'Q','H', Interpreter='latex', Location='EastOutside'); 
        title([num2str(array_k),'$\times K$ and P=',num2str(array_p)],Color='k',Interpreter='latex')   
    
