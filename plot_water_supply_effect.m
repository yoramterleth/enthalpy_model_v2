%% make figure about effect of water input

close all 
clearvars
clc


%% path to draw from %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_paths = {'D:/SEP/no_seasonality_B1_0/','D:/SEP/no_seasonality_B1_10/',...
    'D:/OCT/seasonality_B1_0/','D:/OCT/seasonality_B1_10/'} ; %% 

beta_index = [0, 10, 0, 10]; 

title_list = { 'no seasonality, $\beta_{1}$', 'no seasonality, $\beta_{2}$', ...
    'seasonality, $\beta_{1}$', 'seasonality, $\beta_{2}$'} ; 

xlim_zoom_array = [44,44.5; 44.7,45.2 ; 49,49.5 ; 44,44.5] ; 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define vectors
%a_array = 0:0.1:1 ; 
a_A =  0.6 ; % 0.1 ; 
Ta_A = -8  ; %  -10 ; %-14:1-5 ; 

% window length for pxx smoothing 
wl = 1e3 ; 

% make inset map or not
inset = 0 ; 

leg = [] ; 

fig = figure ; 
tiles = tiledlayout(length(save_paths),4) ; 

for p_nb = 1:length(save_paths)

    save_path = string(save_paths(p_nb)) ; 

for ii = 1:length(a_A)
    for jj = 1:length(Ta_A)

        
        %% select file 
        % filename = [pwd '/seasonality_exploration_august/params_a_' num2str(a_A(ii)) '_Ta_' num2str(Ta_A(jj)) '.mat'] ; 
        %filename = ['D:/saved_runs_mega_4/params_a_' num2str(a_A(ii)) '_Ta_' num2str(Ta_A(jj)) '.mat'] ; 

        filename = strcat(save_path, 'params_a_', num2str(a_A(ii)), '_Ta_', num2str(Ta_A(jj)), '.mat') ;

        disp(filename)
        load(filename) ; 

        % build legend entry 
        leg = [leg ; {['a= ' num2str(a_A(ii)) ', Ta= ' num2str(Ta_A(jj))]}] ; 

        %% recalculate additional values 
    
        Eplus = max(y(:,2)*C.E0,0)/C.E0 ;
        N = min(y(:,1)/C.chi, 1./Eplus) ; 
        
        u = P.slope^(1/C.p) * y(:,1).^(1+(1/C.p)) .* N.^(-C.q/C.p) ;
        
        Phi = min(1,(Eplus./(y(:,1)/C.chi))) ; 
    
    
        Q = (1/P.l)*(P.slope .* Eplus.^(C.alpha) + (Phi .* P.slope^(1/2) .* y(:,3).^(4/3))) ; 

        % this time with seasonality 
        try 
           Ta_plot = interp1(time, TaS,t,'linear') ; 
        catch 
            Ta_plot = ones(length(t),1) * Ta_A ; 
        end 


        % melt %% !! C.DDF2 implements the concentration of all the melt to
        % only when air temperature is >0 deg C! 
        try 
            m_plot = max(0, C.DDF2 * (Ta_plot - C.Tm) / C.a_0) ; 
        catch 
            m_plot = max(0, C.DDF * (Ta_plot - C.Toffset) / C.a_0) ;
        end 
        
        % compute how much melt actually reached the bed 
        if beta_index(p_nb) == 10 
            Beta = min(max(0.1,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 
        else 
            Beta = min(max(0,((u*C.u_0 - P.u1)/(P.u2-P.u1))),1) ; 
        end 
        m_bed_plot = m_plot .* Beta ; 

        %% visualise 
        
        %subplot(length(save_paths),1,p_nb)
        nexttile([1,3])
        xlims =[35,52] ;
        ylims = [-.05,1.1] ; 
        
        % melt 
        plot(t*C.t_0,m_plot/max(m_plot),'--'), hold on 
        area(t*C.t_0,m_bed_plot/max(m_bed_plot),FaceColor=[0 0.4470 0.7410],FaceAlpha=0.3,EdgeColor='none')
        plot(t*C.t_0 , Beta,'r')

        % E 
        E = y(:,2) ; 
        area(t*C.t_0,E/max(E),FaceColor=[0.4940 0.1840 0.5560],FaceAlpha=0.3,EdgeColor='none')

        % velocity 
        plot(t*C.t_0,u/max(u),'LineWidth',2,Color=[0.4660 0.6740 0.1880])
        % find velocity peaks: they need to be at least 30% of the maximal
        % velocities.
%         [pks, idx] = findpeaks(u); 
%         pks(pks/max(u)<0.3)= nan ; 
%         idx(pks/max(u)<0.3)= nan ; 
%         scatter(t(idx)*C.t_0,pks/max(u))

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
        grid on ; 

        if p_nb == length(save_paths)
          xlabel('time (years)',Interpreter='latex')  
%         else
%             h_cur = gca;
%             h_cur.XAxis.Visible = 'off';
        end 

        

        if p_nb == 1
        legend('surf m','bed m','$\beta$', 'E', 'u', 'N', 'S', 'Q','H', Interpreter='latex',Location='SouthWest')
        end 

        %set(gca,'YScale','log')
        %legend('melt', 'enthalpy', 'ice velocity', 'effective pressure', 'channel size', 'discharge')

        

        title(title_list(p_nb),Interpreter='latex',FontWeight='bold')      

        nexttile([1,1])

        xlims_zoom = xlim_zoom_array(p_nb,:) ;
        
        
        % melt 
        plot(t*C.t_0,m_plot/max(m_plot),'--'), hold on 
        area(t*C.t_0,m_bed_plot/max(m_bed_plot),FaceColor=[0 0.4470 0.7410],FaceAlpha=0.3,EdgeColor='none')
        plot(t*C.t_0 , Beta,'r')

        % E 
        E = y(:,2) ; 
        area(t*C.t_0,E/max(E),FaceColor=[0.4940 0.1840 0.5560],FaceAlpha=0.3,EdgeColor='none')

        % velocity 
        plot(t*C.t_0,u/max(u),'LineWidth',2,Color=[0.4660 0.6740 0.1880])
        % find velocity peaks: they need to be at least 30% of the maximal
        % velocities.
%         [pks, idx] = findpeaks(u); 
%         pks(pks/max(u)<0.3)= nan ; 
%         idx(pks/max(u)<0.3)= nan ; 
%         scatter(t(idx)*C.t_0,pks/max(u))

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

        xlim(xlims_zoom)
        ylim(ylims)
       
        h_cur = gca;
        h_cur.YAxis.TickValues = [] ; 
        

        if p_nb == length(save_paths)
          xlabel('time (years)',Interpreter='latex')  
        end 

        if p_nb == 1
        legend('surf m','bed m','$\beta$', 'E', 'u', 'N', 'S', 'Q','H', Interpreter='latex',Location='SouthWest')
        end 

    end 
end 


end 

fontsize(fig, 16,"points")
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';
