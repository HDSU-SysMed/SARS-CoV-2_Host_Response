function fun_plot_ts_host_gene_profile(gname,genes_lfc,m_host_lfc,sem_host_lfc,res)

    fs = 14;
    time_points =    [0     1     2     4     7    12    24    48];
    x = time_points;
    colors = lines(4);
    
    t = log([.5 1 2 4 7 12 24 48])'; % log-scaled time points
                                     % (due to the parameter limits, the
                                     % time of the first point doesn't
                                     % matter)
    tfine = linspace(min(t),max(t),100);
    
    xl = log([.4 55]);
    xlabels =[{'Ctr'};cellstr(num2str([time_points(2:end)']))];

    hf = figure; 
    set(hf,'Position',[100 100 600 400])
    hold on; 

    col_fc = find(strcmp(genes_lfc,gname));
    col_pr = find(strcmp(res(:,1),gname),1,'first');
    if ~isempty(col_fc)
        y = m_host_lfc(col_fc,:);
        s = sem_host_lfc(col_fc,:);
        errorbar(t,y,s,'k','LineWidth',1.5)
        set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs,...
            'XLim',xl,'XTick',t,'XTickLabel',xlabels)
        % plot profile
        beta = res{col_pr,3}; fun_type = res{col_pr,2};
        if fun_type==1 % sigmoidal increase
            b = [beta(1:2) log(beta(3))];
            y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3))));
        elseif fun_type==2 % sigmoidal increase before sigmoidal decrease
            b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
            y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3)))) - ...
                b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));
        elseif fun_type==3 % sigmoidal decrease
            b = [beta(1:2) log(beta(3))];
            y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3))));
        elseif fun_type==4 % sigmoidal decrease before increase        
            b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
            y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3)))) + ...
                b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));        
        end

        plot(tfine,y_mod,'Color',colors(fun_type,:),'LineWidth',1.5)
 
        grid on;
        plot(-.5*ones(1,2),get(gca,'YLim'),'k:')
        xlabel('Time in hours','FontSize',fs)
        ylabel('Log2-fold change','FontSize',fs)
    else
        text(mean(x)-1,.5,sprintf('Data\ninsufficient'))
        set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs,...
            'XLim',xl,'XTick',x,'XTickLabel','')
    end
    title(gname,'FontSize',fs)

end