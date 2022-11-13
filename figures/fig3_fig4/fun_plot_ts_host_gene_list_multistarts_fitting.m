function res = fun_plot_ts_host_gene_list_multistarts_fitting(list_name,glist,genes_host,m_host,sem_host,show_plot,save_fig)

    fs = 12;
    np_per_fig = 24;
    t = log([.5 1 2 4 7 12 24 48])'; % log-scaled time points
                                     % (due to the parameter limits, the
                                     % time of the first point doen't
                                     % matter)
    tfine = linspace(min(t),max(t),100);
    colors = [.8 0    0; 1    0.5    0;...
        .0196    0.04    0.65; 0.2627    0.5765    0.85];   
    xlabels =[{'Ctr'}; cellstr(num2str(exp(t(2:end))))];
    
    sel_genes = unique(glist);
    Ng = numel(sel_genes);
    res = {};
    disp(['Evaluate ' list_name])
    if Ng<np_per_fig
        hf = figure;
        set(hf,'Position',[50 50 1800 940])
        np = ceil(sqrt(numel(sel_genes)));
        for ii=1:Ng
            if show_plot==1
                subplot(np,np,ii); hold on;
            end
            col = find(strcmp(genes_host,sel_genes{ii}));
            if ~isempty(col)
                y = m_host(col,:)';
                s = sem_host(col,:)';
                ind = ~isinf(y);
                if sum(ind)>=6
                    res_fitting = fun_do_multistarts(t(ind,1),y(ind,1),s(ind,1));
                    ind_optimal = find(cell2mat(res_fitting(:,2))==min(cell2mat(res_fitting(:,2))),1,'first');
                    % save fun_type, chisq, parameters
                    res = cat(1,res,[sel_genes(ii,1) {ind_optimal} res_fitting(ind_optimal,:)]);

                    disp_message = fun_report_results(ind_optimal,res_fitting{ind_optimal,1});
                    disp([sel_genes{ii,1} ', fun ' num2str(ind_optimal) ', ' disp_message ...
                        ', ' num2str(ii) '/' num2str(Ng)])

                    if show_plot==1
                        errorbar(t,y,s,'ko','LineWidth',2)
                        % plot alternative fits, indicate optimal fit
                        for fun_type=1:4
                            beta = res_fitting{fun_type,1};
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
                            if fun_type==ind_optimal
                                plot(tfine,y_mod,'Color',colors(fun_type,:),'LineWidth',1.5)
                            else
                                plot(tfine,y_mod,'--','Color',colors(fun_type,:),'LineWidth',1)    
                            end
                        end    
                        yl = [min(y-s)-.5 max(y+s)+.5];
                        set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs,...
                            'XLim',[min(t) max(t)]+[-.1 .1],'YLim',yl,'XTick',t,'XTickLabel',xlabels)
                        grid on;
                        xlabel('Time in hours','FontSize',fs)
                        ylabel(sprintf('log2\nfold change'),'FontSize',fs)
                        title(sel_genes{ii},'FontSize',fs)
                    end
                else
                    if show_plot==1
                        text(0,.5,'Data insufficient','FontSize',fs-3)
                        set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs-3,...
                             'XLim',[min(t) max(t)]+[-.1 .1],'XTick',t,'XTickLabel','','YTickLabel','')
                        title(sel_genes{ii},'FontSize',fs)
                    end
                end
            else
                if show_plot==1
                    text(0,.5,'Data insufficient','FontSize',fs-3)
                    set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs-3,...
                         'XLim',[min(t) max(t)]+[-.1 .1],'XTick',t,'XTickLabel','','YTickLabel','')
                    title(sel_genes{ii},'FontSize',fs)
                end
            end
        end
        % use last tile for legend
        subplot(np,np,Ng+1); hold on;
        plot(0,0,'ko','LineWidth',2)
        for cc=1:4
            
            plot(1,1,'Color',colors(cc,:),'LineWidth',2)
        end
        set(gca,'XLim',[-1 1],'YLim',[-1 1])
        legend({'Exp. data','f_1, sig. increase','f_2, sig. incr./decr.',...
            'f_3, sig. decrease','f_4, sig. decr./incr.'},...
            'Location','North','FontSize',fs)
        set(gca,'Visible','off')

        if (show_plot==1)&&(save_fig==1)
            im_name = ['Fig_host_genes_' list_name];
            set(gcf,'PaperPositionMode','auto')
            print('-dpng','-r200',im_name)
            set(gcf,'PaperSize',[45 30]);
            print(im_name,'-dpdf')
        end
    else 
        nfig = ceil(Ng/np_per_fig);
        for jj=1:nfig
            if show_plot==1
                np = ceil(sqrt(np_per_fig));
                hf = figure;
                set(hf,'Position',[50 50 1800 900])
            end
            ind_p = setdiff(((jj-1)*np_per_fig + 1):(jj*np_per_fig),Ng:(nfig*np_per_fig));
            for ii=ind_p
                if show_plot==1
                    subplot(np,np,ii-(jj-1)*np_per_fig); hold on;
                end
                col = find(strcmp(genes_host,sel_genes{ii}));
                if ~isempty(col)                    
                    y = m_host(col,:)';
                    s = sem_host(col,:)';
                    ind = ~isinf(y);
                    if sum(ind)>=6
                        res_fitting = fun_do_multistarts(t(ind,1),y(ind,1),s(ind,1));
                        ind_optimal = find(cell2mat(res_fitting(:,3))==min(cell2mat(res_fitting(:,3))),1,'first');
                        % save fun_type, chisq, parameters
                        res = cat(1,res,[sel_genes(ii,1) {ind_optimal} res_fitting(ind_optimal,:)]);

                        disp_message = fun_report_results(ind_optimal,res_fitting{ind_optimal,1});
                        disp([sel_genes{ii,1} ', fun ' num2str(ind_optimal) ', ' disp_message ...
                            ', ' num2str(ii) '/' num2str(Ng)])

                        if show_plot==1
                            errorbar(t,y,s,'ko','LineWidth',2)
                            % plot alternative fits, indicate optimal fit
                            for fun_type=1:4
                                beta = res_fitting{fun_type,1};
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
                                if fun_type==ind_optimal
                                    plot(tfine,y_mod,'Color',colors(fun_type,:),'LineWidth',2)
                                else
                                    plot(tfine,y_mod,'--','Color',colors(fun_type,:),'LineWidth',1)    
                                end
                            end

                            % set lower and upper limit to better distinguish
                            % genes with weak and strong effects
                            yl = [min([-1 min(y-s)-.1]) max([1 max(y+s)+.1])];
                            set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs,...
                                'XLim',[min(t) max(t)]+[-.1 .1],'YLim',yl,'XTick',t,'XTickLabel',xlabels)
                            grid on;
                            xlabel('Time in hours','FontSize',fs)
                            ylabel(sprintf('log2\nfold change'),'FontSize',fs)  
                            title(sel_genes{ii},'FontSize',fs)
                        end
                    else % too small number of experimental data points
                        if show_plot==1
                            text(0,.5,'Data insufficient','FontSize',fs-3)
                            set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs-3,...
                                'XLim',[min(t) max(t)]+[-.1 .1],'XTick',t,'XTickLabel','','YTickLabel','')
                            title(sel_genes{ii},'FontSize',fs)
                        end
                    end
                else % no basal value defined
                    if show_plot==1
                        text(0,.5,'Data insufficient','FontSize',fs-3)
                        set(gca,'box','on','LineWidth',1,'TickDir','out','FontSize',fs-3,...
                            'XLim',[min(t) max(t)]+[-.1 .1],'XTick',t,'XTickLabel','','YTickLabel','')
                        title(sel_genes{ii},'FontSize',fs)
                    end
                end  
            end
            if show_plot==1
                % use last tile for legend
                subplot(np,np,np^2); hold on;
                plot(0,0,'ko','LineWidth',2)
                for cc=1:4
                    plot(1,1,'Color',colors(cc,:),'LineWidth',2)
                end
                set(gca,'XLim',[-1 1],'YLim',[-1 1])
                legend({'Exp. data','f_1, sig. increase','f_2, sig. incr./decr.',...
                    'f_3, sig. decrease','f_4, sig. decr./incr.'},...
                    'Location','North','FontSize',fs)
                set(gca,'Visible','off')
            end
        
            if (show_plot==1)&&(save_fig==1)
                im_name = ['Fig_host_genes_' list_name '_' num2str(jj)];
                set(gcf,'PaperPositionMode','auto')
                print('-dpng','-r600',im_name)
                set(gcf,'PaperSize',[45 30]);
                print(im_name,'-dpdf')
            end 
        end
    end
end

function res_fitting = fun_do_multistarts(t,y,s)
    seed = 1;
    mod_var = @(b,X,Y,fun_type) fun_obj_RNAseq_profiles(b,X,Y,fun_type);
    N_runs = 100;
    options = optimoptions('lsqnonlin','Display','none',...
        'MaxIter',1e2,'MaxFunEvals',1e2,'TolFun',1e-6, 'TolX',1e-6);
    rng(seed,'twister')

    res_fitting = cell(4,3);
    fun_AICcorr = @(chisq,n,k) chisq + n*log(2*pi) + 2*k + (2*k^2 + 2*k)/(n-k-1);
    fun_BIC = @(chisq,n,k) chisq + n*log(2*pi) + k*log(n);
    for fun_type=1:4
        if (fun_type==1)||(fun_type==3) % sigmoidal
            x0 = [1 20 7]; lo = [.1 1 1]; hi = [10 100 48];
        elseif (fun_type==2)||(fun_type==4) % two sigmoidals
            x0 = [1 20 4 1 4]; lo = [.1 1 1 .1 4]; hi = [10 100 48 10 48]; 
        end
        problem = createOptimProblem('lsqnonlin','x0',log(x0),'lb',log(lo),'ub',log(hi),...
            'objective',@(beta)mod_var(beta,t,[y s],fun_type),'options',options);
        ms = MultiStart('UseParallel',true,'Display','none');
        [log_b,fval,~,~,~] = run(ms,problem,N_runs);
        res_fitting(fun_type,:)={exp(log_b),fun_AICcorr(fval,numel(y),numel(log_b)),fun_BIC(fval,numel(y),numel(log_b))};
    end    
end

function mess = fun_report_results(ind,res_par)
    mess = [];
    if ind==1
        mess = ['Induction at ' sprintf('%.2g',res_par(3)) 'h'];
    elseif ind==2
        mess = ['Induction at ' sprintf('%.2g',res_par(3)) 'h' ...
            ', inhibition at ' sprintf('%.2g',sum(res_par([3 5]))) 'h'];
    elseif ind==3
        mess = ['Inhibition at ' sprintf('%.2g',res_par(3)) 'h'];
    elseif ind==4
        mess = ['Inhibition at ' sprintf('%.2g',res_par(3)) 'h' ...
            ', induction at ' sprintf('%.2g',sum(res_par([3 5]))) 'h'];
    end
end






