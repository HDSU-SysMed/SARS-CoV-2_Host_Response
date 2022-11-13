function res = fun_line_plot_KEGG_pw(pw_id,alpha,pw_list,res_table,lin_or_log,show_plots,save_fig)

    pw_name = pw_list{pw_id,1};

    % only select genes with log2-f.c. exceeding [-alpha, alpha]
    ind = (cell2mat(res_table(:,1))==pw_id)&(cell2mat(res_table(:,6))>alpha);
    N_tot = sum(cell2mat(res_table(:,1))==pw_id);
    pw_tmp = res_table(ind,:);
    if ~isempty(pw_tmp)
        N_pw = pw_tmp{1,end};
    end
    cm = [.8 0    0; 1    0.5    0;...
        .0196    0.04    0.65; 0.2627    0.5765    0.85];
    fs = 22;
    n_per_h = 10;
    arlen = 48*n_per_h;

    % scale for displaying as line plot
    if lin_or_log==1
        tmin = 0; tmax = 48;
        fscale = @(t,tmin,tmax,arlen) 1+round((arlen-1) * n_t/tmax);
        xt = [0 12 24 36 48]*n_per_h;
        xlabels = {'0','12','24','36','48'};
    elseif lin_or_log==2
        tmin = .5; tmax = 48;
        fscale = @(t,tmin,tmax,arlen) 1 + max(round((arlen-1) * (log(t)-log(tmin))/(log(tmax)-log(tmin))),0);
        xt = fscale([.5 1 2 4 7 12 24 48],tmin,tmax,arlen);
        xlabels = {'0.5' '1' '2' '4' '7' '12' '24' '48'};
    end
    
    if show_plots==1
        hf = figure;
        set(hf,'Position',[50 50 1200 800])
    end
    if ~isempty(pw_tmp)
        Ng = size(pw_tmp,1);
        [~,ord]=sort(cell2mat(pw_tmp(:,4)));
        pw_tmp = pw_tmp(ord,:);
        sel_genes = pw_tmp(:,3);
        fun_types = cell2mat(pw_tmp(:,4));
        im_tmp = zeros(Ng,arlen);
        resp_times = zeros(Ng,1);
        for ii=1:Ng
           tp = fscale(pw_tmp{ii,5},tmin,tmax,arlen);
           tp(2) = min(tp(2),arlen); % only display until tmax

           if fun_types(ii,1)==1
               im_tmp(ii,max(tp):end)=1;
               resp_times(ii,1) = max(pw_tmp{ii,5});
           elseif fun_types(ii,1)==2
               im_tmp(ii,tp(1):tp(2))=2;
               resp_times(ii,1) = pw_tmp{ii,5}(1,1);
           elseif fun_types(ii,1)==3
               im_tmp(ii,1:max(tp))=3;
               resp_times(ii,1) = max(pw_tmp{ii,5});
           elseif fun_types(ii,1)==4
               im_tmp(ii,1:min(tp))=4;       
               im_tmp(ii,max(tp):end)=4;       
               resp_times(ii,1) = pw_tmp{ii,5}(1,1);
           end       
        end

        im_ord = zeros(Ng,arlen);
        y_labels = cell(Ng,1);

        % re-order each block according to response times
        for ii=1:4
            ind = max(im_tmp,[],2)==ii;
            im_block = im_tmp(ind,:);
            yl_block = sel_genes(ind,:);
            [~,ord] = sort(resp_times(ind,1));
            im_ord(ind,:) = im_block(ord,:);
            y_labels(ind,:) = yl_block(ord,:);
        end

        % create lines for pathway to indicate (de)activation as color
        % intensity, store numbers of regulated genes, profile types and
        % response times
        line_pos_tmp = sum((im_tmp==1)|(im_tmp==2),1);
        line_neg_tmp = sum(im_tmp>=3,1);
        num_genes_reg = [sum(fun_types<=2) sum(fun_types>=3) N_tot];
        res = {pw_name, line_pos_tmp,line_neg_tmp,num_genes_reg,[fun_types resp_times],y_labels};
        
        if show_plots==1
            imagesc(im_ord,[0 4])
            colormap([.8*ones(1,3);cm])
            set(gca,'box','on','TickDir','out','XTick',xt,...
                'XTickLabel',xlabels,'FontSize',fs,'XLim',[.5 arlen+.5],...
                'Position',get(gca,'Position').*[1 1 .9 .9]+[.05 .05 0 0])
            for ii=1:Ng
                text(1.01*arlen,ii,y_labels{ii,1},'FontSize',round(fs * (.6 + .4/Ng)),...
                    'HorizontalAlignment','left')
            end

            xlabel('Time in hours','FontSize',fs)
            ylabel('Genes in pathway','FontSize',fs)
            perc_affected = 100*Ng/N_tot;
            text(arlen/2,.5-.1*Ng,pw_name,'HorizontalAlignment','center','FontSize',fs)
            text(arlen/2,.5-.05*Ng,[num2str(Ng) ' of ' num2str(N_tot)...
                ' genes affected (' num2str(round(perc_affected)) '%), ' ...
                num2str(N_pw) ' genes in pathway'],'HorizontalAlignment','center','FontSize',fs-4)
        end
    else
        if show_plots==1
            title([pw_name, sprintf('\ndata insufficient')],'FontSize',fs)
        end
        res = cell(1,6);
    end

    if (show_plots==1)&&(save_fig==1)
        im_name = ['Fig_pw' num2str(pw_id) '_'...
            replace(pw_name,{' ','(',')','-',',','\','/'},'_')];
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r100',im_name)
        close(hf)
    end   
end