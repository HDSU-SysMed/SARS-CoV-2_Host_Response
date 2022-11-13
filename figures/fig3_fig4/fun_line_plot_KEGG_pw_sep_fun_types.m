function fun_line_plot_KEGG_pw_sep_fun_types(pw_id,alpha,pw_list,res_table,~,...
    ~,lin_or_log,save_fig)

    pw_name = pw_list{pw_id,1};

    % only select genes with log2-f.c. exceeding [-alpha, alpha]
    ind = (cell2mat(res_table(:,1))==pw_id)&(cell2mat(res_table(:,6))>alpha);
    N_tot = sum(cell2mat(res_table(:,1))==pw_id);
    pw_tmp = res_table(ind,:);
    cm = [.8 0    0; 1    0.5    0;...
        .0196    0.04    0.65; 0.2627    0.5765    0.85];
    fs1 = 22;
    fs2 = 14;
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
    
    hf = figure;
    set(hf,'Position',[50 50 1700 950])

    if ~isempty(pw_tmp)
        Ng = size(pw_tmp,1);
        [~,ord]=sort(cell2mat(pw_tmp(:,4)));
        pw_tmp = pw_tmp(ord,:);
        sel_genes = pw_tmp(:,3);
        fun_types = cell2mat(pw_tmp(:,4));
        im_tmp = zeros(Ng,arlen);
        resp_times = zeros(Ng,2);
        for ii=1:Ng
           tp = fscale(pw_tmp{ii,5},tmin,tmax,arlen);
           tp(2) = min(tp(2),arlen); % only display until tmax
          
           % visualize expression profile by color
           if fun_types(ii,1)==1
               im_tmp(ii,max(tp):end)=1;
               resp_times(ii,1) = max(pw_tmp{ii,5});
           elseif fun_types(ii,1)==2
               im_tmp(ii,tp(1):tp(2))=2;
               resp_times(ii,:) = pw_tmp{ii,5}(1,:);
           elseif fun_types(ii,1)==3
               im_tmp(ii,1:max(tp))=3;
               resp_times(ii,1) = max(pw_tmp{ii,5});
           elseif fun_types(ii,1)==4
               im_tmp(ii,1:min(tp))=4;       
               im_tmp(ii,max(tp):end)=4;       
               resp_times(ii,:) = pw_tmp{ii,5}(1,:);
           end       
        end
        sp_num = zeros(1,4);
        for ii=1:4
            subplot(1,4,ii); hold on
            % re-order each block according to response times     
            ind = max(im_tmp,[],2)==ii;
            im_block = im_tmp(ind,:);
            yl_block = sel_genes(ind,:);
            Nb = size(im_block,1);
            sp_num(1,ii)=Nb;
            [~,ord] = sort(resp_times(ind,1));
            im_ord = im_block(ord,:);
            y_labels = yl_block(ord,:);
        
            imagesc(im_ord,[0 4])
            colormap([.8*ones(1,3);cm])
            if Nb>0
                set(gca,'box','on','TickDir','out','YTick',0:Nb,'YTickLabel','',...
                    'YLim',[.5 Nb+.5],'XTick',xt,'XTickLabel',xlabels,'FontSize',fs1,'XLim',[1 arlen+.5])
                for jj=1:Nb
                    text(1.02*arlen,jj,y_labels{jj,1},'FontSize',fs2,...
                        'HorizontalAlignment','left')
                end
            end
            axis ij
            xlabel('Time in hours','FontSize',fs1)
        end
    else
        title([pw_name, sprintf('\ndata insufficient')],'FontSize',fs1)
    end
    for ii=1:4
        % reduce dimensions according to number of genes
        subplot(1,4,ii)
        pos = get(gca,'Position');
        dy_new = sp_num(1,ii)/max(sp_num);
        set(gca,'Position',pos.*[1 1 1 1.01*dy_new])
    end
    if save_fig==1
        im_name = ['Fig_pw' num2str(pw_id) '_'...
            replace(pw_name,{' ','(',')','-',',','\','/'},'_'),'_sep_fun_types'];
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r100',im_name)
        close(hf)
    end   
end