%% evaluation of count matrices, plot results
% ska, Nov 2022

%% pre-processing data evaluation steps - when set to 'false', pre-processed datasets are loaded
task_read_in_count_matrix = true;  % read in data
task_perform_normalization = false;  % subtract log-scaled linearly interpolated mock values
task_estimate_parameters = false; % estimate parameters for all pathways
task_extract_inflection_points = false; % extract inflection points for line plots
task_create_line_plots = false; % create lines plots for all KEGG pathways
task_print_images = false;  % save figures to png and pdf

% settings for Figs. 4C and 4D
N_best = 25;        % number of pathways displayed
Nmin_pw_map = 10;   % min. number of regulated genes in pathway
% displayed font size
fs = 22;

%% host data
% order of time points
time_points =    [0     1     2     4     7    12    24    48];
x = [-1 log(time_points(2:end))];
xl = [-1.5 max(x)+.5];

if task_read_in_count_matrix == true
    % read in count matrix
    % adapt to directory containing count matrices obtained from
    % pre-processing using DeSeq2
    path_tmp = pwd; % enter path here (e.g., 'C:\ [...] \res_Exp');
    files_tmp = dir(path_tmp);
    res_l2fc = [];
    err_l2fc = [];
    res_l2fc_m = [];
    err_l2fc_m = [];
    time_points_ns = [];
    time_points_ns_m = [];
    all_genes = [];
    for jj=1:size(files_tmp,1)
        if ~isempty(strfind(files_tmp(jj,1).name,'timepoints_X'))
            fname = files_tmp(jj,1).name;
            time_points_ns = cat(1,time_points_ns,str2double(extractBetween(fname,'timepoints_X','h_vs_mock')));
            D = readcell(fullfile(path_tmp,fname));
            if isempty(all_genes)
                all_genes = D(2:end,1);
            end
            res_l2fc = cat(2,res_l2fc,cell2mat(D(2:end,3)));
            err_l2fc = cat(2,err_l2fc,cell2mat(D(2:end,4)));
        end
        if ~isempty(strfind(files_tmp(jj,1).name,'timepoints_m'))
            fname = files_tmp(jj,1).name;
            time_points_ns_m = cat(1,time_points_ns_m,str2double(extractBetween(fname,'timepoints_m','h_vs_mock')));
            D = readcell(fullfile(path_tmp,fname));
            res_l2fc_m = cat(2,res_l2fc_m,cell2mat(D(2:end,3)));
            err_l2fc_m = cat(2,err_l2fc_m,cell2mat(D(2:end,4)));
        end
    end
    [time_points,ord] = sort(time_points_ns);
    res_l2fc = res_l2fc(:,ord); err_l2fc = err_l2fc(:,ord);
    [time_points_m,ord] = sort(time_points_ns_m);
    res_l2fc_m = res_l2fc_m(:,ord); err_l2fc_m = err_l2fc_m(:,ord);
    save('Data_host_transcriptome_log2FC','all_genes','res_l2fc','err_l2fc','time_points',...
        'res_l2fc_m','err_l2fc_m','time_points_m');
else
    load('Data_host_transcriptome_log2FC.mat');
end

%% perform background corraction by subtracting mock values (obtained for 1, 2, 7, 48h by linear interpolation)
if task_perform_normalization==true
    m_bg = interp1([0 4 12 48]', [res_l2fc(:,1) res_l2fc_m]',[0 1 2 4 7 12 24 48]')';
    s_bg_old = interp1([0 4 12 48]', [err_l2fc(:,1) err_l2fc_m]',[0 1 2 4 7 12 24 48]')';
    xi = [0 1 2 4 7 12 24 48]';

    % fit error model to estimate error for all time points
    seed = 1; % always use the same rng seed
    y_bg = [res_l2fc(:,1) res_l2fc_m]';
    err_y = [err_l2fc(:,1) err_l2fc_m]';
    err_fun = @(b,X,Y) fun_obj_linear_err(b,X,Y);
    N_runs = 100;
    options = optimoptions('lsqnonlin','Display','none','MaxIter',1e2,'MaxFunEvals',1e2,'TolFun',1e-6, 'TolX',1e-6);
    rng(seed,'twister')
    x0_err = [.05 .05]; lo_err = [1e-5 1e-5]; hi_err = [1 1];
    N_g = size(res_l2fc,1);
    s_bg = zeros(N_g,size(res_l2fc,2));
    
    % for interpolated data as displayed in Fig. 2:
    y_fine = interp1(xi,res_l2fc',(0:48)')';
    m_bg_fine = interp1([0 4 12 48]', [res_l2fc(:,1) res_l2fc_m]',(0:48)')';
    s_bg_fine_d = zeros(N_g,49);
    s_bg_fine_m = zeros(N_g,49);

    fvals_m = zeros(N_g,1);
    fvals_d = zeros(N_g,1);
    cvals = zeros(N_g,2);
    dvals = zeros(N_g,2);

    for ii=1:N_g
        ym = m_bg(ii,:); ym_fine = m_bg_fine(ii,:);
        y = y_bg(:,ii); e = err_y(:,ii); 
        yd = res_l2fc(ii,:)'; ed = err_l2fc(ii,:)'; yd_fine = y_fine(ii,:)';

        % fit error model to mock data that is used to estimate
        % 1. error of interpolated mock values at eight measured time points
        % 2. error of interpolated mock values at 1-hour steps between 0 and 48h
        % coefficients in 'cvals', residuals in 'fvals_m'

        problem = createOptimProblem('lsqnonlin','x0',log(x0_err),'lb',log(lo_err),'ub',log(hi_err),...
            'objective',@(beta)err_fun(beta,y,e),'options',options);      
        ms = MultiStart('UseParallel',true,'Display','none');
        [log_c,fval,~,~,~] = run(ms,problem,N_runs);
        c = exp(log_c);
        err_i = c(1)*abs(ym) + c(2)*max(abs(ym));
        err_i_fine = c(1)*abs(ym_fine) + c(2)*max(abs(ym_fine));
        fvals_m(ii,1)=fval;
        cvals(ii,:)=c;
        s_bg(ii,:)=err_i;
        s_bg_fine_m(ii,:)=err_i_fine;

        % fit error model to data of infected cells that is used to estimate
        % error of interpolated measurements in infected cells at 1-hour steps between 0 and 48h
        % coefficients in 'dvals', residuals in 'fvals_d'
        problem = createOptimProblem('lsqnonlin','x0',log(x0_err),'lb',log(lo_err),'ub',log(hi_err),...
            'objective',@(beta)err_fun(beta,yd,ed),'options',options);
        ms = MultiStart('UseParallel',true,'Display','none'); % set to 'false' if parallel computing toolbox n. a.
        [log_d,fval,~,~,~] = run(ms,problem,N_runs);
        d = exp(log_d);
        err_i_fine = d(1)*abs(yd_fine) + d(2)*max(abs(yd_fine));
        s_bg_fine_d(ii,:)=err_i_fine;
        fvals_d(ii,1)=fval;
        dvals(ii,:)=d;
        if mod(ii,100)==0
            disp([num2str(ii) '/' num2str(N_g)])
        end
    end

    save('Data_mean_bg_interpolation','m_bg','s_bg','y_fine','m_bg_fine','s_bg_fine_m','s_bg_fine_d',...
        'fvals_m','cvals','fvals_d','dvals')
else
    load('Data_mean_bg_interpolation.mat')
end

res_l2fc_corr = res_l2fc - m_bg;
err_l2fc_corr = sqrt(err_l2fc.^2 + s_bg.^2);
res_l2fc_corr_fine = y_fine - m_bg_fine;
err_l2fc_corr_fine = sqrt(s_bg_fine_d.^2 + s_bg_fine_m.^2);

%% estimate parameters for all pathways
if task_estimate_parameters == true
    % for saving time, create list of unique entries and then associate entries
    % with pathways
    res = fun_plot_ts_host_gene_list_multistarts_fitting('all genes',all_genes,all_genes,...
            res_l2fc_corr,err_l2fc_corr,0,0);
    save('Results_profiles_all_genes_interpolation','all_genes','res')
else
    load('Results_profiles_all_genes_interpolation.mat')
end

% load genes and pathways from KEGG database (December 2020)
load('Database_KEGG_pw_genes.mat')

%% plot exemplary figures
% examples for figure:
% fun_type 1, ATF3
% fun_type 2, SHOX2 
% fun_type 3, NPAS2
% fun_type 4, RFX6
test_genes = {'ATF3','SHOX2','NPAS2','RFX6'};
for ii=1:numel(test_genes)
    fun_plot_ts_host_gene_profile(test_genes{ii},all_genes,res_l2fc_corr,err_l2fc_corr,res)
end

%% extract inflection points for line plots
if task_extract_inflection_points == true
    res_table = {};
    t = log([.5 1 2 4 7 12 24 48])'; % log-scaled time points
                                     % (due to the parameter limits, the
                                     % time of the first point doen't
                                     % matter)
    tfine = linspace(min(t),max(t),1000);
    for ii=1:size(pw_genes,1)
        list_tmp = pw_genes{ii,1};
        for jj=1:size(list_tmp,1)
            tmp_name = list_tmp{jj,1};
            ind = strcmp(tmp_name,res(:,1));
            if sum(ind)>0
                fun_type = res{ind,2};
                beta = res{ind,3};
                reg_interval = [0 48];
                if (fun_type==1)||(fun_type==3)
                    reg_interval_coeff = beta(3);
                    amplitude_coeff = beta(1); 
                    % calculate maximal amplitude within 48 hours
                    if fun_type==1 % sigmoidal increase
                        b = [beta(1:2) log(beta(3))];
                        y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3))));
                        amplitude = max(y_mod); 
                        ind = find(y_mod>.5*amplitude,1,'first');
                        resp_time = exp(tfine(ind));
                        reg_interval = [0 resp_time];
                    elseif fun_type==3 % sigmoidal decrease
                        b = [beta(1:2) log(beta(3))];
                        y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3))));
                        amplitude = - min(y_mod); 
                        ind = find(y_mod<-.5*amplitude,1,'first');
                        resp_time = exp(tfine(ind));
                        reg_interval = [0 resp_time];
                    end            
                elseif (fun_type==2)||(fun_type==4)
                    reg_interval_coeff = [beta(3) exp(log(beta(3))+log(beta(5)))];
                    amplitude_coeff = max(beta([1 4])); % maximum of amplitudes                              
                    if fun_type==2 % sigmoidal increase before sigmoidal decrease
                        b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
                        y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3)))) - ...
                            b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));
                        amplitude = max(y_mod) - min(y_mod); 
                        ind = find(y_mod>.5*max(y_mod),1,'first');
                        resp_time_1 = exp(tfine(ind));
                        ind = find(y_mod>(max(y_mod) + y_mod(end))*.5,1,'last');
                        resp_time_2 = exp(tfine(ind));
                        reg_interval = [resp_time_1 resp_time_2];
                    elseif fun_type==4 % sigmoidal decrease before increase        
                        b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
                        y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3)))) + ...
                            b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));   
                        amplitude = max(y_mod) - min(y_mod);
                        ind = find(y_mod<.5*min(y_mod),1,'first');
                        resp_time_1 = exp(tfine(ind));
                        ind = find(y_mod<(min(y_mod) + y_mod(end))*.5,1,'last');
                        resp_time_2 = exp(tfine(ind));
                        reg_interval = [resp_time_1 resp_time_2];
                    end
                end

                res_table = cat(1,res_table,[{ii} {jj} tmp_name fun_type reg_interval amplitude {beta} numel(list_tmp)]);
            end
        end
    end
    save('Results_table_profiles_interpolation','res_table')
else
    load('Results_table_profiles_interpolation.mat')
end

%% create line plots for KEGG pathways and extract inflection points
alpha = 1; 
lin_or_log = 2; 
show_plots = 1; 
save_fig = 1;
if task_create_line_plots == true
    res_lineplots_resptimes = cell(size(pw_list,1),6);
    for ii=1:size(pw_list,1)
        res_lineplots_resptimes(ii,:) = fun_line_plot_KEGG_pw(ii,alpha,pw_list,res_table,lin_or_log,show_plots,save_fig);
    end
    save(['Results_table_lineplots_resptimes_interpolation_alpha_' strrep(num2str(alpha),'.','_')],'res_lineplots_resptimes')
else
    load(['Results_table_lineplots_resptimes_interpolation_alpha_' strrep(num2str(alpha),'.','_') '.mat'])
end
%% separately create figure for 'Transcription factors'
fun_line_plot_KEGG_pw_sep_fun_types(354,alpha,pw_list,res_table,...
    all_genes,res_l2fc_corr,lin_or_log,task_print_images);

%% create scatter plot of up- or downregulated pathways (total number or fraction)
xy = []; xy_names = {}; N_min = 10; m_color = .7*ones(1,3);
for ii=1:size(res_lineplots_resptimes,1)
    if ~isempty(res_lineplots_resptimes{ii,4})&&(res_lineplots_resptimes{ii,4}(1,3)>=N_min)
       xy = cat(1,xy,res_lineplots_resptimes{ii,4}); 
       xy_names = cat(1,xy_names,res_lineplots_resptimes{ii,1});
    end
end

% marker size for scatter plots
ms = 7;
% display zero entries separately
factor_sep = 10.^[-.5 -.25];
xy_disp = xy;
xy_disp(xy(:,1)==0,1)=factor_sep(1);
xy_disp(xy(:,2)==0,2)=factor_sep(1);
xy_perc = 100*bsxfun(@times,xy(:,1:2),1./xy(:,3));
xy_perc_disp = xy_perc;
xy_perc_min = .33;
xy_perc_disp(xy_perc(:,1)==0,1)=xy_perc_min*factor_sep(1);
xy_perc_disp(xy_perc(:,2)==0,2)=xy_perc_min*factor_sep(1);

% only include pathways with more than 10 regulated genes
hf = figure(100); clf(hf); set(hf,'Position',[50 50 600 450]); 
hold on;
plot(xy_disp(:,2),xy_disp(:,1),'o','MarkerFaceColor',m_color,'MarkerEdgeColor','k','MarkerSize',ms)
set(gca,'XScale','log','YScale','log')
xyl = get(gca,'YLim').*[0 3]+[.2 0];
plot(ones(2,1)*factor_sep(2),xyl,'k--')
plot(xyl,ones(2,1)*factor_sep(2),'k--')
plot(xyl,xyl,'k:')
xyt = [factor_sep(1) 1 10 100];
xytl = {'0' '1' '10' '100'};
set(gca,'LineWidth',1,'box','on','TickDir','out','XLim',xyl,'XTick',xyt,'XTickLabel',xytl,...
    'YLim',xyl,'YTick',xyt,'YTickLabel',xytl,'FontSize',fs,'XMinorTick','off','YMinorTick','off')
xlabel('Number of downregulated genes','FontSize',fs-5)
ylabel('Number of upregulated genes','FontSize',fs-5)
title([sprintf('Number of regulated genes\n') 'with log2 f.c.\geq1'],'FontSize',fs-5)

% for labelling selected pathways
[~,ord]=sort(xy(:,1).*xy(:,2),'descend');
top_diag_xy = [xy_names(ord),num2cell(xy(ord,:))];
annotation('arrow',[.645 .7183],[.76 .76]);
text(1.742,108.458,'Transcription factors','FontSize',fs-9)
annotation('arrow',[.75 .75],[.8067 .7489]);
text(10.7,210.2,'Membrane trafficking','FontSize',fs-9)
annotation('arrow',[.6567 .704],[.722 .722]);
text(0.607,69.301,'Chromosome and assoc. prot.','FontSize',fs-9)
annotation('arrow',[.815 .7682],[.5911 .6556]);
text(83.361,14.011,'Exosome','FontSize',fs-9)
annotation('arrow',[.7433 .7433],[.433 .6067]);
text(32.66,1.945,sprintf('Peptidases\nand inhibitors'),'FontSize',fs-9)

hf = figure(101); clf(hf); set(hf,'Position',[550 50 600 450]);
hold on;
plot(xy_perc_disp(:,2),xy_perc_disp(:,1),'o','MarkerFaceColor',m_color,'MarkerEdgeColor','k','MarkerSize',ms)
set(gca,'XScale','log','YScale','log')
xyl = [.06,180];
xyt = [xy_perc_min*factor_sep(1) 1 3 10 30];
xytl = {'0' '1' '3' '10' '30'};
plot(ones(2,1)*xy_perc_min*factor_sep(2),xyl,'k--')
plot(xyl,ones(2,1)*xy_perc_min*factor_sep(2),'k--')
plot(xyl,xyl,'k:')
set(gca,'LineWidth',1,'box','on','TickDir','out','XLim',xyl,'XTick',xyt,'XTickLabel',xytl,...
    'YLim',xyl,'YTick',xyt,'YTickLabel',xytl,'FontSize',fs,'XMinorTick','off','YMinorTick','off')
xlabel('Downregulated genes, %','FontSize',fs-5)
ylabel('Upregulated genes, %','FontSize',fs-5)
title([sprintf('Percentage of regulated genes\n') ' with log2 f.c.\geq1'],'FontSize',fs-5)

% for labelling selected pathways
[~,ord]=sort(xy_perc_disp(:,1).*xy_perc_disp(:,2),'descend');
top_diag_xy_perc = [xy_names(ord),num2cell(xy_perc_disp(ord,:)) num2cell(xy(ord,3))];

text(168,37.7,sprintf('Cytochrome\nP450'),'HorizontalAlignment','right','FontSize',fs-9)
annotation('arrow',[.79 .775],[.7222 .6267]);

text(13.0273,91.95,'Retinol metabolism','FontSize',fs-9)
annotation('arrow',[.7117 .7583],[.7778 .6467]);

text(0.78,62.8,sprintf('Ribosome'),'HorizontalAlignment','right','FontSize',fs-9)
annotation('arrow',[0.335,0.375],[0.744,0.707]);

text(0.221,128.1,sprintf('Cytokines and growth factors'),'FontSize',fs-9)
annotation('arrow',[0.572,0.653],[0.804,0.711]);

text(1.1,75,sprintf('COVID-19'),'FontSize',fs-9)
annotation('arrow',[0.51,0.572],[0.762,0.698]);

%% rank pathways according to percentage of types 1 and 2:
Nl = size(res_lineplots_resptimes,1);
mat_t12_t34_Ng = zeros(Nl,3);
Tab_Nup_Ndown_Nexp_Ntot = zeros(Nl,4);
% include only pathways with >=Nmin_pw_map regulated genes
for ii=1:Nl
    tmp = res_lineplots_resptimes{ii,4};
    % find total number of genes/expressed genes in pw
    ind = find(cell2mat(res_table(:,1))==ii,1,'first');
    if (~isempty(ind))&&(~isempty(tmp))
        N_tot = res_table{ind,end};
        Tab_Nup_Ndown_Nexp_Ntot(ii,:) = [tmp N_tot];
        if (~isempty(tmp))&&(sum(tmp(1:2))>=Nmin_pw_map)
            mat_t12_t34_Ng(ii,:) = tmp;
        end
    end
end

% rank according to percentage
perc_t12 = mat_t12_t34_Ng(:,1)./mat_t12_t34_Ng(:,3);
perc_t12(isnan(perc_t12),1)=0;
[perc_t12_ord,ord_pos] = sort(perc_t12,'descend');
res_lineplots_resptimes_ord_pos = res_lineplots_resptimes(ord_pos,:);

% additionally highlight N_best up-regulated pathways in Figs. 4A and B:
xy_N_best = []; xy_names_N_best = {}; 
res_N_best = res_lineplots_resptimes_ord_pos(1:(N_best+1),:);
for ii=1:size(res_N_best,1)
       xy_N_best = cat(1,xy_N_best,res_lineplots_resptimes_ord_pos{ii,4}); 
       xy_names_N_best = cat(1,xy_names_N_best,res_lineplots_resptimes_ord_pos{ii,1});
end
xy_perc_disp = 100*bsxfun(@times,xy_N_best(:,1:2),1./xy_N_best(:,3));

figure(100) % Fig. 4A
% display zero entries separately
xy_N_best(xy_N_best(:,1)==0,1)=factor_sep(1);
xy_N_best(xy_N_best(:,2)==0,2)=factor_sep(1);
plot(xy_N_best(:,2),xy_N_best(:,1),'o','MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','k','MarkerSize',ms)

figure(101) % Fig. 4B
% display zero entries separately
xy_perc_disp(xy_N_best(:,1)==0,1)=xy_perc_min*factor_sep(1);
xy_perc_disp(xy_N_best(:,2)==0,2)=xy_perc_min*factor_sep(1);
plot(xy_perc_disp(:,2),xy_perc_disp(:,1),'o','MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','k','MarkerSize',ms)

ind_half = zeros((N_best+1),2); mat_tmp = [];
for ii=1:(N_best+1)
    tmp_line = res_lineplots_resptimes_ord_pos{ii,2};
    mat_tmp = cat(1,mat_tmp,tmp_line/max(tmp_line));
    ind_half(ii,:) = fun_find_first_and_second_thd(tmp_line/max(tmp_line));
end

Tab_Nup_Ndown_Nexp_Ntot_ord = Tab_Nup_Ndown_Nexp_Ntot(ord_pos(1:(N_best+1)),:);
pw_names = res_lineplots_resptimes_ord_pos(1:(N_best+1),1);

% display shorter name
ind = find(strcmp(pw_names,'Viral protein interaction with cytokine and cytokine receptor'));
if ~isempty(ind)
    pw_names{ind,1}='Viral prot. interaction cytokine/cytokine rec.';
end

% re-order according to time points
[~,ord_pw] = sort(ind_half(:,1),'descend');
ind_half_reord=ind_half(ord_pw,:);
mat_tmp_reord = mat_tmp(ord_pw,:)*100;
pw_names_reord = pw_names(ord_pw,1);
Tab_num_sort = Tab_Nup_Ndown_Nexp_Ntot_ord(ord_pw,:);

% remove double entry 'Ribosome'
ind = find(strcmp(pw_names_reord,'Ribosome'));
num_ind = Tab_num_sort(ind,3);
ind_keep = setdiff(1:(N_best+1),ind(num_ind==max(num_ind)));
ind_half_reord=ind_half_reord(ind_keep,:);
mat_tmp_reord = mat_tmp_reord(ind_keep,:);
pw_names_reord = pw_names_reord(ind_keep,1);
Tab_num_sort = Tab_num_sort(ind_keep,:);

hf = figure(1); clf(hf); 
set(hf,'Position',[50 50 1200 38*N_best])

n_per_h = 10;
arlen = 48*n_per_h;
tmin = .5; tmax = 48;
fscale = @(t,tmin,tmax,arlen) 1 + max(round((arlen-1) * (log(t)-log(tmin))/(log(tmax)-log(tmin))),0);
xt = fscale([.5 1 2 4 7 12 24 48],tmin,tmax,arlen);

rbcm = redbluecmap;
rbcm_up= [interp1((1:6)',rbcm(6:end,1),(1:.01:6)'), interp1((1:6)',rbcm(6:end,2),(1:.01:6)'), ...
    interp1((1:6)',rbcm(6:end,3),(1:.01:6)')];
rbcm_down= flipud([interp1((1:6)',rbcm(1:6,1),(1:.01:6)'), interp1((1:6)',rbcm(1:6,2),(1:.01:6)'), ...
    interp1((1:6)',rbcm(1:6,3),(1:.01:6)')]);

xlabels = {'0.5' '1' '2' '4' '7' '12' '24' '48'};
imagesc(mat_tmp_reord)
set(gca,'box','on','TickDir','out','YTick',0:N_best,'YTickLabel','',...
    'YLim',[.5 N_best+.5],'XTick',xt,'XTickLabel',xlabels,'FontSize',fs,'XLim',[.5 arlen+.5])
colormap(gca,rbcm_up); hold on;

% display time points of 50% regulation
for jj=1:N_best
    plot(ind_half_reord(jj,1)*ones(2,1),jj+[-.5 .5],'Color','y','LineWidth',3) 
    if ~isnan(ind_half_reord(jj,2))
        plot(ind_half_reord(jj,2)*ones(2,1),jj+[-.5 .5],'Color','y','LineWidth',3) 
    end
end

% legend
text('FontSize',17,'String','Time point at which 50% of genes are up/downregulated','Position',[576 29.9]);
annotation('rectangle',[0.5166 0.8684 0.0042 0.0212],'FaceColor',[1 1 0]);
text('FontSize',17,'String','Upregulation of mitochondrial genes',...
    'Position',[576 28.81]);
annotation('rectangle',[0.5058 0.8457 0.0158 0.016],'Color',[0 0.6 0],'FaceColor',[0 0.6 0]);
text('FontSize',17,'String','Upregulation of genes involved in ribosomes',...
    'Position',[576 27.72]);
annotation('rectangle',[0.5058 0.8230 0.0158 0.016],'Color',[.2 .2 .9],'FaceColor',[.2 .2 .9]);

% table, gene counts
text(1.9*arlen,1.05*N_best,'(N_{up}, N_{expr}, N_{tot})','FontSize',fs-10)

green_labels = {'Cardiac muscle contraction';'Non-alcoholic fatty liver disease';'Prion disease';'Parkinson disease';...
    'Oxidative phosphorylation';'Retrograde endocannabinoid signaling';'Thermogenesis'};
blue_labels = {'Ribosome'; 'Coronavirus disease - COVID-19'};

for jj=1:N_best
    if sum(strcmp(green_labels,pw_names_reord{jj,1}))>0
        text(1.02*arlen,jj,pw_names_reord{jj,1},'FontSize',fs-5,'HorizontalAlignment','left','Color',[0 0.6 0]) 
    elseif sum(strcmp(blue_labels,pw_names_reord{jj,1}))>0
        text(1.02*arlen,jj,pw_names_reord{jj,1},'FontSize',fs-5,'HorizontalAlignment','left','Color',[.2 .2 .9]) 
    else
        text(1.02*arlen,jj,pw_names_reord{jj,1},'FontSize',fs-5,'HorizontalAlignment','left')         
    end
    text(2.03*arlen,jj,['(' num2str(Tab_num_sort(jj,1)) ', ' num2str(Tab_num_sort(jj,3)) ', '...
        num2str(Tab_num_sort(jj,4)) ')'],'FontSize',fs-10,'HorizontalAlignment','center') 
end
xlabel('Time in hours','FontSize',fs)
hc = colorbar('FontSize',fs,'Location','northoutside');
set(hc,'YTick',0:25:100,'YLim',[0 100],'Position',[.21 .8 .23 .03],'FontSize',fs-5)
zlab = get(hc,'ylabel');
set(zlab,'String','Percentage of upregulated genes','FontSize',fs-5)

set(gca,'Position',[.05 .22 .4 .56]); axis xy;     
im_name = ['Fig_pws_incr_genes_top' num2str(N_best) '_perc_interp'];        

if task_print_images==true
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r200',im_name)
    set(gcf,'PaperSize',[45 30]);
    print(im_name,'-dpdf')
end   

%% rank pathways according to percentage of types 3 and 4 profiles (cont. or transient decrease):
perc_t34 = mat_t12_t34_Ng(:,2)./mat_t12_t34_Ng(:,3);
perc_t34(isnan(perc_t34),1)=0;
[perc_t34_ord,ord_neg] = sort(perc_t34,'descend');

mat_tmp = [];
res_lineplots_resptimes_ord_neg = res_lineplots_resptimes(ord_neg,:);
Tab_Nup_Ndown_Nexp_Ntot_ord = Tab_Nup_Ndown_Nexp_Ntot(ord_neg(1:N_best),:);
pw_names = res_lineplots_resptimes_ord_neg(1:N_best,1);

% additionally highlight N_best down-regulated pathways in Figs. 4A and B:
xy_N_best = []; xy_names_N_best = {}; 
res_N_best = res_lineplots_resptimes_ord_neg(1:(N_best+1),:);
for ii=1:size(res_N_best,1)
       xy_N_best = cat(1,xy_N_best,res_lineplots_resptimes_ord_neg{ii,4}); 
       xy_names_N_best = cat(1,xy_names_N_best,res_lineplots_resptimes_ord_neg{ii,1});
end
xy_perc_disp = 100*bsxfun(@times,xy_N_best(:,1:2),1./xy_N_best(:,3));

figure(100) % Fig. 4A
% display zero entries separately
xy_N_best(xy_N_best(:,1)==0,1)=factor_sep(1);
xy_N_best(xy_N_best(:,2)==0,2)=factor_sep(1);
plot(xy_N_best(:,2),xy_N_best(:,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k','MarkerSize',ms)
if task_print_images==true
    im_name = 'Fig_scatter_KEGG_pw_number';
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r200',im_name)
    set(gcf,'PaperSize',[30 20]);
    print(im_name,'-dpdf')
end  

figure(101) % Fig. 4B
% display zero entries separately
xy_perc_disp(xy_N_best(:,1)==0,1)=xy_perc_min*factor_sep(1);
xy_perc_disp(xy_N_best(:,2)==0,2)=xy_perc_min*factor_sep(1);
plot(xy_perc_disp(:,2),xy_perc_disp(:,1),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k','MarkerSize',ms)
if task_print_images==true
    im_name = 'Fig_scatter_KEGG_pw_percentage';
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r200',im_name)
    set(gcf,'PaperSize',[30 20]);
    print(im_name,'-dpdf')
end  

% display shorter name
ind = find(strcmp(pw_names,'Viral protein interaction with cytokine and cytokine receptor'));
if ~isempty(ind)
    pw_names{ind,1}='Viral prot. interaction cytokine/cytokine rec.';
end
ind_half = zeros(N_best,2);
for ii=1:N_best
    tmp_line = res_lineplots_resptimes_ord_neg{ii,3};
    tmp_scaled = tmp_line/max(tmp_line);
    mat_tmp = cat(1,mat_tmp,tmp_scaled);
    ind_half(ii,:) = fun_find_first_and_second_thd(1-tmp_scaled);
end

% re-order according to time points
[~,ord_pw] = sort(ind_half(:,1),'descend');
ind_half_reord=ind_half(ord_pw,:);
mat_tmp_reord = mat_tmp(ord_pw,:);
mat_tmp_reord_disp = (1-mat_tmp_reord)*100;
pw_names_reord = pw_names(ord_pw,1);
Tab_num_sort = Tab_Nup_Ndown_Nexp_Ntot_ord(ord_pw,:);

hf = figure(2); clf(hf); 
set(hf,'Position',[50 50 1200 38*N_best])

n_per_h = 10;
arlen = 48*n_per_h;
tmin = .5; tmax = 48;
fscale = @(t,tmin,tmax,arlen) 1 + max(round((arlen-1) * (log(t)-log(tmin))/(log(tmax)-log(tmin))),0);
xt = fscale([.5 1 2 4 7 12 24 48],tmin,tmax,arlen);
        
xlabels = {'0.5' '1' '2' '4' '7' '12' '24' '48'};
imagesc(mat_tmp_reord_disp)
set(gca,'box','on','TickDir','out','YTick',0:N_best,'YTickLabel','',...
    'YLim',[.5 N_best+.5],'XTick',xt,'XTickLabel',xlabels,'FontSize',fs,'XLim',[1 arlen])
colormap(gca,rbcm_down); hold on;
% display time points of 50% regulation
for jj=1:N_best
    plot(ind_half_reord(jj,1)*ones(2,1),jj+[-.5 .5],'Color','y','LineWidth',3) 
    if ~isnan(ind_half_reord(jj,2))
        plot(ind_half_reord(jj,2)*ones(2,1),jj+[-.5 .5],'Color','y','LineWidth',3) 
    end
end

% legend
text('FontSize',17,'String','Time point at which 50% of genes are up/downregulated','Position',[621 29.9]);
annotation('rectangle',[0.5533 0.8684 0.0042 0.0212],'FaceColor',[1 1 0]);

text(2*arlen,1.06*N_best,'(N_{down},N_{expr},N_{tot})','FontSize',fs-10)
for jj=1:N_best
    text(1.02*arlen,jj,pw_names_reord{jj,1},'FontSize',fs-5,'HorizontalAlignment','left') 
    text(2.13*arlen,jj,['(' num2str(Tab_num_sort(jj,2)) ', ' num2str(Tab_num_sort(jj,3)) ', '...
        num2str(Tab_num_sort(jj,4)) ')'],'FontSize',fs-10,'HorizontalAlignment','center') 
end
xlabel('Time in hours','FontSize',fs)

hc = colorbar('FontSize',fs,'Location','northoutside');
set(hc,'YTick',0:25:100,'YLim',[0 100],'Position',[.21 .8 .23 .03],'FontSize',fs-5)
zlab = get(hc,'ylabel');
set(zlab,'String','Percentage of downregulated genes','FontSize',fs-5)

set(gca,'Position',[.05 .22 .4 .56]); axis xy; 
im_name = ['Fig_pws_decr_genes_top' num2str(N_best) '_perc_interp'];        
if task_print_images==true
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r200',im_name)
    set(gcf,'PaperSize',[45 30]);
    print(im_name,'-dpdf')
end   

%% create line plot for GO term "defense response to virus", GO:0051607
D2 = readcell('QuickGO_defense_response_to_virus_GO0051607.txt');
genes_antivir = upper(D2(2:end,3));
g_list = unique(genes_antivir);
res_list = {};
t = log([.5 1 2 4 7 12 24 48])'; % log-scaled time points
                                 % (due to the parameter limits, the
                                 % time of the first point doen't
                                 % matter)
tfine = linspace(min(t),max(t),1000);
antivir_profiles = [];
lin_or_log = 2;
alpha = log2(1.5); 
for jj=1:size(g_list,1)
    tmp_name = g_list{jj,1};
    ind = strcmp(tmp_name,res(:,1));
    if sum(ind)>0
        fun_type = res{ind,2};
        beta = res{ind,3};
        reg_interval = [0 48];
        if (fun_type==1)||(fun_type==3)
            reg_interval_coeff = beta(3);
            amplitude_coeff = beta(1); 
            % calculate maximal amplitude within 48 hours
            if fun_type==1 % sigmoidal increase
                b = [beta(1:2) log(beta(3))];
                y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3))));
                y_prof = b(1)./(1 + exp(-b(2)*(t-b(3))));
                amplitude = max(y_mod); 
                if amplitude>=alpha
                    antivir_profiles = cat(1,antivir_profiles,y_prof');
                end
                ind = find(y_mod>.5*amplitude,1,'first');
                resp_time = exp(tfine(ind));
                reg_interval = [0 resp_time];
            elseif fun_type==3 % sigmoidal decrease
                b = [beta(1:2) log(beta(3))];
                y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3))));
                amplitude = - min(y_mod); 
                ind = find(y_mod<-.5*amplitude,1,'first');
                resp_time = exp(tfine(ind));
                reg_interval = [0 resp_time];
            end            
        elseif (fun_type==2)||(fun_type==4)
            reg_interval_coeff = [beta(3) exp(log(beta(3))+log(beta(5)))];
            amplitude_coeff = max(beta([1 4])); % maximum of amplitudes                              
            if fun_type==2 % sigmoidal increase before sigmoidal decrease
                b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
                y_mod = b(1)./(1 + exp(-b(2)*(tfine-b(3)))) - ...
                    b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));
                amplitude = max(y_mod) - min(y_mod); 
                y_prof = b(1)./(1 + exp(-b(2)*(t-b(3)))) - ...
                    b(4)./(1 + exp(-b(2)*(t-(b(3)+b(5)))));
                if amplitude>=alpha
                    antivir_profiles = cat(1,antivir_profiles,y_prof');
                end
                ind = find(y_mod>.5*max(y_mod),1,'first');
                resp_time_1 = exp(tfine(ind));
                ind = find(y_mod>(max(y_mod) + y_mod(end))*.5,1,'last');
                resp_time_2 = exp(tfine(ind));
                reg_interval = [resp_time_1 resp_time_2];
            elseif fun_type==4 % sigmoidal decrease before increase        
                b = [beta(1:2) log(beta(3)) beta(4) log(beta(5))];
                y_mod = -b(1)./(1 + exp(-b(2)*(tfine-b(3)))) + ...
                    b(4)./(1 + exp(-b(2)*(tfine-(b(3)+b(5)))));  
                amplitude = max(y_mod) - min(y_mod);
                y_prof = -b(1)./(1 + exp(-b(2)*(t-b(3)))) + ...
                    b(4)./(1 + exp(-b(2)*(t-(b(3)+b(5)))));  
                ind = find(y_mod<.5*min(y_mod),1,'first');
                resp_time_1 = exp(tfine(ind));
                ind = find(y_mod<(min(y_mod) + y_mod(end))*.5,1,'last');
                resp_time_2 = exp(tfine(ind));
                reg_interval = [resp_time_1 resp_time_2];
            end
        end
        res_list = cat(1,res_list,[{1} {jj} tmp_name fun_type reg_interval amplitude {beta} numel(g_list)]);
    end
end

fun_line_plot_KEGG_pw_sep_fun_types(1,alpha,{'Defense response to virus'},res_list,...
     all_genes,res_l2fc_corr,lin_or_log,task_print_images);

pos_prof_norm = bsxfun(@times,antivir_profiles,1./max(antivir_profiles,[],2));
m_pos_profiles = mean(pos_prof_norm,1);
sem_pos_profiles = std(pos_prof_norm,[],1)/sqrt(size(pos_prof_norm,1));
data_for_fitting = [exp(t),m_pos_profiles'/max(m_pos_profiles),sem_pos_profiles'/max(m_pos_profiles)];

hf = figure(102); clf(hf); set(hf,'Position',[50 50 800 500]); hold on;
errorbar(t,m_pos_profiles/max(m_pos_profiles),sem_pos_profiles/max(m_pos_profiles),'o','LineWidth',1.5)
set(gca,'TickDir','out','LineWidth',1,'box','on','FontSize',fs-10,'XTick',t,'XTickLabel',xlabels,'YLim',[-.1 1.4])
xlabel('Time in hours'); grid on;
ylabel('Average profile of upregulated genes')
title('Upregulated antiviral response genes')

im_name = 'Fig_Defense_response_to_virus_LinePlot';
if task_print_images==true
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r200',im_name)
    set(gcf,'PaperSize',[45 30]);
    print(im_name,'-dpdf')
end 

