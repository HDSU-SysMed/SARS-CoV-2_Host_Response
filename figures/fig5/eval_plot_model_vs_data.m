%% Plot model vs. data, simulate model species (Fig. 5B and C), 
% simulate changes in V_0 (fig. S10)

global k % model parameters

task_print_figures = false;
fs = 18; lw = 1.5;

% results of 0.2% (n=10) best of N=5000 fits
load('Res_model_fitting.mat')

%% perform model simulations, assess effect of changing in initial amount of virus V_0
% create fig. S10
hf = figure(100); clf(hf);
set(hf,'Position',[100,100,1000,700])
xIDs = {'V','P','M','A'}; titles = {'Virus RNA, V', 'Virus protein, P',...
    'mRNA of anti-viral proteins, M','Anti-viral proteins, A'};
nsteps = 5;
colors = redbluecmap(nsteps);
colors(round(nsteps/2),:)=[0 0 0];
factors = logspace(-2,2,nsteps);
lt = {'V_0/100','V_0/10','V_0','10\cdotV_0','100\cdotV_0'};

k = cell2mat(k_s{1,1}(:,2));
x0 = zeros(1,4);
for jj=1:nsteps   
    x0(1) =   k(end) * factors(jj); % V_0
    [t,x] = ode15s(@sim_SARS_CoV2_model_ODEs,[0 48], x0);
    hold on;
    n = ceil(sqrt(length(x0)));
    for ii=1:4
        subplot(2,2,ii); hold on;
        h = plot(t,x(:,ii),'Color',colors(jj,:),'LineWidth',lw);
    end
end
positions = [.13 .63 .28 .28; .51 .63 .28 .28;    
    .13 .15 .28 .28; .51 .15 .28 .28];
ylu = [1e2,1.1,1.1,1.1]; 

for ii=1:4
    subplot(2,2,ii);
    title(titles{ii},'FontSize',fs)
    yl = get(gca,'YLim');
    set(gca,'FontSize',fs,'box','on','LineWidth',1,'TickDir','out',...
        'XTick',0:12:48,'XLim',[0 48],'YLim',[yl(1) ylu(ii)])
    grid on
    if ii==1
       set(gca,'YScale','log','YTick',logspace(-6,2,5))
    end
    if ii==4
        hl = legend(lt,'Location','EastOutside','FontSize',fs);
        legend('boxoff');  hl.Position = [0.8395    0.3979    0.0978    0.2593];
    end
    set(gca,'Position',positions(ii,:))
    ylabel([xIDs{ii} ' in a. u.'])
    xlabel('Time in hours')
end
if task_print_figures==1
    im_name = 'Fig_S10_sim_changes_in_V0';
    set(gcf,'PaperPositionMode','auto')
    print(im_name,'-dpng','-r300')
    set(gcf,'PaperSize',[30 20]);
    print(im_name,'-dpdf')
end

%% determine 1-sigma C.I. for 10 best fits
y_all = []; y_best = y_tr_all{1,1}{1,1};
x_all = []; x_best = x_tr_all{1,1};
fluxes_all = []; Nbr = 10;
for jj=1:Nbr
    y_all = cat(3,y_all,y_tr_all{jj,1}{1,1});
    x_tmp = x_tr_all{jj,1}; % V P M A
    x_all = cat(3,x_all,x_tmp);
    k = k_s{jj,1};
end
y_s = []; x_s = []; y_m = []; x_m = [];
for ii=1:size(y_all,2)
    tmp = squeeze(y_all(:,ii,:));
    CI_tmp = std(tmp,[],2);
    m_tmp = mean(tmp,2);
    y_s = cat(2,y_s,CI_tmp);
    y_m = cat(2,y_m,m_tmp);
end
for ii=1:size(x_all,2)
    tmp = squeeze(x_all(:,ii,:));
    CI_tmp = std(tmp,[],2);
    m_tmp = mean(tmp,2);
    x_s = cat(2,x_s,CI_tmp);
    x_m = cat(2,x_m,m_tmp);
end

%% create Fig. 5B
xtext = 'Time in hours';
obs_labels = {sprintf('SARS-CoV-2\ntranscripts in a. u.');...
    sprintf('SARS-CoV-2\nproteins in a. u.');...
    sprintf('Anti-viral response\ntranscripts in a. u.')};
lloc = {'NorthEast';'SouthEast';'SouthEast'};
hf = figure(1);
clf(hf)
set(hf,'Position',[100 100 1300 320]);
plot_data = y_tr_all{1,1};

positions = [.1 .277 .189 .6578; .38 .277 .189 .6578; .66 .277 .189 .6578];
ylims = [-.05 0.55;-.1 1.4; -.1 1.4]; 
for ii=1:3
    subplot(1,3,ii)
    hold on
    xtmp = plot_data{1,1}(:,1);
    ymtmp= plot_data{1,1}(:,1+ii);
    xdtmp= plot_data{1,2}(:,1); 
    ydtmp= plot_data{1,2}(:,1+ii);
    stmp = plot_data{1,2}(:,4+ii);

    ind_tmax = find(~isnan(ydtmp),1,'last');
    set(gca,'FontSize',fs,'box','on','LineWidth',1,'TickDir','out','XLim',...
        [0 50],'XTick',0:12:48,'YLim',ylims(ii,:))
    ind = ~isnan(ydtmp);
    ha = errorbar(xdtmp(ind,1),ydtmp(ind,1),stmp(ind,1),'o','LineWidth',lw,'Color','k');
    plot(xtmp,ymtmp,'LineWidth',lw,'Color',[.7 0 0])
    xlabel(xtext,'FontSize',fs)
    ylabel(obs_labels(ii,1),'FontSize',fs)
    grid on
    if ii==3
        legend({'Exp. data','Model fit'},'Location','EastOutside','FontSize',fs-2)
    end
    set(gca,'Position',positions(ii,:))
end
if task_print_figures==1
    im_name = 'Fig_5B_model_vs_data';
    set(gcf,'PaperPositionMode','auto')
    print(im_name,'-dpng','-r300')
    set(gcf,'PaperSize',[30 20]);
    print(im_name,'-dpdf')
end

%% create Fig. 5C
hf = figure(2);
clf(hf)
set(hf,'Position',[100 100 570 320]);
t = x_tr_all{1,1}(:,1);
x = x_tr_all{1,1}(:,2:5); % V P M A
colors = lines(4);

hold on;
for kk=1:4
    plot(-10,0,'Color',colors(kk,:),'LineWidth',1.5)
end
for kk=1:4
    sf = 1/max(x(:,kk));
    patch_data = [[x_m(:,1);flipud(x_m(:,1))],[x_m(:,1+kk)+x_s(:,1+kk);flipud(x_m(:,1+kk)-x_s(:,1+kk))]*sf];
    hp = patch(patch_data(:,1),patch_data(:,2),...
        colors(kk,:),'EdgeColor','none','FaceAlpha',.2);
    plot(x_m(:,1),x_m(:,1+kk)*sf,'Color',colors(kk,:),'LineWidth',1.5)
end
set(gca,'FontSize',fs,'box','on','TickDir','out','LineWidth',1,...
    'YLim',[-0.05 1.08],'XLim',[0 50],'XTick',0:12:48)
xlabel(xtext,'FontSize',fs)
ylabel(sprintf('Simulation, c. in a. u.'),'FontSize',fs)
legend({sprintf('Virus\ntranscripts'),sprintf('Virus\nprotein'),sprintf('mRNA,\nantivir. prot.'),...
    sprintf('Antiviral\nproteins')},'Location','EastOutside','FontSize',fs-4)
grid on

if task_print_figures==1
    im_name = 'Fig_5C_simul';
    set(gcf,'PaperPositionMode','auto')
    print(im_name,'-dpng','-r300')
    set(gcf,'PaperSize',[30 20]);
    print(im_name,'-dpdf')
end

freude = 1;    