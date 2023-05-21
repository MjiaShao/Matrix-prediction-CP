clear
font_size = 20;
MarkerSize = 15;  LineWidth = 1.2;
Methodname = ["Our method","softImpute"];
datasetname = ["SZ","NC"];
misproall = [300, 1000];
typemissing = "maxmissing";
m = 246;
linewidth1 = 1.2;
Nodeseq = 1:246;
colors = [	[1 0 0]; 	[0 0 1]; 	[0.3010, 0.7450, 0.9330]* 0.8;  	[0, 0.5, 0]; [0.4940, 0.1840, 0.5560];   	[0.9290, 0.6940, 0.1250]* 0.8];
fig = figure('Visible','on','Position',[100,100,1160,650],'Units','centimeters');      
tlo = tiledlayout(2,4,'TileSpacing','none','Padding','compact');  
timeaveour = [];
timeavesoft = [];
for inx = 1:length(misproall)
    misper = misproall(inx);
    switch misper
        case 300
            powerq = 3;
        case 1000
            powerq = 3;
    end
    outputall = [];
for k = 1:length(datasetname)
    dataname = datasetname(k);
    switch dataname 
        case "SZ"
            totalboot = 103;
            labelname = "Patients";
            widthbox = 130;
        case "NC"
            totalboot = 124;
            labelname = "Healthy";
           widthbox = 122;
    end
    coversofttable = [];
    CIlensofttable = [];
    timesofttable = [];
    coverourtable = [];
    CIlenourtable = [];
    timeourtable = [];
    for seed = 1:8
    load(sprintf("./result/data_ourmethod_%s_random_mislarge_%d_%s_%d_newmethod_%d_svd",dataname,misper, typemissing,seed,powerq));
    coversoftall = fscanf(fopen(sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_cover.txt",dataname,misper, typemissing, seed),'r'),'%f',[totalboot,m])';
    CIlensoftall =  fscanf(fopen(sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_CIlen.txt",dataname,misper, typemissing, seed),'r'),'%f',[totalboot,m])';
    timesoftall =  fscanf(fopen(sprintf("./result/data_softImpute_%s_random_mislarge_%d_%s_%d_timecost.txt",dataname,misper, typemissing, seed),'r'),'%f',[totalboot,m])';
    coversofttable = [coversofttable, coversoftall];
    CIlensofttable = [CIlensofttable, CIlensoftall];
    timesofttable = [timesofttable, timesoftall];
    coverourtable = [coverourtable, coverall];
    CIlenourtable = [CIlenourtable, CIlenall];
    timeourtable = [timeourtable, timeall];
    end
    ax = nexttile((inx-1)*length(datasetname)+k);
    plot1 = shadedErrorBar(1:length(Nodeseq),mean(coverourtable,2),sqrt(mean(coverourtable,2).*(1-mean(coverourtable,2))/length(mean(coverourtable,2))),{'color',colors(1,:),'linewidth',linewidth1},1);
    hold on 
    plot2 = shadedErrorBar(1:length(Nodeseq),mean(coversofttable,2),sqrt(mean(coversofttable,2).*(1-mean(coversofttable,2))/length(mean(coversofttable,2))),{'color',colors(4,:),'linewidth',linewidth1},1);
    hold off
    xlim([1-0.1,246+0.1])
    ylim([0.7,1])
    line(xlim,[0.9 0.9], 'Color', 'black', 'LineWidth', 2);
    text(4,0.72,labelname,'fontsize',font_size+6,'interpreter','latex')
    hold on 
    rectangle('position',[1+0.2 0.7 widthbox 0.05],'LineWidth', 1.2)
    hold off 
    title(sprintf('Miss. = %d',misper),'interpreter','latex');

    set(gca,'fontsize',font_size);
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    pbaspect([1.1 1 1])
    if inx == 1 && k  ==1
        ylabel('Cover. prob.','fontsize',font_size+8);
    else
        yticklabels([]);
        ylabel('');
    end
     xticklabels([]);


    ax = nexttile((inx-1)*length(datasetname)+k+4*1);
    plot1 = shadedErrorBar(1:length(Nodeseq),nanmean(CIlenourtable,2),nanstd(CIlenourtable,0,2)/sqrt(length(mean(CIlenourtable,2))),{'color',colors(1,:),'linewidth',linewidth1},1);	
    hold on
    plot2 = shadedErrorBar(1:length(Nodeseq),nanmean(CIlensofttable,2),nanstd(CIlensofttable,0,2)/sqrt(length(mean(CIlensofttable,2))),{'color',colors(4,:),'linewidth',linewidth1},1);
    hold off
    xlim([1-0.1,246+0.1])
    ylim([0.2,1.2])

    set(gca, 'FontSize', font_size)
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    pbaspect([1.1 1 1])
    if inx == 1 && k  ==1
        ylabel('CI length','fontsize',font_size+8);
    else
        yticklabels([]); 
        ylabel('');
    end
    if inx == 1 && k == 1
    legend([plot1.mainLine,plot2.mainLine],'Our Alg.1','softImpute','location','northeast','fontsize', font_size+4); 
    end
    xlabel('Node index', 'FontSize', font_size+12);
    xlim([1-0.1,246+0.1])
    xticks([100,200])
    xticklabels({'100','200'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',font_size+2)

timeaveour(inx,k) = mean(timeourtable(:));
timeavesoft(inx,k) = mean(timesofttable(:));
timestdour(inx,k) = std(timeourtable(:));
timestdsoft(inx,k) = std(timesofttable(:));
end
end
saveas(fig,sprintf("./figure/data_missing_comparison.png"))
