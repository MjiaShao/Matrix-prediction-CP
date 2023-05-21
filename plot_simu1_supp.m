clear
font_size = 20;
MarkerSize = 15;  LineWidth = 2;
con_u= [1:1:9];
Nall = [100,400];
graphonall = ["f1","f2", "f3"];
q = 5;
Methodname = ["Our Alg.1", "Our Alg.2","missMDA", "softImpute"];
fig = figure('Visible','on','Position',[100,100,1150*2-50,680],'Units','centimeters');      
tlo = tiledlayout(2,6,'TileSpacing','none','Padding','compact');  
for i = 2
    for j = 6
        ax = nexttile;
    end
end

colors = [	[1 0 0]; 	[0 0 1]; 	[0.3010, 0.7450, 0.9330]* 0.8;  	[0, 0.5, 0]; [0.4940, 0.1840, 0.5560];   	[0.9290, 0.6940, 0.1250]* 0.8];
markers = ['o','+','*','v','X','s'];
for idgra = 1:length(graphonall)
graphonname = graphonall(idgra);
switch graphonname 
    case "f2"
        graphonoutput = "$f_2$ (high-rank)";
    case "f3"
        graphonoutput = "$f_3$ (non-smooth)";
    case "f1"
        graphonoutput = "$f_1$ (smooth)";
end
switch graphonname 
    case "f2"
        cutoffCI = 5.2;
    case "f3"
        cutoffCI = 4.5;
    case "f1"
        cutoffCI = 5.2;
end
switch graphonname 
    case "f2"
        widthbox = 0.74;
    case "f3"
        widthbox = 0.84;
    case "f1"
        widthbox = 0.62;
end
for inx = 1:length(Nall)

    n = Nall(inx);
    coverall = [];
    CIlenall = [];
    timeall = [];
    coverstdall = [];
    CIlenstdall = [];
    timestdall = [];
    timecombineall1 = [];
    timecombineall2 = [];
    timecombineall3 = [];
    timecombineall6 = [];

    
    for k = 1:length(con_u)
        u = con_u(k);
		load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_5_con.mat",n,u,graphonname))
         coverall(k,1) = nanmean(cover);
        CIlenall(k,1) = nanmean(CIlen);
        timeall(k,1) = nanmean(timerecord);
        CIlenstdall(k,1) = nanstd(CIlen)/sqrt(length(CIlen));
        timestdall(k,1) = nanstd(timerecord)/sqrt(length(timerecord));
        coverstdall(k,1) = sqrt(coverall(k,1)*(1-coverall(k,1))/length(cover));
        timecombineall1 = [timecombineall1,timerecord];
		load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_rand_5_svd_onemissing.mat",n,u,graphonname))        
        coverall(k,6) = nanmean(cover);
        CIlenall(k,6) = nanmean(CIlen);
        timeall(k,6) = nanmean(timerecord);
        CIlenstdall(k,6) = nanstd(CIlen)/sqrt(length(CIlen));
        timestdall(k,6) = nanstd(timerecord)/sqrt(length(timerecord));
        coverstdall(k,6) = sqrt(coverall(k,1)*(1-coverall(k,1))/length(cover));
        timecombineall6 = [timecombineall6,timerecord];
        if n == 100
            data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_5.csv",n,u,graphonname));
           coverall(k,2) = nanmean(data.cover);
            CIlenall(k,2) = nanmean(data.CIlen);
            timeall(k,2) = nanmean(data.timerecord);
            coverstdall(k,2) = sqrt(coverall(k,2)*(1-coverall(k,2))/length(data.cover));
            CIlenstdall(k,2) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            timestdall(k,2) = nanstd(data.timerecord)/sqrt(length(data.timerecord));
            timecombineall2 = [timecombineall2,data.timerecord'];
        elseif n == 400
            coverall(k,2) = NaN;
            CIlenall(k,2) = NaN;
            timeall(k,2)  = NaN;
            coverstdall(k,2) = NaN;
            CIlenstdall(k,2) = NaN;
            timestdall(k,2)  = NaN; 
            timecombineall2  = NaN;
        end
		data = readtable(sprintf("./result/conf_net_soft_%d_%d_%s_con_5.csv",n,u,graphonname));
        coverall(k,3) = nanmean(data.cover);
        CIlenall(k,3) = nanmean(data.CIlen);
        timeall(k,3) = nanmean(data.timerecord);
        coverstdall(k,3) = sqrt(coverall(k,3)*(1-coverall(k,3))/length(data.cover));
        CIlenstdall(k,3) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
        timestdall(k,3) = nanstd(data.timerecord)/sqrt(length(data.timerecord)); 
        timecombineall3 = [timecombineall3,data.timerecord'];
    end
    ax = nexttile((inx-1)*length(graphonall)+idgra);
    box on;
    ax.LineWidth = 2;
    plot1 =  shadedErrorBar(con_u/10,coverall(:,6),coverstdall(:,6),{'color',colors(1, :), 'linewidth',LineWidth,'Marker',markers(1),'MarkerSize',MarkerSize},1);
    hold on 
    plot2 =  shadedErrorBar(con_u/10,coverall(:,1),coverstdall(:,1),{'color',colors(2, :), 'linewidth',LineWidth,'Marker',markers(2),'MarkerSize',MarkerSize},1);
    hold on 
    plot3 =  shadedErrorBar(con_u/10,coverall(:,2),coverstdall(:,2),{'color',colors(3, :),'linewidth',LineWidth,'Marker',markers(3),'MarkerSize',MarkerSize},1);
    hold on
    plot4 =  shadedErrorBar(con_u/10,coverall(:,3),coverstdall(:,3),{'color',colors(4, :),'linewidth',LineWidth,'Marker',markers(4),'MarkerSize',MarkerSize},1);
    hold on
    line(xlim,[0.9 0.9], 'Color', 'black', 'LineWidth', 2);

    hold off
     xticklabels([]); 
    title(sprintf("$n = %d$", n),'interpreter','latex');
    ylim([0.4,1]);
    if inx==1  && idgra==1
        ylabel('Coverage prob.','fontsize',font_size+2);  
    else
        if inx == length(Nall) && idgra==length(graphonall)
            yticklabels([]); 
            ylabel('');
        else
        yticklabels([]); 
        ylabel('');
        end
    end
    text(0.01,0.4+0.05,graphonoutput,'fontsize',font_size+2,'interpreter','latex')
    hold on 
    rectangle('position',[0 0.4 widthbox 0.1],'LineWidth', 1.2)
    hold off 
    set(gca,'fontsize',font_size+4);
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.3;
    pbaspect([1.1 1 1])

    ax = nexttile((inx-1)*length(graphonall)+idgra+6*1);
    box on;
    plot1 = shadedErrorBar(con_u/10,CIlenall(:,6),CIlenstdall(:,6),{'color',colors(1, :), 'linewidth',LineWidth,'Marker',markers(1),'MarkerSize',MarkerSize},1);
    hold on 
    plot2 = shadedErrorBar(con_u/10,CIlenall(:,1),CIlenstdall(:,1),{'color',colors(2, :), 'linewidth',LineWidth,'Marker',markers(2),'MarkerSize',MarkerSize},1);
    hold on 
    plot3 =  shadedErrorBar(con_u/10,CIlenall(:,2),CIlenstdall(:,2),{'color',colors(3, :),'linewidth',LineWidth,'Marker',markers(3),'MarkerSize',MarkerSize},1);
    hold on 
    plot4 = shadedErrorBar(con_u/10,CIlenall(:,3),CIlenstdall(:,3),{'color',colors(4, :),'linewidth',LineWidth,'Marker',markers(4),'MarkerSize',MarkerSize},1);
    hold on
   line(xlim,[cutoffCI cutoffCI],'LineStyle','--', 'Color', 'black', 'LineWidth', 2);
    hold off
    xticks([0:0.1:1]);  xticklabels({'0','','','','','0.5','','','',''}); 
    xtickangle(0);
    xlabel('$\xi_{n+1}$', 'interpreter','latex', 'FontSize', font_size);
    ylim([0,6])
    if inx==1 && idgra ==1
        ylabel('CI length','fontsize',font_size+2);
    else
        yticklabels([]); 
        ylabel('');
    end
    if inx == length(Nall) && idgra == length(graphonall)
        hL = legend([plot1.mainLine,plot2.mainLine,plot3.mainLine,plot4.mainLine], Methodname,'NumColumns',2,'location','north','fontsize', font_size);
        set(hL, 'position',[0.558 0.61 0.2 0.12])
    end
   set(gca,'fontsize',font_size+4);
   ax = gca;
   ax.TitleFontSizeMultiplier = 1.5;
   pbaspect([1.1 1 1])

end

end
saveas(fig,sprintf("./figure/simu1_total_allgraphon_supp.png"))





