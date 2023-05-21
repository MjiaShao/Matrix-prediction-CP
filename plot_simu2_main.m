clear
font_size = 20;
MarkerSize = 15;  LineWidth = 2;
misperall= [0:5:20];
Nall = [50,200];
graphonall = ["f1","f2", "f3"];
colors = [	[1 0 0]; 	[0 0 1]; 	[0.3010, 0.7450, 0.9330]* 0.8;  	[0, 0.5, 0]; [0.4940, 0.1840, 0.5560];   	[0.9290, 0.6940, 0.1250]* 0.8];
markers = ['o','+','*','v','X','s'];

fig = figure('Visible','on','Position',[100,100,1150*2-50,670],'Units','centimeters');      
tlo = tiledlayout(2,6,'TileSpacing','none','Padding','compact');  
for i = 2
    for j = 6
        ax = nexttile;
    end
end
for idgra = 1:length(graphonall)
graphonname = graphonall(idgra);
switch graphonname
    case "f2"
        u = 7;
    case "f1"
        u = 9;
    case "f3"
        u = 6;
end
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
        widthbox = 0.20;
    case "f3"
        widthbox = 0.23;
    case "f1"
        widthbox = 0.17;
end
Methodname = {"Our Alg.1","Our Alg.2","missMDA", "softImpute","PPCA","mice"};
typemissing = "maxmissing";
for inx = 1:length(Nall)
    n = Nall(inx);
    coverall = [];
    CIlenall = [];
    timeall = [];
    coverstdall = [];
    CIlenstdall = [];
    timestdall = [];
    for k = 1:length(misperall)
        misper = misperall(k);
        if misper == 0
            load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_5_con.mat",n,u,graphonname))           
            coverall(k,1) = nanmean(cover);
            coverstdall(k,1) = sqrt(coverall(k,1)*(1-coverall(k,1))/length(cover));
            CIlenall(k,1) = nanmean(CIlen);
            CIlenstdall(k,1) = nanstd(CIlen)/sqrt(length(CIlen));
            timeall(k,1) = nanmean(timerecord);
            timestdall(k,1) = nanstd(timerecord)/sqrt(length(timerecord));
            load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_rand_5_svd_onemissing.mat",n,u,graphonname))        
            coverall(k,6) = nanmean(cover);        
            CIlenall(k,6) = nanmean(CIlen);
            timeall(k,6) = nanmean(timerecord);
            coverstdall(k,6) = sqrt(coverall(k,6)*(1-coverall(k,6))/length(cover));
            CIlenstdall(k,6) = nanstd(CIlen)/sqrt(length(CIlen));
             timestdall(k,6) = nanstd(timerecord)/sqrt(length(timerecord));


            data = readtable(sprintf("./result/conf_net_soft_%d_%d_%s_con_5.csv",n,u,graphonname));
            coverall(k,3) = nanmean(data.cover);
            CIlenall(k,3) = nanmean(data.CIlen);
            timeall(k,3) = nanmean(data.timerecord);
            coverstdall(k,3) = sqrt(coverall(k,3)*(1-coverall(k,3))/length(data.cover));
            CIlenstdall(k,3) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            timestdall(k,3) = nanstd(data.timerecord)/sqrt(length(data.timerecord));     
            
            if n<=50
                data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_5.csv",n,u,graphonname));
                coverall(k,2) = nanmean(data.cover);
                CIlenall(k,2) = nanmean(data.CIlen);
                timeall(k,2) = nanmean(data.timerecord);
                coverstdall(k,2) = sqrt(coverall(k,2)*(1-coverall(k,2))/length(data.cover));
                CIlenstdall(k,2) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
                timestdall(k,2) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
                data = readtable(sprintf("./result/conf_net_ppca_%d_%d_%s_con_5.csv",n,u,graphonname));
                coverall(k,4) = nanmean(data.cover);
                coverstdall(k,4) = sqrt(coverall(k,4)*(1-coverall(k,4))/length(data.cover));
                timeall(k,4) = nanmean(data.timerecord); 
                timestdall(k,4) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
                if abs(nanmean(data.CIlen))<cutoffCI
                     CIlenall(k,4)= min(abs(nanmean(data.CIlen)),cutoffCI);
                     CIlenstdall(k,4) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
                else
                    CIlenall(k,4) = min(abs(nanmean(data.CIlen)),cutoffCI);
                    CIlenstdall(k,4) = 0;
                end

                data = readtable(sprintf("./result/conf_net_mice_%d_%d_%s_con_5.csv",n,u,graphonname));
                coverall(k,5) = nanmean(data.cover);
                CIlenall(k,5) = nanmean(data.CIlen);
                timeall(k,5) = nanmean(data.timerecord);
                coverstdall(k,5) = sqrt(coverall(k,5)*(1-coverall(k,5))/length(data.cover));
                CIlenstdall(k,5) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
                timestdall(k,5) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
            elseif n == 200
                data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_5.csv",n,u,graphonname));
                coverall(k,2) = nanmean(data.cover);
                CIlenall(k,2) = nanmean(data.CIlen);
                timeall(k,2) = nanmean(data.timerecord);
                coverstdall(k,2) = sqrt(coverall(k,2)*(1-coverall(k,2))/length(data.cover));
                CIlenstdall(k,2) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
                timestdall(k,2) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
                coverall(k,4) = NaN;
                CIlenall(k,4) = NaN;
                timeall(k,4)  = NaN;        
                coverall(k,5) = NaN;
                CIlenall(k,5) = NaN;
                timeall(k,5)  = NaN;
                coverstdall(k,4) = NaN;
                CIlenstdall(k,4) = NaN;
                timestdall(k,4) = NaN;  
                coverstdall(k,5) = NaN;
                CIlenstdall(k,5) = NaN;
                timestdall(k,5) = NaN;  
            else
                coverall(k,2) = NaN;
                CIlenall(k,2) = NaN;
                timeall(k,2)  = NaN;  
                coverall(k,4) = NaN;
                CIlenall(k,4) = NaN;
                timeall(k,4)  = NaN;        
                coverall(k,5) = NaN;
                CIlenall(k,5) = NaN;
                timeall(k,5)  = NaN;
                coverstdall(k,2) = NaN;
                CIlenstdall(k,2) = NaN;
                timestdall(k,2)  = NaN; 
                coverstdall(k,4) = NaN;
                CIlenstdall(k,4) = NaN;
                timestdall(k,4) = NaN;  
                coverstdall(k,5) = NaN;
                CIlenstdall(k,5) = NaN;
                timestdall(k,5) = NaN;  
            end
        else
        load(sprintf("./result/conf_net_stable_two_%d_sit_%d_2_%s_%d_rand_%s_5.mat",n,u,graphonname,misper,typemissing))
        coverall(k,1) = nanmean(cover);
        timeall(k,1) = nanmean(timerecord);
        coverstdall(k,1) = sqrt(coverall(k,1)*(1-coverall(k,1))/length(cover));
         timestdall(k,1) = nanstd(timerecord)/sqrt(length(timerecord));
        if abs(nanmean(CIlen))<cutoffCI
             CIlenall(k,1)= min(abs(nanmean(CIlen)),cutoffCI);
             CIlenstdall(k,1) = nanstd(CIlen)/sqrt(length(CIlen));
        else
            CIlenall(k,1) = min(abs(nanmean(CIlen)),cutoffCI);
            CIlenstdall(k,1) = 0;
        end

        data = readtable(sprintf("./result/conf_net_soft_%d_%d_%s_con_%d_rand_%s_5.csv",n,u,graphonname,misper,typemissing));
        coverall(k,3) = nanmean(data.cover);
        CIlenall(k,3) = nanmean(data.CIlen);
        timeall(k,3) = nanmean(data.timerecord);
        coverstdall(k,3) = sqrt(coverall(k,3)*(1-coverall(k,3))/length(data.cover));
        CIlenstdall(k,3) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
        timestdall(k,3) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
        load(sprintf("./result/conf_net_stable_two_%d_sit_%d_2_%s_%d_rand_%s_5_svd.mat",n,u,graphonname,misper,typemissing))
        coverall(k,6) = nanmean(cover);
        CIlenall(k,6) = nanmean(CIlen);
        timeall(k,6) = nanmean(timerecord);
        coverstdall(k,6) = sqrt(coverall(k,6)*(1-coverall(k,6))/length(cover));
        CIlenstdall(k,6) = nanstd(CIlen)/sqrt(length(CIlen));
         timestdall(k,6) = nanstd(timerecord)/sqrt(length(timerecord));
        if n<=50
            data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_%d_rand_%s_5.csv",n,u,graphonname,misper,typemissing));
            coverall(k,2) = nanmean(data.cover);
            CIlenall(k,2) = nanmean(data.CIlen);
            timeall(k,2) = nanmean(data.timerecord);
            coverstdall(k,2) = sqrt(coverall(k,2)*(1-coverall(k,2))/length(data.cover));
             CIlenstdall(k,2) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            timestdall(k,2) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
            data = readtable(sprintf("./result/conf_net_ppca_%d_%d_%s_%d_con_new_rand_%s_5.csv",n,u,graphonname,misper,typemissing));
            coverall(k,4) = nanmean(data.cover);
            coverstdall(k,4) = sqrt(coverall(k,4)*(1-coverall(k,4))/length(data.cover));
            timeall(k,4) = nanmean(data.timerecord); 
            timestdall(k,4) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
            if abs(nanmean(data.CIlen))<cutoffCI
                 CIlenall(k,4)= min(abs(nanmean(data.CIlen)),cutoffCI);
                 CIlenstdall(k,4) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            else
                CIlenall(k,4) = min(abs(nanmean(data.CIlen)),cutoffCI);
                CIlenstdall(k,4) = 0;
            end
            data = readtable(sprintf("./result/conf_net_mice_%d_%d_%s_con_%d_rand_%s_5.csv",n,u,graphonname,misper,typemissing));
            coverall(k,5) = nanmean(data.cover);
            CIlenall(k,5) = nanmean(data.CIlen);
            timeall(k,5) = nanmean(data.timerecord);
            coverstdall(k,5) = sqrt(coverall(k,5)*(1-coverall(k,5))/length(data.cover));
            CIlenstdall(k,5) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            timestdall(k,5) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
        elseif n == 200
            data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_%d_rand_%s_5.csv",n,u,graphonname,misper,typemissing));
            coverall(k,2) = nanmean(data.cover);
            CIlenall(k,2) = nanmean(data.CIlen);
            timeall(k,2) = nanmean(data.timerecord);
            coverstdall(k,2) = sqrt(coverall(k,2)*(1-coverall(k,2))/length(data.cover));
             CIlenstdall(k,2) = nanstd(data.CIlen)/sqrt(length(data.CIlen));
            timestdall(k,2) = nanstd(data.timerecord)/sqrt(length(data.timerecord));  
            coverall(k,4) = NaN;
            CIlenall(k,4) = NaN;
            timeall(k,4)  = NaN;        
            coverall(k,5) = NaN;
            CIlenall(k,5) = NaN;
            timeall(k,5)  = NaN;
            coverstdall(k,4) = NaN;
            CIlenstdall(k,4) = NaN;
            timestdall(k,4) = NaN;  
            coverstdall(k,5) = NaN;
            CIlenstdall(k,5) = NaN;
            timestdall(k,5) = NaN; 
        else
            coverall(k,2) = NaN;
            CIlenall(k,2) = NaN;
            timeall(k,2)  = NaN;
            coverall(k,4) = NaN;
            CIlenall(k,4) = NaN;
            timeall(k,4)  = NaN;        
            coverall(k,5) = NaN;
            CIlenall(k,5) = NaN;
            timeall(k,5)  = NaN;
            coverstdall(k,2) = NaN;
            CIlenstdall(k,2) = NaN;
            timestdall(k,2)  = NaN; 
            coverstdall(k,4) = NaN;
            CIlenstdall(k,4) = NaN;
            timestdall(k,4) = NaN;  
            coverstdall(k,5) = NaN;
            CIlenstdall(k,5) = NaN;
            timestdall(k,5) = NaN; 
        end
        end
    end
    timeall = log(timeall);
    ax = nexttile((inx-1)*length(graphonall)+idgra);
    plot1 = shadedErrorBar(misperall/100,coverall(:,6),coverstdall(:,6),{'color',colors(1, :),'linewidth',LineWidth,'Marker',markers(1),'MarkerSize',MarkerSize},1);
    hold on 
    plot2 = shadedErrorBar(misperall/100,coverall(:,1),coverstdall(:,1),{'color',colors(2, :), 'linewidth',LineWidth,'Marker',markers(2),'MarkerSize',MarkerSize},1);
    hold on
    plot3 = shadedErrorBar(misperall/100,coverall(:,2),coverstdall(:,2),{'color',colors(3, :),'linewidth',LineWidth,'Marker',markers(3),'MarkerSize',MarkerSize},1);
    hold on 
    plot4 = shadedErrorBar(misperall/100,coverall(:,3),coverstdall(:,3),{'color',colors(4, :),'linewidth',LineWidth,'Marker',markers(4),'MarkerSize',MarkerSize},1);
    hold on
    plot5 = shadedErrorBar(misperall/100,coverall(:,4),coverstdall(:,4),{'color',colors(5, :),'linewidth',LineWidth,'Marker',markers(5),'MarkerSize',MarkerSize},1);
    hold on 
    plot6 = shadedErrorBar(misperall/100,coverall(:,5),coverstdall(:,5),{'color',colors(6, :),'linewidth',LineWidth,'Marker',markers(6),'MarkerSize',MarkerSize},1);   
    hold on
    xlim([0-0.02,0.2+0.02]);
    line(xlim,[0.9 0.9], 'Color', 'black', 'LineWidth', 2);
    hold on 
    text(0.005-0.02,0.27,graphonoutput,'fontsize',font_size+4,'interpreter','latex')
    hold on 
    rectangle('position',[0-0.02 0.2 widthbox 0.15],'LineWidth', 1.2)

    hold off 
     xticklabels([]);
    title(sprintf("$n =%d$", n),'interpreter','latex');
    ylim([0.2,1]);


    if n==50  && idgra==1
        yticks([0.2:0.2:1])
        yticklabels([0.2:0.2:1])
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'fontsize',font_size+2)
        ylabel('Cover. prob.','fontsize',font_size+8);  
    else
        yticklabels([]);
        ylabel('');
    end
    set(gca,'fontsize',font_size+4);
     ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    pbaspect([1.1 1 1])
    

    ax = nexttile((inx-1)*length(graphonall)+idgra+6*1);
    plot1 = shadedErrorBar(misperall/100,CIlenall(:,6),CIlenstdall(:,6),{'color',colors(1, :),'linewidth',LineWidth,'Marker',markers(1),'MarkerSize',MarkerSize},1);
    hold on 
    plot2 = shadedErrorBar(misperall/100,CIlenall(:,1),CIlenstdall(:,1),{'color',colors(2, :), 'linewidth',LineWidth,'Marker',markers(2),'MarkerSize',MarkerSize},1);
    hold on 
    plot3 = shadedErrorBar(misperall/100,CIlenall(:,2),CIlenstdall(:,2),{'color',colors(3, :),'linewidth',LineWidth,'Marker',markers(3),'MarkerSize',MarkerSize},1);
    hold on 
    plot4 = shadedErrorBar(misperall/100,CIlenall(:,3),CIlenstdall(:,3),{'color',colors(4, :),'linewidth',LineWidth,'Marker',markers(4),'MarkerSize',MarkerSize},1);
    hold on
    plot5 = shadedErrorBar(misperall/100,CIlenall(:,4),CIlenstdall(:,4),{'color',colors(5, :),'linewidth',LineWidth,'Marker',markers(5),'MarkerSize',MarkerSize},1);
    hold on 
    plot6 = shadedErrorBar(misperall/100,CIlenall(:,5),CIlenstdall(:,5),{'color',colors(6, :),'linewidth',LineWidth,'Marker',markers(6),'MarkerSize',MarkerSize},1);
    hold on
    xlim([0-0.02,0.2+0.02]);
    line(xlim,[cutoffCI cutoffCI],'LineStyle','--', 'Color', 'black', 'LineWidth', 2);
    hold off
	ylim([0,6]);


    if inx==1 && idgra ==1
        yticks([0:2:6])
        yticklabels([0:2:6])
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'fontsize',font_size+8)
        ylabel('CI length','fontsize',28);
    else
        yticklabels([]); 
        ylabel('');
    end

    if n == 200 && idgra == length(graphonall)
        hL = legend([plot1.mainLine,plot2.mainLine,plot3.mainLine,plot4.mainLine,plot5.mainLine,plot6.mainLine], Methodname,'NumColumns',2,'location','north','fontsize', font_size);
        set(hL, 'position',[0.558 0.61 0.2 0.12])
    end
    xticks(misperall/100);
    xticklabels(ceil(misperall/100*n))
    xtickangle(0);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',font_size+2)
    xlim([0-0.02,0.2+0.02]);
    xlabel('$\#$ miss. $A_{i,j}$''s', 'interpreter','latex', 'FontSize', font_size+10); 
    set(gca,'fontsize',font_size+4);
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    pbaspect([1.1 1 1])

end
end


saveas(fig,sprintf("./figure/simu1_prop_total_allgraphon.png"))







