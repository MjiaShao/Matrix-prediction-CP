clear
font_size = 20+2;
MarkerSize = 15;  LineWidth = 2;
con_u= [1:1:9];
Nall = [50,100,200, 400];
graphonall = ["f1","f2", "f3"];
q = 5;
Methodname = ["Our Alg.1", "Our Alg.2","missMDA", "softImpute","PPCA","mice"];
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
        widthbox = 3.5;
    case "f3"
        widthbox = 3.9;
    case "f1"
        widthbox = 2.8;
end
for inx = 1:length(Nall)

    n = Nall(inx);
    timecombineall1 = [];
    timecombineall2 = [];
    timecombineall3 = [];
    timecombineall4 = [];
    timecombineall5 = [];
    timecombineall6 = [];

    
    for k = 1:length(con_u)
        u = con_u(k);
		load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_5_con.mat",n,u,graphonname))    
        timecombineall1 = [timecombineall1,timerecord];

	    load(sprintf("./result/conf_net_stable_two_%d_sit_%d_0_%s_rand_5_svd_onemissing.mat",n,u,graphonname))        
        timecombineall6 = [timecombineall6,timerecord];
        if n <= 200    
            data = readtable(sprintf("./result/conf_net_missMDA_%d_%d_%s_con_5.csv",n,u,graphonname));
            timecombineall2 = [timecombineall2,data.timerecord'];
        else
            timecombineall2  = NaN;
        end

		data = readtable(sprintf("./result/conf_net_soft_%d_%d_%s_con_5.csv",n,u,graphonname));
        timecombineall3 = [timecombineall3,data.timerecord'];


        if n<=50
            data = readtable(sprintf("./result/conf_net_ppca_%d_%d_%s_con_5.csv",n,u,graphonname)); 
            timecombineall4 = [timecombineall4,data.timerecord'];
 
            data = readtable(sprintf("./result/conf_net_mice_%d_%d_%s_con_5.csv",n,u,graphonname));
            timecombineall5 = [timecombineall5,data.timerecord'];
        else
                  
            timecombineall4 = NaN;
            timecombineall5 = NaN;
        end
    end
      fig = figure("visible","on");
    timecombinetotal = [timecombineall6,timecombineall1,timecombineall2,timecombineall3,timecombineall4,timecombineall5]';
    timecombinetotal = log(timecombinetotal);
    g =  [zeros(length(timecombineall6), 1); ones(length(timecombineall1), 1); 2*ones(length(timecombineall2), 1);...
        3*ones(length(timecombineall3), 1); 4*ones(length(timecombineall4), 1);5*ones(length(timecombineall5), 1);];
    h = boxplot(timecombinetotal, g,'PlotStyle','traditional','Colors','k','Symbol',".");
    set(h,{'linew'},{1.3})
    hold on 
    text(0.5,5.3,graphonoutput,'fontsize',font_size+8,'interpreter','latex')
    hold on 
    rectangle('position',[0.5 4.6 widthbox 1.4],'LineWidth', 1.2)
    hold off
    xticklabels(Methodname);
    a = get(gca, 'xticklabels');
    set(gca,'xticklabels',a,'fontsize',font_size+8)
    xtickangle(25)
    title(strcat("Time cost, n = ", string(n)),'interpreter','latex');
    ylim([-5,6]);
    ylabel('Log time cost','fontsize',font_size+12);  
    set(gca,'fontsize',font_size+4);
    ax = gca;
    ax.TitleFontSizeMultiplier = 1.5;
    pbaspect([1.2 1 1])
    set(gca,'units','centimeters')
    set(gcf,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gca, 'Position',[2.4 2.8 pos(3)+ti(1)+ti(3)-1.5 pos(4)+ti(2)+ti(4)-1]);
     set(gcf, 'Units', 'Inches', 'Position', [0.5 0.1 5.8 5.7], 'PaperUnits', 'Inches', 'PaperSize', [5.5, 5.5])

    saveas(fig,sprintf("./figure/simu1_time_%d_%s_boxplot.png",n,graphonname))

end

end





