

clear
rng(1)
missingcase = 2;
sparsity_parameters_a = 1;
scaledata = 5;
grid_factor=1.2;
scalen = 1000;
grid_num = 2*scalen;
typemissing = "maxmissing";
boot = 1000;
Nall = [50,100, 200,400];
graphonnameall = ["f1","f2","f3"];
misperall = 0.05:0.05:0.20;

for m = Nall
for GraphonName1 = graphonnameall
for misper = misperall

switch GraphonName1
    case "f1"
        fixed_u = 0.9;
	case "f2"
        fixed_u = 0.7;
    case "f3"
        fixed_u = 0.6;
end

CIlen = [];
cover = [];
timerecord = [];
q = 1/m^(1/2);
con_row = m;
misn = ceil(m*misper);

for inx = 1:boot
    X = rand(m,1);
    X(con_row,1) = fixed_u;
    W = graphon(X,X,sparsity_parameters_a,GraphonName1)*scaledata;  
    W = generate_randW(W,scaledata);
    C_0_setting = (max(max(abs(W))))*grid_factor;
    Wcomp = W;
    tri_vals = Wcomp(triu(true(m), 1));
    switch typemissing 
        case "maxmissing"
            sortlist = sort(tri_vals, 'descend');
            cutoff = sortlist(misn);
            Wcopy = Wcomp + diag(ones(m,1)*-10*C_0_setting);
            Wcopy(Wcopy>=cutoff) = NaN; 
        case "minmissing"
            sortlist = sort(tri_vals, 'ascend'); 
            cutoff = sortlist(misn);
            Wcopy = Wcomp + diag(ones(m,1)*10*C_0_setting);
            Wcopy(Wcopy<=cutoff) = NaN; 
    end
    Wcopy = Wcopy-diag(diag(Wcopy));
    selected_row = m;
    selected_col = m-1;
    pred = Wcomp(selected_row,selected_col);
    Wdata = Wcomp;
    switch missingcase
        case 0
        case 1
            Wdata(1,1:misn) = NaN;
            Wdata(1:misn,1) = NaN;
        case 2
             Wdata = Wcopy;
    end
    Wdata(selected_row,selected_col) = NaN;
    Wdata(selected_col,selected_row) = NaN;
    tic;
    midpoint = (max(max(Wdata))+min(min(Wdata)))/2;
    C_0 = (max(max(Wdata))- midpoint)*grid_factor;
    Wdata = Wdata - midpoint;
    pred = pred - midpoint;
    boundx = 1;
    boundy = 2*C_0;  
    mask_cols = [1:selected_row-1,selected_row+1:m];
    mask_rows = [1:selected_row-1,selected_row+1:m];
    maxy = C_0;
    miny = -C_0;
    grid_int = (maxy-miny)/(grid_num);
    yseq = [miny:  grid_int: maxy];
    WdataY = Wdata(selected_row, mask_cols);
    WdataX = Wdata(mask_rows, mask_cols);
    MY = isnan(WdataY);
    MX = isnan(WdataX);
    WdataX(MX) = C_0;
    WdataY(MY) = C_0;
    [D,hlist] = NeighborhoodSmoothing(WdataX,q);
    DY_temp = kernel_Y(WdataY);
    DY_temp(MY,:) = ones(sum(MY==1),1)*max(C_0-WdataY,WdataY+C_0);
    DY_temp(:,MY) = max(C_0-WdataY,WdataY+C_0)'*ones(sum(MY==1),1)';
    DY_temp = DY_temp - diag(diag(DY_temp));
    D_guessed = D;
    D_guesnozero = D_guessed ~= 0;
    tau = zeros(m-1,1);
    recordboundx = zeros(m-1);
    for ink = 1:m-1
        if ink == selected_col
            sumlist = zeros(m-1,1);
            seqlist = 1:(m-1);
            nozeroinx = seqlist(D_guesnozero(ink,:));
            for jnkindex = 1:length(nozeroinx)
                jnk = nozeroinx(jnkindex);
                if jnk == selected_col
                    sumlist(selected_col,1) = 0;
                    kernelboundX = 0;
                else
                firstterm1 = 2*C_0*(MX(ink,:))*abs(WdataX);
                firstterm2 = 2*C_0*(MX(jnk,:))*abs(WdataX);
                firstterm = sum(firstterm1) + sum(firstterm2) - firstterm1(ink) - firstterm1(jnk) -firstterm2(ink) -firstterm2(jnk);
                secondterm1 = 2*C_0*MX*abs(WdataX(:,ink)-WdataX(:,jnk));
                secondterm = sum(secondterm1) - secondterm1(ink) - secondterm1(jnk);
                kernelboundX = min(1/hlist(ink)/(m-1)/(m-3)*(firstterm+secondterm),boundx);
                sumlist(jnk,1) = D(ink,jnk)*min((MY(jnk))*2*C_0,2*C_0)+...
                    DY_temp(ink,jnk)*kernelboundX;
            
                end
                 recordboundx(ink,jnk) = kernelboundX;
            end
            tau(selected_col,1) = sum(sumlist);
        else
            sumlist = zeros(m-1,1);
            seqlist = 1:(m-1);
            nozeroinx = seqlist(D_guesnozero(ink,:));
            for jnkindex = 1:length(nozeroinx)
                
                jnk = nozeroinx(jnkindex);
                kernelboundX = 0;
                if jnk ~= ink && jnk ~= selected_col
                    firstterm1 = 2*C_0*(MX(ink,:))*abs(WdataX);
                    firstterm2 = 2*C_0*(MX(jnk,:))*abs(WdataX);
                    firstterm = sum(firstterm1) + sum(firstterm2) - firstterm1(ink) - firstterm1(jnk) -firstterm2(ink) -firstterm2(jnk);
                    secondterm1 = 2*C_0*MX*abs(WdataX(:,ink)-WdataX(:,jnk));
                    secondterm = sum(secondterm1) - secondterm1(ink) - secondterm1(jnk);
                    kernelboundX = min(1/hlist(ink)/(m-1)/(m-3)*(firstterm+secondterm),boundx);
                    sumlist(jnk,1) = D(ink,jnk)*min((MY(ink)+MY(jnk))*2*C_0,2*C_0)+...
                        DY_temp(ink,jnk)*kernelboundX;
                elseif jnk == selected_col
                    firstterm1 = 2*C_0*(MX(ink,:))*abs(WdataX);
                    firstterm2 = 2*C_0*(MX(jnk,:))*abs(WdataX);
                    firstterm = sum(firstterm1) + sum(firstterm2) - firstterm1(ink) - firstterm1(jnk) -firstterm2(ink) -firstterm2(jnk);
                    secondterm1 = 2*C_0*MX*abs(WdataX(:,ink)-WdataX(:,jnk));
                    secondterm = sum(secondterm1) - secondterm1(ink) - secondterm1(jnk);
                    kernelboundX = min(1/hlist(ink)/(m-1)/(m-3)*(firstterm+secondterm),boundx);
                    sumlist(selected_col,1) = D(ink,selected_col)*min((MY(ink))*2*C_0,2*C_0)+...
                     DY_temp(ink,jnk)*kernelboundX;

                end
                recordboundx(ink,jnk) = kernelboundX;
            end
          
            tau(ink,1) = sum(sumlist);
        end
    end
    probupper = [];
    CIset = [];
    D_pred = [];
    record = [];
    ciind =0;  
    DY = kernel_Y(WdataY);
    for linx = 1:length(yseq)
        y = yseq(linx);
        WdataY(selected_col) =y;
        Dpred = abs(y-WdataY);
        DY(selected_col,:) = Dpred;
        DY(:,selected_col) = Dpred;
        D_total_new = sum(D.*DY,2)';
        Un1 = D_total_new(selected_col)-tau(selected_col);
        L = D_total_new'+tau;
        probupper(linx) = mean(Un1<L);
        if probupper(linx)>=0.1
           ciind = ciind +1;
           CIset(ciind) = y;
        end
    end
    if isempty(CIset)
        CIlen(inx) = NaN;
        cover(inx) = 0;
    else
        cover(inx) = any(ismember(round(CIset/C_0, 3), round(pred/C_0, 3)));
        posiind = CIset>0;
        negiind = ~posiind;
        posiset = CIset(posiind);
        negiset = CIset(negiind);
        CIlen(inx) = length(posiset)/scalen*C_0+length(negiset)/scalen*C_0;
    end
    timerecord(inx) = toc;

end


save(sprintf("./result/conf_net_stable_two_%d_sit_%d_%d_%s_%d_rand_%s_%d",m, floor(fixed_u*10), missingcase,GraphonName1,floor(misper*100),typemissing,floor(scaledata)),"CIlen","cover","timerecord")
end
end
end
