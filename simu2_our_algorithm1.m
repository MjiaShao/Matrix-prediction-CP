
clear
rng(1)
missingcase = 2;
sparsity_parameters_a = 1;
boot = 1000;
scaledata = 5;
grid_factor=1.2;
scalen = 1;
grid_num = 2*scalen;
typemissing = "maxmissing";
nsamplenan = 10;

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
q = 1/m^(1/3);
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
    %Wdata = Wdata - diag(diag(Wdata));
    pred = pred - midpoint;
    mask_cols = [1:selected_row-1,selected_row+1:m];
    mask_rows = [1:selected_row-1,selected_row+1:m];
    maxy = C_0;
    miny = -C_0;
    grid_int = (maxy-miny)/(grid_num);
    yseq = [miny:  grid_int: maxy];

    tau = zeros(m-1,1);
    CIsetall = [];
    missingnumber = sum(sum(isnan(Wdata)==1));
    Wdatarecord = Wdata;
    for inputnumber = 1: nsamplenan 
        switch inputnumber 
            case 1
                Wdata = Wdatarecord;
                Wdata(isnan(Wdata))  = C_0;
            case 2 
                Wdata = Wdatarecord;
                Wdata(isnan(Wdata)) = -C_0;
            case 3
                Wdata = Wdatarecord;
                allnonmissing = [ones(1,missingnumber)*C_0,ones(1,missingnumber)*-C_0];
                Wdata(isnan(Wdata)) = datasample(allnonmissing,missingnumber, replace =true);
            case 4
                Wdata = Wdatarecord;
                allnonmissing = Wdata(~isnan(Wdata));
                allnonmissing = [allnonmissing', ones(1,missingnumber)*C_0];
                Wdata(isnan(Wdata)) = datasample(allnonmissing,missingnumber, replace =true);
            case 5
                Wdata = Wdatarecord;
                allnonmissing = Wdata(~isnan(Wdata));
                allnonmissing = [allnonmissing', ones(1,missingnumber)*-C_0];
                Wdata(isnan(Wdata)) = datasample(allnonmissing,missingnumber, replace =true);
            otherwise
                Wdata = Wdatarecord;
                allnonmissing = Wdata(~isnan(Wdata));
                Wdata(isnan(Wdata)) = datasample(allnonmissing,missingnumber, replace = true);
        end
        Wdata = triu(Wdata);
        Wdata = Wdata+triu(Wdata,1)';

        WdataY = Wdata(selected_row, mask_cols);
        WdataX = Wdata(mask_rows, mask_cols);
        [u_0,s_0,v_0] = svds(WdataX,ceil(1/q));
        
        CIsetindi = [];
        ciind =0;  

        for linx = 1:length(yseq)
            y = yseq(linx);
            WdataY(selected_col) =y;
            v_0 = inv( s_0.*s_0 ) * s_0 * u_0' *  WdataY';
            WdataYhat = u_0 * s_0 * v_0;
            D_total_new = abs(WdataY'- WdataYhat);
            Un1 = D_total_new(selected_col);
            L = D_total_new';
            probupper(linx) = mean(Un1<L);
            if probupper(linx) >=0.1
               ciind = ciind +1;
               CIsetindi(ciind) = y;
            end
        end
        CIsetall = [CIsetall,CIsetindi];
    end
    CIset = unique(CIsetall);
    
    if isempty(CIset)
        CIlen(inx) = NaN;
        cover(inx) = 0;
    else
        cover(inx) = any(ismember(round(CIset/C_0, 3), round(pred/C_0, 3)));
        posiind = CIset>0;
        negiind = ~posiind;
        posiset = CIset(posiind);
        negiset = CIset(negiind);
        CIlen(inx) = length(posiset)/scalen*C_0+length(negiset)/scalen*C_0 -1/scalen *C_0;
    end


    timerecord(inx) = toc;

end


save(sprintf("./result/conf_net_stable_two_%d_sit_%d_%d_%s_%d_rand_%s_%d_svd",m, floor(fixed_u*10), missingcase,GraphonName1,floor(misper*100),typemissing,floor(scaledata)),"CIlen","cover","timerecord")
end
end
end
