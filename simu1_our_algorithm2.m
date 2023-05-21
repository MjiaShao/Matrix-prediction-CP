clear
rng(1)
missingcase = 0;
sparsity_parameters_a = 1;
boot = 1000;
grid_factor=1;
scaledata = 5;
grid_num = 2*1000;

Nall = [50,100, 200,400];
xi_n = 0.1:0.1:0.9;
graphonnameall = ["f1","f2","f3"];

for m = Nall
for GraphonName1 = graphonnameall
for fixed_u = xi_n

CIlen = [];
cover = [];
timerecord = [];
q = 1/sqrt(m);
C_0 = 1*scaledata;
con_row = m;

for inx = 1:boot
    X = rand(m,1);
    X(con_row,1) = fixed_u;
    W = graphon(X,X,sparsity_parameters_a,GraphonName1)*scaledata;  
    W = generate_randW(W,scaledata);
    Wcomp = W;
    selected_row = m;
    selected_col = m-1;
    pred = Wcomp(selected_row,selected_col);
    Wdata = Wcomp;
    Wdata(selected_row,selected_col) = NaN;
    Wdata(selected_col,selected_row) = NaN;
    tic;
    mask_cols = [1:selected_row-1,selected_row+1:m];
    mask_rows = [1:selected_row-1,selected_row+1:m];
    maxy = C_0;
    miny = -C_0;
    grid_int = (maxy-miny)/(grid_num);
    yseq = [miny:  grid_int: maxy];
    WdataY = Wdata(selected_row, mask_cols);
    WdataX = Wdata(mask_rows, mask_cols);
    D = NeighborhoodSmoothing(WdataX,q);
    tau = zeros(m-1,1);
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
        cover(inx) = any(ismember(round(CIset/scaledata, 3), round(pred/scaledata, 3)));
        posiind = CIset>0;
        negiind = ~posiind;
        posiset = CIset(posiind);
        negiset = CIset(negiind);
        CIlen(inx) = length(posiset)/1000*scaledata+length(negiset)/1000*scaledata;
    end
    timerecord(inx) = toc;

end

save(sprintf("./result/conf_net_stable_two_%d_sit_%d_%d_%s_%d_con",m, floor(fixed_u*10), missingcase,GraphonName1,floor(scaledata)),"CIlen","cover","timerecord")
end
end
end
