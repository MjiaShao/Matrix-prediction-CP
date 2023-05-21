load('./data/Final_data.mat');
datatypeall = ["SZ", "NC"];
nummissenall = [300,1000];
seedall = 1:1:8;

for datatype = datatypeall
    for seed = seedall
        for nummissen = nummissenall
            
switch datatype 
case "SZ"
    data = cor_SZ_fisherZ_correct;
case "NC"
    data = cor_NC_fisherZ_correct;
end

rng(seed)
boot = size(data,3);
m = size(data,1);
missingcase = 2;
CIlenmean = [];
covermean = [];
timemean = [];
grid_num = 2*1000;
switch nummissen 
    case 300
        powerinq = 3;
    case 1000
        powerinq = 3;
end

q = 1/m^(1/powerinq);
grid_factor=2.5;
nsamplenan = 6;


for rowinx = 1:m
    colindices = 1:m;
    colindices(rowinx) = [];
    CIlen = [];
    cover = [];
    timerecord = [];
    
    for inx = 1:boot
        W = data(:,:,inx);
        C_0_setting = (max(max(abs(W))))*grid_factor;
        Wcomp = W;
        tri_vals = Wcomp(triu(true(m), 1));
        switch typemissing 
            case"maxmissing"
                sortlist = sort(tri_vals, 'descend'); % Sort the entries in the upper triangle in descending order
                cutoff = sortlist(nummissen);
                Wcopy = Wcomp + diag(ones(m,1)*-10*C_0_setting);
                Wcopy(Wcopy>=cutoff) = NaN; % Set the selected entries to zero
            case "minmissing"
                sortlist = sort(tri_vals, 'ascend'); % Sort the entries in the upper triangle in descending order
                cutoff = sortlist(nummissen);
                Wcopy = Wcomp + diag(ones(m,1)*10*C_0_setting);
                Wcopy(Wcopy<=cutoff) = NaN; % Set the selected entries to zero
        end
        Wcopy = Wcopy-diag(diag(Wcopy));

        selected_row = rowinx;
        selected_col = datasample(colindices,1);
        [selected_row,selected_col]=deal(max(selected_row, selected_col), min(selected_row, selected_col));
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
        mask_cols = [1:selected_row-1,selected_row+1:m];
        mask_rows = [1:selected_row-1,selected_row+1:m];
        maxy = C_0;
        miny = -C_0;
        grid_int = (maxy-miny)/(grid_num);
        yseq = [miny:  grid_int: maxy];
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
                Dpred = abs(y-WdataY);
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
            CIlen(inx) = length(posiset)/1000*C_0+length(negiset)/1000*C_0-0.001*C_0;
        end
        timerecord(inx) = toc;
    
    end
    coverall(rowinx, :) = cover(:)';
    CIlenall(rowinx, :) = CIlen(:)';
    timeall(rowinx,:) = timerecord(:)';
    covermean(rowinx,1) = nanmean(cover);
    CIlenmean(rowinx,1) = nanmean(CIlen);
    timemean(rowinx,1) = nanmean(timerecord);
end


save(sprintf("./result/data_ourmethod_%s_random_mislarge_%d_%s_%d_newmethod_%d_svd",datatype,nummissen,typemissing,floor(seed),floor(powerinq)),"CIlenall","coverall","timeall","covermean","CIlenmean","timemean")

    end
    end
end
