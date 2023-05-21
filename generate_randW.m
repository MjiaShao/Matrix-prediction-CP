
function W = generate_randW(W,scaledata)
	
	n = size(W, 1);
    er = rand(n,n)*0.04*scaledata-0.02*scaledata;
    er = triu(er);
    W = triu(W);
    W = W+er;
    W = W +triu(W,1)';
    W = W - diag(diag(W));
	
end