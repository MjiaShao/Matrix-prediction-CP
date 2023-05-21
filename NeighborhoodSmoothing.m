function [Kernel_mat,hlist, D] = NeighborhoodSmoothing(Adata,q)
	
	N = size(Adata, 1);

	D = zeros(N, N);
	A_sq = Adata*Adata;
	for i = 1:N
		for j = i:N
			if(j>=i)
                distance_all = abs(A_sq(i, :) - A_sq(j, :));
                distance_all(i) = 0;
                distance_all(j) = 0;
                D(i, j) = abs(sum(distance_all)/N/(N-2));
				D(j, i) = D(i, j);
			end
		end
    end
	for i = 1:N
		Kernel_mat(i, :) = max(1-D(i,:)/quantile(D(i, :),q),0);
        hlist(i) = quantile(D(i, :),q);
    end
	
end