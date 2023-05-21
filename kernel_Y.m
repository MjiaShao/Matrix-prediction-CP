function [DY] = kernel_Y(Ydata)
	N = length(Ydata);
	DY = zeros(N, N);
	for i = 1:N
		for j = i:N
			if(j>=i)
				DY(i, j) =abs(Ydata(i)-Ydata(j)); 
				DY(j, i) = DY(i, j);
			end
		end
	end
	
	
end