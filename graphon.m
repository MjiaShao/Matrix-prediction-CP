
function z = graphon(x,y,sparsity_parameters,GraphonName)
	
	x = x(:);  y = y(:);
	
	n = length(x);
    rho = sparsity_parameters;

	switch GraphonName
	case 'f1'		
		XX = x*ones(1,n);  YY = ones(n,1)*y';
        z = (XX+YY)./2 - 0.15;		
	case 'f2'
		XX = x*ones(1,n);  YY = ones(n,1)*y';
		z = 0.5*(cos(0.1./((XX-0.5).^3+(YY-0.5).^3+0.01))).*(max(XX.^(2/3),YY.^(2/3))) + 0.4;
	case 'f3'		
        XX = x*ones(1,n);  YY = ones(n,1)*y';
		z = ((XX.^2+YY.^2)/3 .* cos(1./(XX.^4+YY.^4)) + 0.15);		
	end
	
	z = z*rho;
	
end
