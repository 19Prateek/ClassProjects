function [Z,nbr_roots]=cubic_roots(coeff)
% cubic_roots returns the sorted real roots of the cubic equation with coefficents 
% given by the argument. The input must be a 1x4 vector.
total_roots = roots(coeff);
if (isreal(total_roots) == 0)
    nbr_roots=1;
    for i=1:3
        if (isreal(total_roots(i)) ==1)
            Z = total_roots(i);            
        end
    end
else
    nbr_roots = 3; % all roots are real
    Z_temp = sort(total_roots);
    if (Z_temp(1) < 0)
        Z = Z_temp(3);% if the lowest root is less than zero 
        nbr_roots = 1;% there is only one positive root
    else
        Z=Z_temp;
    end
end
