function [basis_fun] = local_basis_fun_1D(x, vertices, basis_type, basis_index, basis_der_x)
%%
% 100: 1D P0, piecewise constant
% 101:1D linear
% 102:1D quadratic
%%
if vertices(2) < vertices(1)
    temp = vertices(1);
    vertices(1) = vertices(2);
    vertices(2) = temp;
end
h = vertices(2) - vertices(1); % 判断哪一个大


if basis_type == 100
    basis_fun = 1;
elseif basis_type == 101
    if basis_index == 1 % the first basis function
        if basis_der_x == 0
            basis_fun = (vertices(2)-x)/h;
        elseif basis_der_x == 1
            basis_fun = -1/h;
        elseif basis_der_x >= 2 % 判断整数
            basis_fun = 0;
        else
            error('wrong input for basis derivative order!');
        end
    elseif basis_index == 2
        if basis_der_x == 0
            basis_fun = (x-vertices(1))/h;
        elseif basis_der_x == 1
            basis_fun = 1/h;
        elseif basis_der_x >= 2 % 判断整数
            basis_fun = 0;
        else
            error('wrong input for basis derivative order!');
        end
    else
        error('wrong input for basis index!');
    end
elseif basis_type == 102
    error('wrong input for basis type!');
end
end

