function [result] = Gauss_quadrature_for_line_integral_2D(coefficient_function_1, coefficient_function_2, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, end_point_1, end_point_2, normal_vector)
%% Use Gauss quadrature to compute a line integral on an edge for a matrix.
%We will use "FE" to replace "finite element" in the comments.
%The integrand of the volume integral must be in the following format:
%a coefficient function * a trial FE basis function(or its derivatives) * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%Gauss_coefficient_reference_1D,Gauss_point_reference_1D:the Gauss coefficients and Gauss points on the reference interval [-1,1]
%end_point_1,end_point_2:the coordinates of the end points of the local interval on which we are computing the line integral
%trial_vertices: the coordinates of all vertices of the triangular element for trial functions.
%trial_basis_type:the type of the trial FE basis function.
%trial_basis_type=1:2D linear FE.  
%trial_basis_type=2:2D Lagrange quadratic FE.
%trial_basis_index: the index of trial FE basis function to specify which trial FE basis function we want to use.
%trial_derivative_degree_x:the derivative degree of the trial FE basis function with respect to x.
%trial_derivative_degree_y:the derivative degree of the trial FE basis function with respect to y.
%test_vertices: the coordinates of all vertices of the triangular element for test functions.
%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=2:2D Lagrange quadratic FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%Gauss_coefficient_local_1D,Gauss_point_local_1D:the Gauss coefficients and Gauss points on the local interval.
%% 
Gpn = length(Gauss_weights_reference_1D); % Gpn: the number of the Gauss points of the Gauss formula we are using.
result = 0;

if end_point_1(2) == end_point_2(2) % The line is horizontal.

    lower_bound = min(end_point_1(1),end_point_2(1));
    upper_bound = max(end_point_1(1),end_point_2(1));
    [Gauss_weights_local_1D, Gauss_nodes_local_1D] = generate_Gauss_local_1D(Gauss_weights_reference_1D, Gauss_nodes_reference_1D, lower_bound, upper_bound);
    for i=1:Gpn
         result = result + Gauss_weights_local_1D(i) * ( feval(coefficient_function_1, Gauss_nodes_local_1D(i), end_point_1(2)) * normal_vector(1) + feval(coefficient_function_2, Gauss_nodes_local_1D(i), end_point_1(2)) * normal_vector(2));
    end    

elseif end_point_1(1)==end_point_2(1) % The line is vertical.

    lower_bound = min(end_point_1(2),end_point_2(2));
    upper_bound = max(end_point_1(2),end_point_2(2));
    [Gauss_weights_local_1D, Gauss_nodes_local_1D] = generate_Gauss_local_1D(Gauss_weights_reference_1D, Gauss_nodes_reference_1D, lower_bound, upper_bound);
    for i=1:Gpn
         result = result + Gauss_weights_local_1D(i) * ( feval(coefficient_function_1, end_point_1(1), Gauss_nodes_local_1D(i)) * normal_vector(1) + feval(coefficient_function_2, end_point_1(1), Gauss_nodes_local_1D(i)) * normal_vector(2));
    end
    
else % The slope of the edge is in (0,infinity).[δŪ��]

    lower_bound = min(end_point_1(1),end_point_2(1));
    upper_bound = max(end_point_1(1),end_point_2(1));
    [Gauss_weights_local_1D, Gauss_nodes_local_1D] = generate_Gauss_local_1D(Gauss_weights_reference_1D, Gauss_nodes_reference_1D, lower_bound, upper_bound);
    slope = (end_point_2(2)-end_point_1(2))/(end_point_2(1)-end_point_1(1));
    Jacobi = sqrt(1+slope^2);
    for i=1:Gpn
         x = Gauss_nodes_local_1D(i);
         y = slope*(x-end_point_1(1))+end_point_1(2);
         result = result + Gauss_weights_local_1D(i) * Jacobi * feval(coefficient_function_1, x, y) ...
             * dot(triangular_RT0_local_basis(x, y, trial_vertices, trial_basis_index, trial_derivative_degree) , normal_vector) ...
             * dot(triangular_RT0_local_basis(x, y, test_vertices, test_basis_index, test_derivative_degree) , normal_vector);
    end

end