function [Gauss_weights_local_1D,Gauss_nodes_local_1D] = generate_Gauss_local_1D(Gauss_weights_reference_1D, Gauss_nodes_reference_1D, lower_bound, upper_bound)
%% Generate the Gauss coefficients and Gauss points on an arbitrary interval [lower_bound,upper_bound] by using affine tranformation.
% Gauss_coefficient_local_1D,Gauss_point_local_1D:the Gauss coefficients and Gauss points on an arbitrary interval.
% Gauss_coefficient_reference,Gauss_point_reference: the Gauss coefficients and Gauss points on the reference interval [-1,1].
%%

Gauss_weights_local_1D = (upper_bound-lower_bound)*Gauss_weights_reference_1D/2;
Gauss_nodes_local_1D = (upper_bound-lower_bound)*Gauss_nodes_reference_1D/2 + (upper_bound+lower_bound)/2;

end