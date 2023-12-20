function [result] = triangular_reference_basis(x,y,basis_type,basis_index,derivative_degree_x,derivative_degree_y)
%% the reference FE basis function on triangle ABC where A=(0,0), B=(1,0) and C=(0,1).
% x,y: the coordinates of the point where we want to evaluate the reference FE basis function.
% basis_type: the type of the FE.
%basis_type=0: 2D constant FE.
% basis_type=1: 2D linear FE.
% basis_type=2: 2D Lagrange quadratic FE.
% basis_type=10: 2D Crouzeix-Raviart FE.
% basis_index: the index of FE basis function to specify which FE basis function we want to use.
% derivative_degree_x: the derivative degree of the FE basis function with respect to x.
% derivative_degree_y: the derivative degree of the FE basis function with respect to y.
%%

if basis_type==0
    
    if derivative_degree_x==0&&derivative_degree_y==0
        result = 1;
    end

elseif basis_type==1

    if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            result=1-x-y;
        elseif basis_index==2
            result=x;
        elseif basis_index==3
            result=y;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==0
        
        if basis_index==1
            result=-1;
        elseif basis_index==2
            result=1;
        elseif basis_index==3
            result=0;
        end

    elseif derivative_degree_x==0&&derivative_degree_y==1
        
        if basis_index==1
            result=-1;
        elseif basis_index==2
            result=0;
        elseif basis_index==3
            result=1;
        end
        
    end
    
elseif basis_type==2
    
    if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            result=1-3*x-3*y+2*x.^2+2*y.^2+4*x.*y;
        elseif basis_index==2
            result=2*x.^2-x;
        elseif basis_index==3
            result=2*y.^2-y;
        elseif basis_index==4
            result=4*x-4*x.^2-4*x.*y;
        elseif basis_index==5
            result=4*x.*y;
        elseif basis_index==6
            result=4*y-4*y.^2-4*x.*y;
        end
             
    elseif derivative_degree_x==1&&derivative_degree_y==0
 
        if basis_index==1
            result=-3+4*x+4*y;
        elseif basis_index==2
            result=4*x-1;
        elseif basis_index==3
            result=0;
        elseif basis_index==4
            result=4-8*x-4*y;
        elseif basis_index==5
            result=4*y;
        elseif basis_index==6
            result=-4*y;
        end           

                      
    elseif derivative_degree_x==0&&derivative_degree_y==1
            
        if basis_index==1
            result=-3+4*y+4*x;
        elseif basis_index==2
            result=0;
        elseif basis_index==3
            result=4*y-1;
        elseif basis_index==4
            result=-4*x;
        elseif basis_index==5
            result=4*x;
        elseif basis_index==6
            result=4-8*y-4*x;
        end
      
    elseif derivative_degree_x==2&&derivative_degree_y==0  
        
        if basis_index==1
            result=4;
        elseif basis_index==2
            result=4;
        elseif basis_index==3
            result=0;
        elseif basis_index==4
            result=-8;
        elseif basis_index==5
            result=0;
        elseif basis_index==6
            result=0;
        end

    elseif derivative_degree_x==0&&derivative_degree_y==2 

        if basis_index==1
            result=4;
        elseif basis_index==2
            result=0;
        elseif basis_index==3
            result=4;
        elseif basis_index==4
            result=0;
        elseif basis_index==5
            result=0;
        elseif basis_index==6
            result=-8;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==1 

        if basis_index==1
            result=4;
        elseif basis_index==2
            result=0;
        elseif basis_index==3
            result=0;
        elseif basis_index==4
            result=-4;
        elseif basis_index==5
            result=4;
        elseif basis_index==6
            result=-4;
        end 
      
    end

    
elseif basis_type==10

    if derivative_degree_x==0&&derivative_degree_y==0
        
        if basis_index==1
            result=1-2*y;
        elseif basis_index==2
            result=-1+2*x+2*y;
        elseif basis_index==3
            result=1-2*x;
        end

    elseif derivative_degree_x==1&&derivative_degree_y==0
        
        if basis_index==1
            result=0;
        elseif basis_index==2
            result=2;
        elseif basis_index==3
            result=-2;
        end

    elseif derivative_degree_x==0&&derivative_degree_y==1
        
        if basis_index==1
            result=-2;
        elseif basis_index==2
            result=2;
        elseif basis_index==3
            result=0;
        end
        
    end
       
end
end