clc;
clear;
close all;

global nq xq wq
nq = 2;
gauss();
T=10;
u_init=1;
v_init =1;
which_rec = "normal_one"; 

% for i =1:length(u_init) 
% 
%         var_out = fn_lxf_shw();
% 
% 
% 
% end

xlim_p = T;
ylim_p1 = 1e-10;
ylim_p2 = 5e-7;
EI_plot_y_lim = 10;
R_plot_lims = [1e-11, 2e-7];
if which_rec == "normal_one"
    for i =1:length(u_init)      

        var_out = fn_plot_shw_dambreak(1,1,ylim_p1, ylim_p2, xlim_p , EI_plot_y_lim, R_plot_lims);

    end
elseif which_rec == "reverse"
    for i =1:length(u_init)       
        var_out = fn_plot_shw_dambreak(1,1,ylim_p1, ylim_p2, xlim_p , EI_plot_y_lim, R_plot_lims);
    end
end

function [] = gauss
% For n point gauss quadrature, return evaluation points and weights for
% gauss quadrature over [-1,1]
global nq xq wq
if nq == 1
    xq = 0;
    wq = 1;
elseif nq == 2
    xq = [-0.57735026918962576451, 0.57735026918962576451];
    wq = [0.5, 0.5];
elseif nq == 3
    xq = [-0.7745966692414834, 0, 0.7745966692414834];
    wq = 0.5*[0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
elseif nq > 3
    fprintf(1,'No Gauss quadrature implemented of this degree\n')
    return
end
end