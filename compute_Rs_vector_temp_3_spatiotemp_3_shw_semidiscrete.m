function [L2Rt,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs_vector_temp_3_spatiotemp_3_shw_semidiscrete(x,dist_x_pl,dist_x_min,h_old,h_new,evalt,tj,dt,f_h_old,f_h_new,spatial_disc,flux_fn)
global nq xq wq
g=9.81;
dx= x(2)-x(1);
L2Rt = 0;
L2Rt_arr = zeros(size(x));
c_0_coeff_arr_new= zeros(size(x));
c_0_coeff_arr_old= zeros(size(x));

c_0_ts = h_old;
c_1_ts = f_h_old;
c_2_ts = -(1/dx)*(f_h_new-f_h_old) + (3/(dx^2))*(h_new - h_old - dx*f_h_old);
c_3_ts =  (1/(dx^2))*(f_h_new-f_h_old) - (2/(dx^3))*(h_new - h_old - dx*f_h_old);

% c_0_ts = h_new;
% c_1_ts = 1/(2*dx)*[0,h_new(1,3:end)-hnew(1,1:end-2),0;0,h_new(2,3:end)-hnew(1,2:end-2),0]%flux_fn(dist_x_pl,dt, h_new);
% c_2_ts = -(1/(dx)) * (circshift(c_1_ts,-1,2)-c_1_ts) + (3/(dx^2))*(circshift(h_new,-1,2) - h_new - dx*c_1_ts);
% c_3_ts = (1/(dx^2)) * (circshift(c_1_ts,-1,2)-c_1_ts) - (2/(dx^3))*(circshift(h_new,-1,2) - h_new - dx*c_1_ts);
% 
% 

%Now we separate the two reconstructions
c_0_ts_h = c_0_ts(1,:);
c_0_ts_hv = c_0_ts(2,:);

c_1_ts_h = c_1_ts(1,:);
c_1_ts_hv = c_1_ts(2,:);

c_2_ts_h = c_2_ts(1,:);
c_2_ts_hv = c_2_ts(2,:);

c_3_ts_h = c_3_ts(1,:);
c_3_ts_hv = c_3_ts(2,:);




for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU_h = c_0_ts_h + c_1_ts_h .* diff_x + c_2_ts_h .* diff_x .^2 +c_3_ts_h.*diff_x.^3;
    IU_hv = c_0_ts_hv + c_1_ts_hv .* diff_x + c_2_ts_hv .* diff_x .^2 +c_3_ts_hv.*diff_x.^3;

    IUx_h = ( c_1_ts_h + 2* c_2_ts_h .*diff_x +  3* c_2_ts_h .*diff_x.^2);
    IUx_hv = ( c_1_ts_hv + 2* c_2_ts_hv .*diff_x +  3* c_2_ts_hv .*diff_x.^2);
    
    IUt_h = -f_h_new(1,:);
    IUt_hv = -f_h_new(2,:);
    
    R_h = IUt_h +IUx_hv;
    R_hv = IUt_hv + ((IU_h).^(-2)).*(2*IU_hv.*IUx_hv.*IU_h - (IU_hv.^2).*IUx_h) + .5*g*IU_h.^2;
    
    
    integral_txq = (wq(iq)*dist_x_pl.*(R_h).^2) + (wq(iq)*dist_x_pl.*(R_hv).^2);
    
   
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_arr = L2Rt_arr + (integral_txq);

    
    
end


end