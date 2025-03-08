function [dr_LVLH, dv_LVLH] = rotate_ECI2LVLH(r_B, v_B, r_A, v_A)
    
    L_x = r_A/(norm(r_A));
    L_z = (cross(r_A,v_A)/ norm(cross(r_A,v_A)));
    L_y = cross(L_z,L_x);

    Om_A = cross(r_A,v_A)/(norm(r_A)^2);
    R_i = [L_x L_y L_z]';
    dr_i = r_B - r_A;
    dv_i = v_B - v_A;

    
    dr_LVLH = R_i * dr_i;
    dv_LVLH = (dv_i - (cross(Om_A, dr_i)));
end