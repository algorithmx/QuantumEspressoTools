_FDNM__ps_kp_ct_kBar_dg_(ps,kp,ct,kBar,dg) = 
    "ps___$(pseudo_mode_name(ps))___kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_cutoff_$(ct[1]),$(ct[2])_kBar_$(kBar)_dg_$(dg)"

_FDNM__ps_kp_ct_(ps,kp,ct) = 
    "ps___$(pseudo_mode_name(ps))___kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_cutoff_$(ct[1]),$(ct[2])"

_FDNM__ps_kp_kp_ct_(ps,kp,kp1,ct) = 
    "ps___$(pseudo_mode_name(ps))___kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_kp_$(kp1[1]),$(kp1[2]),$(kp1[3]),$(kp1[4]),$(kp1[5]),$(kp1[6])_cutoff_$(ct[1]),$(ct[2])"
