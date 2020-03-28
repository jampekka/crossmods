function [pdf, noncross] = vddm(tau, distance, tau_dot, ehmi, p)
    prior_tau = distance/(50/3.6);
    mangled_tau = p.dist_coeff*(prior_tau - tau) + tau + p.dot_coeff*(tau_dot + 1) + p.ehmi_coeff*ehmi;
    [pdf, noncross] = vddm_tau(mangled_tau, p);
end
