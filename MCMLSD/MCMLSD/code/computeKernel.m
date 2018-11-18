function [kern,ind]=computeKernel(e,r_max,SIGMA_X, SIGMA_TH, r_res, th_res)
%computing the kernel for the probabilistic hough tranform
        r_th = e;
        C =[SIGMA_X^2 + (r_th^2) * (SIGMA_TH^2), r_th * (SIGMA_TH^2); r_th * (SIGMA_TH^2), SIGMA_TH^2];

        C_inv = inv(C);
        C_det = det(C);
        N = 1 / (2*pi * sqrt(C_det));
        R = -round2frac(3 *sqrt(C(1,1)), r_res):r_res:round2frac(3 *sqrt(C(1,1)), r_res);
        TH = -round2frac(3 * SIGMA_TH, th_res):th_res:round2frac(3 * SIGMA_TH, th_res);
        c=C_inv(1,1);
        d=C_inv(1,2);
        g=C_inv(2,1);
        f=C_inv(2,2);
       [a,b]=meshgrid(R,TH);
        vals=N * exp(-0.5*((a*c+b*g).*a+(a*d+b*f).*b));
        [r, c] = find(vals > 0.05*max(vals(:)));
        m_c = min(c) - 1;
        m_r = min(r) - 1;
        linearInd = sub2ind([length(TH),length(R)], r, c);
        kern=[c- m_c,r- m_r, vals(linearInd)];
        ind = e + r_max + 1;
        
end