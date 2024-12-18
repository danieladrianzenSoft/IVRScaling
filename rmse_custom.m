function r = rmse_custom (t,values_a,values_b,varargin)
    p = inputParser;
    p.addRequired('t', @(x) length(x)>=1);
    p.addRequired('values_a', @(x) length(x)>=1);
    p.addRequired('values_b', @(x) length(x)>=1);
    p.addOptional('time_pt', length(t), @isscalar);
    p.parse(t,values_a,values_b,varargin{:});

    time_pt = p.Results.time_pt;
    
    i_time_pt = find(t >= time_pt,1,'first');
    t_vec = t(1:i_time_pt);
    a_vec = values_a(1:i_time_pt);
    b_vec = values_b(1:i_time_pt);
    
    r = sqrt(trapz(t_vec,(a_vec(:)-b_vec(:)).^2/time_pt));
end 