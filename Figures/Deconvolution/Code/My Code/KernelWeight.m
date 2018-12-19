function weight = KernelWeight(weight_type, x)
    length_x = length(x);
    switch weight_type
        case 'Epanechnikov'
            sig_x=-x(1)/2;
            weight=0.75/(2*sig_x)*(1-(x/(2*sig_x)).^2);
        case 'Uniform'
            max_x=-x(1);
            weight = zeros(1,length_x) + 1/(2*max_x);
        case 'Triangular'
            max_x = -x(1);
            weight = -1/(max_x^2)*(abs(x)-max_x);
        case 'Triweight'
            max_x = -x(1);
            weight = 35/(32*max_x)*(1 - (x/max_x).^2).^3;
    end
end