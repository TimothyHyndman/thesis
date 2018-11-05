function Z = calculate_flag_graph_n2(phi, width, res)

    y1 = linspace(-width,width,res);
    y2 = linspace(-width,width,res);

    y1 = [y1,y1(res)+y1(res) - y1(res-1)];  %*awful hack so we display final row and column of data in plot
    y2 = [y2,y2(res)+y2(res) - y2(res-1)];  %*
    [Y1,Y2] = meshgrid(y1,y2);

    Z = zeros(res+1)+3;   %*
    Z_likelihood = zeros(res+1);
    predictor = Z_likelihood;

    for i = 1:res
    % a = 1:res;
    % ind = randperm(length(a));
    % for i = a(ind)
        for j = 1:min(i,res - i+1)    %For n=3 problems
    %     for j = 1:i %For n=4 problems
    %     b = 1:i;
    %     ind = randperm(length(b));
    %     for j = b(ind)
            
            %For n = 3;
            Y = [Y1(i,j);Y2(i,j)];
            if 1
    %             [solution_grid_points,solution_masses,likelihood] = MixtureLikelihoodMovingMasses(phi,Y);
                [Q,likelihood] = MixtureLikelihoodMovingMasses2(phi,Y);
                solution_masses = Q.ProbWeights;
                solution_grid_points = Q.Support;

                Z(i,j) = size(solution_masses,2);
                Z_likelihood(i,j) = likelihood;
            end
            
            %For n = 2;
    %         Y = [Y1(i,j);Y2(i,j)];
            if 0
                if abs(Y(1) - Y(2)) < 2*sigma
                    Z(i,j) = 1;
                else
                    Z(i,j) = 2;
                end
            end
            
    %         F(i,j) = exp(-3*log(2*sigma) - (1/(sigma^2))*(sum((Y - mean(Y)).^2)));
            
    %         if sum((Y - mean(Y)).^2) < 3
    %             predictor(i,j) = 1;
    %         else
    %             predictor(i,j) = 0;
    %         end
    %         if max(Y) - min(Y) < 2
    %             Z(i,j) = 1;
    %         else
    %             Z(i,j) = 2;
    %         end

    %         Z(i,j) = max(Y) - min(Y);

            Z(j,i) = Z(i,j);
            Z(res - i + 1,res - j + 1) = Z(i,j);
            Z(res - j + 1,res - i + 1) = Z(i,j);
            
    %         Z_likelihood(j,i) = Z_likelihood(i,j);
    %         Z_likelihood(res - i + 1,res - j + 1) = Z_likelihood(i,j);
    %         Z_likelihood(res - j + 1,res - i + 1) = Z_likelihood(i,j);
            
    %         predictor(j,i) = predictor(i,j);
    %         predictor(res - i + 1,res - j + 1) = predictor(i,j);
    %         predictor(res - j + 1,res - i + 1) = predictor(i,j);
        end
        if mod(i,4) == 1
            pcolor(Y1,Y2,Z)
            shading flat
    %         pcolor(Y1,Y2,Z_likelihood)
            drawnow
        end
    end
end