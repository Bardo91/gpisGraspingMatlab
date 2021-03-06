function covMatFull = DiagCovMatFullTP(R,X)
%% Description
% Computes the covariance matrix between function values and 1st order
% derivatives between input locations X1 and X2
% given GP parameters sigma, gamma
% X1 = X2 returns the covariance matrix between function values and 1st order
% derivatives at X1

%% Inputs
% R
% X

%% Outputs
% KderX1X2 - Covariance Matrix


[D,N1] = size(X);

covMatFull = zeros(N1 * (D + 1));

for n1 = 1 : N1
   n2 = n1;
        dist = X(:,n1) - X(:,n2);
        distNorm = norm(dist);
        for d1 = 0 : D
            d2 = d1;
                if ((n1 - 1) * (D + 1) + d1 + 1 <= ...
                        (n2 - 1) * (D + 1) + d2 + 1)
                    if (d1 == 0) && (d2 == 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1) = ...
                            2 * distNorm^3 ...
                            - 3 * R * distNorm^2 + R^3;
                    elseif (d1 > 0) && (d2 == 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1) = ...
                            6 * dist(d1) * (distNorm - R);
                    elseif (d1 == 0) && (d2 > 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1) = ...
                            - 6 * dist(d2) * (distNorm - R);
                    else
                        if n1 == n2
                            covMatFull((n1 - 1) * (D + 1) + d1 + 1) = ...
                                - 6 * ((d1 == d2) * (distNorm - R));
                        else
                            covMatFull((n1 - 1) * (D + 1) + d1 + 1) = ...
                                - 6 * (dist(d1) * dist(d2)/distNorm ...
                                + (d1 == d2) * (distNorm - R));
                        end
                    end
        end
    end
end
covMatFull = covMatFull + covMatFull' - diag(diag(covMatFull));

