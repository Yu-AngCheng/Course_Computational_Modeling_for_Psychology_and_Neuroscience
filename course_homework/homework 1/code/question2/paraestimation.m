function paramsEst = paraestimation(data)
LB = [-Inf, -Inf];
UB = [Inf, Inf];
x0 = [45, 0.1];
intensity_level=1:100;
paramsEst = fminsearchbnd(@(params)mll(params, intensity_level, data), x0, LB, UB);

end

