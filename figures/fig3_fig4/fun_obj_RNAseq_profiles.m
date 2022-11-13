function r = fun_obj_RNAseq_profiles(log_b,X,Y,fun_type)  

    t = X;
    y_obs = Y(:,1);
    e_obs = Y(:,2);

    if fun_type==1 % sigmoidal increase
        b = [exp(log_b(1:2)) log_b(3)];
        y_mod = b(1)./(1 + exp(-b(2)*(t-b(3))));
    elseif fun_type==2 % sigmoidal increase before sigmoidal decrease
        b = [exp(log_b(1:2)) log_b(3) exp(log_b(4)) log_b(5)];
        y_mod = b(1)./(1 + exp(-b(2)*(t-b(3)))) - ...
            b(4)./(1 + exp(-b(2)*(t-(b(3)+b(5)))));
    elseif fun_type==3 % sigmoidal decrease
        b = [exp(log_b(1:2)) log_b(3)];
        y_mod = -b(1)./(1 + exp(-b(2)*(t-b(3))));
    elseif fun_type==4 % sigmoidal decrease before increase        
        b = [exp(log_b(1:2)) log_b(3) exp(log_b(4)) log_b(5)];
        y_mod = -b(1)./(1 + exp(-b(2)*(t-b(3)))) + ...
            b(4)./(1 + exp(-b(2)*(t-(b(3)+b(5)))));
    end
    
    % lsqnonlin:
    r = (y_obs - y_mod)./e_obs; 
end