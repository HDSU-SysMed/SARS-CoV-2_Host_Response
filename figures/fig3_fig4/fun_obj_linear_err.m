function r = fun_obj_linear_err(log_b,X,Y)  

    b = exp(log_b);
   
    y_mod = b(1) * abs(X) + b(2)*max(abs(X));
    
    % lsqnonlin:
    r = Y - y_mod; 
end