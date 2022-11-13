function ind = fun_find_first_and_second_thd(gene_count_norm)
   
    ind_2=nan;
    ind_x_max = find(gene_count_norm==max(gene_count_norm),1,'first');
    part_1 = gene_count_norm(1:ind_x_max);
    part_2 = gene_count_norm((ind_x_max+1):end);
    
    % ind smaller than maximum
    ind1_eq_50 = find(part_1==.5);
    ind1_ex_50 = find(part_1>.5,1,'first')-.5;
    if isempty(ind1_eq_50)
        ind_1 = ind1_ex_50;
    else
        ind_1 = mean(ind1_eq_50);
    end

    % ind larger than maximum
    part_2_inv = fliplr(part_2);
    ind2_eq_50 = find(part_2_inv==.5);
    ind2_ex_50 = find(part_2_inv>.5,1,'first')-.5;
    if isempty(ind1_eq_50)&&(ind2_ex_50>0.5)
        ind_2 = ind2_ex_50;
    elseif ~isempty(ind1_eq_50)&&(ind2_ex_50>0.5)
        ind_2 = mean(ind2_eq_50);
    end
    ind = [ind_1, length(gene_count_norm)-ind_2+1];
    
end
