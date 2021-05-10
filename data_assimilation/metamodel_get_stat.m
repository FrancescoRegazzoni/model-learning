function [mean_val,stdv_val] = metamodel_get_stat(metamod,min_field,max_field,get_random_field,num)

    if isfield(metamod,min_field) && isfield(metamod,max_field) 
        mean_val = (eval(['metamod.' max_field]) + eval(['metamod.' min_field]))/2;
        stdv_val = (eval(['metamod.' max_field]) - eval(['metamod.' min_field]))/sqrt(12); % std of uniform distribution
    else
        num_samples = 1e2;
        X = zeros(num,num_samples);
        for i = 1:num_samples
            X(:,i) = eval(['metamod.' get_random_field '()']);
        end
        mean_val = mean(X,2);
        stdv_val = std(X,1,2);
    end
    
    mean_val = adapt_dimension(mean_val,num);
    stdv_val = adapt_dimension(stdv_val,num);
    
end