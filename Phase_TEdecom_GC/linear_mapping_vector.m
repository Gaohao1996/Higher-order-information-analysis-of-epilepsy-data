function X_2D = linear_mapping_vector(X_1D)
    X1 = X_1D.*2-1;
    X2 = X_1D.*2;
    X_2D =  sort([X1 X2]);
end

