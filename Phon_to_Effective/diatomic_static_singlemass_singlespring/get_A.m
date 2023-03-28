function A = get_A(num_nodes)
% Get matrix A to construct final Stiffness Matrix

A = zeros(num_nodes,num_nodes);
for i = 1:num_nodes
    for j = 1:num_nodes
        if i==j
            A(i,j) = 1;
        end
        if j==i-1
            A(i,j) = -1;
        end
    end
end

end


