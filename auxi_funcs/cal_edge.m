function [edgeWeights] = cal_edge(sup_t1,lambda,t1_feature,t2_feature)
[h,w]   = size(sup_t1);
nbr_sp  = max(sup_t1(:));
idx_t1 = label2idx(sup_t1);
for i = 1:nbr_sp
    index_vector = idx_t1{i};
    [location_x location_y] = ind2sub(size(sup_t1),index_vector);
    location_center(i,:) = [round(mean(location_x)) round(mean(location_y))];
end
adj_mat = zeros(nbr_sp);
for i=2:h-1
    for j=2:w-1       
        label = sup_t1(i,j);        
        if (label ~= sup_t1(i+1,j-1))
            adj_mat(label, sup_t1(i+1,j-1)) = 1;
        end
        if (label ~= sup_t1(i,j+1))
            adj_mat(label, sup_t1(i,j+1)) = 1;
        end
        if (label ~= sup_t1(i+1,j))
            adj_mat(label, sup_t1(i+1,j)) = 1;
        end
        if (label ~= sup_t1(i+1,j+1))
            adj_mat(label, sup_t1(i+1,j+1)) = 1;
        end      
    end
end
adj_mat_1 = double((adj_mat + adj_mat')>0);
R = 2*round(sqrt(h*w/nbr_sp));
adj_mat = zeros(nbr_sp);
for i=1:nbr_sp
    for j = i:nbr_sp
    if ((location_center(i,1) - location_center(j,1))^2 + (location_center(i,2) - location_center(j,2))^2 < R^2)
     adj_mat (i,j) = 1;  
    end
    end
end
adj_mat = double((adj_mat + adj_mat')>0);
adj_mat_2 = adj_mat - eye(nbr_sp);
adj_mat = adj_mat_1|adj_mat_2;
%% edgeWeights
edgeWeights = zeros(sum(adj_mat(:)),4);
[node_x node_y] = find(adj_mat ==1);
edgeWeights(:,1) = node_x; % index of node 1
edgeWeights(:,2) = node_y; % index of node 2
for i = 1:sum(adj_mat(:))
    index_node_x = edgeWeights(i,1); 
    index_node_y = edgeWeights(i,2); 
    feature_t1_x = t1_feature(:,index_node_x);
    feature_t1_y = t1_feature(:,index_node_y);
    feature_t2_x = t2_feature(:,index_node_x);
    feature_t2_y = t2_feature(:,index_node_y);
    Dpq_t1(i) = norm(feature_t1_x-feature_t1_y,2)^2;
    Dpq_t2(i) = norm(feature_t2_x-feature_t2_y,2)^2;
    dist(i) = max(norm(location_center(index_node_x,:)-location_center(index_node_y,:),2),1);
end
sigma_t1 = mean(Dpq_t1);
sigma_t2 = mean(Dpq_t2);
for i =  1:sum(adj_mat(:))
    if Dpq_t1(i) <= sigma_t1 && Dpq_t2(i) <= sigma_t2
        Vpq(i) = exp(-Dpq_t1(i)/(2*sigma_t1)-Dpq_t2(i)/(2*sigma_t2));
    elseif Dpq_t1(i) <= sigma_t1 && Dpq_t2(i) > sigma_t2
        Vpq(i) = exp(Dpq_t1(i)/(2*sigma_t1)-1-Dpq_t2(i)/(2*sigma_t2));
    elseif Dpq_t1(i) > sigma_t1 && Dpq_t2(i) <= sigma_t2
        Vpq(i) = exp(-Dpq_t1(i)/(2*sigma_t1)+Dpq_t2(i)/(2*sigma_t2)-1);
    else
        Vpq(i) = exp(-1);
    end
end
Vpq = Vpq ./dist;
edgeWeights(:,3) = (1 - lambda)*Vpq;                  % node 1 ---> node 2
edgeWeights(:,4) = (1 - lambda)*Vpq;                  % node 2 ---> node 1