

%% Calculate an indirect structural connectivity matrix for each subject, 
% which is defined for our purposes as the greatest minimum weight in all 
% available two-step chains:


% Load structural connectivity of 19 healthy adults
str_data_dir = dir('connectomes-data/Task2Data/*_WFA_68.csv');

for t = 1:size(str_data_dir,1) %1~19
    s_matrix = readmatrix(str_data_dir(t).name);
    str_data{t} = s_matrix;
    len = length(s_matrix); %68
    t_matrix = zeros(len);
    for i = 1:len
        for j = 1:len 
            
            min_list = [];
            for k = 1:len
                % if s_ik, s_kj are BOTH nonzero 
                if s_matrix(i,k)~=0 & s_matrix(k,j)~=0
                    min_list = [min_list, min(s_matrix(i,k),s_matrix(k,j))];
                end
            end
            
            if length(min_list)~=0
                t_matrix(i,j) = max(min_list);  
            else
                %indirect structural connection doesn't exist
                t_matrix(i,j) = 0;
            end
        end
    end
    ind_str_data{t} = t_matrix;
end


%% Load functional data
f_data_dir = dir('connectomes-data/Task2Data/*_rsfMRI_68.csv');

for t = 1:size(f_data_dir,1) %1~19
    f_data{t} = readmatrix(f_data_dir(t).name);
end



%% Fit each of the following linear models to the data in turn.


% First model: f = a + B*s

a{1} = zeros(len);
B{1} = zeros(len);
num_measurements{1} = zeros(len);

for i = 1:len
    for j = 1:len
             
        func = [];
        structural = [];
        
        for t = 1:size(f_data,2); %1~19
            if str_data{t}(i,j) ~= 0
                if f_data{t}(i,j) ~=0
                    func(end+1) = f_data{t}(i,j);
                    structural(end+1) = str_data{t}(i,j);
                end
            end
        end
        
        if length(func) ~=0
            veclen = length(structural);
            design_mat = [ones(veclen,1),structural'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{1}(i,j) = X(1);
            B{1}(i,j) = X(2);
            num_measurements{1}(i,j) = veclen;

        %if either structural or functional doesn't exist, a{1}(i,j) and
        %b{1}(i,j) = 0

        end
        
    end
end

% calculate residuals
mod1_sum = zeros(68);

for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{1};
    B_val = B{1};
    %FILTER
    f_store(find(~str_data{t})) = 0;
    a_val(find(~str_data{t})) = 0;
    B_val(find(~str_data{t})) = 0;
    
    store_res = f_store - (a_val+B_val.*str_data{t});
    store_res = (store_res.^2);
    mod1_sum = mod1_sum+store_res;
    
end


%%
% Second model: f = a + B*s + Y*s^2

a{2} = zeros(len);
B{2} = zeros(len);
Y{2} = zeros(len);
num_measurements{2} = zeros(len);


for i = 1:len
    for j = 1:len
        func = [];
        structural = [];
        
        for t = 1:size(f_data,2);
            if str_data{t}(i,j) ~= 0
                if f_data{t}(i,j) ~=0
                    func(end+1) = f_data{t}(i,j);
                    structural(end+1) = str_data{t}(i,j);
                end
            end
        end
        
        if length(func)~=0
            veclen = length(structural);
            design_mat = [ones(veclen,1),structural',(structural.^2)'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{2}(i,j) = X(1);
            B{2}(i,j) = X(2);
            Y{2}(i,j) = X(3);
            num_measurements{2}(i,j) = veclen;
        end
    end
end

% calculate residuals
mod2_sum = zeros(68);
for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{2};
    B_val = B{2};
    Y_val = Y{2};
    %FILTER
    f_store(find(~str_data{t})) = 0;
    a_val(find(~str_data{t})) = 0;
    B_val(find(~str_data{t})) = 0;
    Y_val(find(~str_data{t})) = 0;
    
    store_res = f_store - (a_val + (B_val.*str_data{t}) + (Y_val.*(str_data{t}.^2)));
    store_res=(store_res.^2);
    mod2_sum = mod2_sum+store_res;
    
end


%% third model

% f = a + B*t

a{3} = zeros(len);
B{3} = zeros(len);
num_measurements{3} = zeros(len);

for i = 1:len
    for j = 1:len
        func = [];
        indstr = [];
        
        
        for t = 1:size(f_data,2);
            if ind_str_data{t}(i,j)~=0
                if f_data{t}(i,j)~=0
                    func(end+1) = f_data{t}(i,j);
                    indstr(end+1) = ind_str_data{t}(i,j);
                end
            end
        end
        
        if length(func)~=0
            veclen = length(indstr);
            design_mat = [ones(veclen,1),indstr'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{3}(i,j) = X(1);
            B{3}(i,j) = X(2);
            num_measurements{3}(i,j) = veclen;
        end
        
    end
end


% calculate residuals
mod3_sum = zeros(68);
for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{3};
    B_val = B{3};
    %FILTER
    f_store(find(~ind_str_data{t})) = 0;
    a_val(find(~ind_str_data{t})) = 0;
    B_val(find(~ind_str_data{t})) = 0;    
    
    store_res = f_store - (a_val + B_val.*ind_str_data{t});
    store_res = (store_res.^2);
    mod3_sum = mod3_sum+store_res;
end


%% Fourth model

% f = a + B*t + Y*t^2

a{4} = zeros(len);
B{4} = zeros(len);
Y{4} = zeros(len);
num_measurements{4} = zeros(len);

for i = 1:len
    for j = 1:len
        
        func = [];
        indstr = [];
        
        
        for t = 1:size(f_data,2);
            if ind_str_data{t}(i,j)~=0
                if f_data{t}(i,j)~=0
                    func(end+1) = f_data{t}(i,j);
                    indstr(end+1) = ind_str_data{t}(i,j);
                end
            end
            
        end
        
        if length(indstr)~=0
            veclen = length(indstr);
            design_mat = [ones(veclen,1),indstr',(indstr.^2)'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{4}(i,j) = X(1);
            B{4}(i,j) = X(2);
            Y{4}(i,j) = X(3);
            num_measurements{4}(i,j) = veclen;
            
        end
    end
end


% calculate residuals
mod4_sum = zeros(68);

for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{4};
    B_val = B{4};
    Y_val = Y{4};
    %FILTER
    f_store(find(~ind_str_data{t})) = 0;
    a_val(find(~ind_str_data{t})) = 0;
    B_val(find(~ind_str_data{t})) = 0;    
    Y_val(find(~ind_str_data{t})) = 0;
    
    store_res = f_store - (a_val + B_val.*ind_str_data{t} + Y_val.*(ind_str_data{t}.^2));
    store_res = (store_res.^2);
    mod4_sum = mod4_sum+store_res;
end

%% Fifth model

% f = a + B*s + Y*t

a{5} = zeros(len);
B{5} = zeros(len);
Y{5} = zeros(len);
num_measurements{5} = zeros(len);

for i = 1:len
    for j = 1:len
        
        func = [];
        indstr = [];
        structural = [];
        
        for t = 1:size(f_data,2);
            if str_data{t}(i,j)~=0
                if f_data{t}(i,j)~=0
                    if ind_str_data{t}(i,j)~=0
                        func(end+1) = f_data{t}(i,j);
                        structural(end+1) = str_data{t}(i,j);
                        indstr(end+1) = ind_str_data{t}(i,j);
                    end
                end
            end
        end
        
        if length(func)~=0
            veclen = length(func);
            design_mat = [ones(veclen,1),structural',indstr'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{5}(i,j) = X(1);
            B{5}(i,j) = X(2);
            Y{5}(i,j) = X(3);
            num_measurements{5}(i,j) = veclen;
            
        end
    end
end


% calculate residuals
mod5_sum = zeros(68);

for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{5};
    B_val = B{5};
    Y_val = Y{5};
    %FILTER
    f_store(find(~ind_str_data{t})) = 0;
    a_val(find(~ind_str_data{t})) = 0;
    B_val(find(~ind_str_data{t})) = 0;    
    Y_val(find(~ind_str_data{t})) = 0;   
    
    store_res = f_store - (a_val + B_val.*str_data{t} + Y_val.*ind_str_data{t});
    store_res = (store_res.^2);
    mod5_sum = mod5_sum +store_res;
end




%% For each model, evaluate the Akaike and Bayesian information criteria 
% (AIC and BIC), and perform leave-one-out (or k-fold) cross-validation to 
% evaluate its predictive value in terms of the sum-of-squared-errors.


% aic formula
aic{1} = 2*3 + (num_measurements{1} .* log(  (1./num_measurements{1}).*mod1_sum));
%Exclude cells which would output undefined
aic{1}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;

%Exclude diagonal entries
aic{1}(find(eye(68)))=NaN;
aic_val{1} = nanmean(aic{1},'all');


aic{2} = 2*4 + (num_measurements{2} .* log(  (1./num_measurements{2}).*mod2_sum));
aic{2}(find(num_measurements{2} == 0 | num_measurements{2}==1))=NaN;
aic{2}(find(eye(68)))=NaN;
aic_val{2} = nanmean(aic{2},'all');


aic{3} = 2*3 + (num_measurements{3} .* log(  (1./num_measurements{3}).*mod3_sum));
aic{3}(find(num_measurements{3} == 0 | num_measurements{3}==1))=NaN;
% These indices have mod3_sum = 0 (residual=0). Therefore output -inf.
aic{3}(find(~isfinite(aic{3})))=NaN;
aic{3}(find(eye(68)))=NaN;
aic_val{3} = nanmean(aic{3},'all');

aic{4} = 2*4 + (num_measurements{4} .* log(  (1./num_measurements{4}).*mod4_sum));
aic{4}(find(num_measurements{4} == 0 | num_measurements{4}==1))=NaN;
aic{4}(find(eye(68)))=NaN;
aic_val{4} = nanmean(aic{4},'all');

aic{5} = 2*4 + (num_measurements{5} .* log(  (1./num_measurements{5}).*mod5_sum));
aic{5}(find(num_measurements{5} == 0 | num_measurements{5}==1))=NaN;
aic{5}(find(eye(68)))=NaN;
aic_val{5} = nanmean(aic{5},'all')

%%
%bic (lowest is best)




bic{1} = 3*log(num_measurements{1}) + num_measurements{1}.*log((1./num_measurements{1}).*mod1_sum);
%Exclude cells which would output undefined
bic{1}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;
bic{1}(find(eye(68)))=NaN;
bic_val{1} = nanmean(bic{1},'all');

bic{2} = 4*log(num_measurements{2}) + num_measurements{2}.*log((1./num_measurements{2}).*mod2_sum);
bic{2}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;
bic{2}(find(eye(68)))=NaN;
bic_val{2} = nanmean(bic{2},'all');

bic{3} = 3*log(num_measurements{3}) + num_measurements{3}.*log((1./num_measurements{3}).*mod3_sum);
bic{3}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;
bic{3}(find(eye(68)))=NaN;
bic_val{3} = nanmean(bic{3},'all');

bic{4} = 4*log(num_measurements{4}) + num_measurements{4}.*log((1./num_measurements{4}).*mod4_sum);
bic{4}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;
bic{4}(find(eye(68)))=NaN;
bic_val{4} = nanmean(bic{4},'all');


bic{5} = 4*log(num_measurements{5}) + num_measurements{5}.*log((1./num_measurements{5}).*mod5_sum);
bic{5}(find(num_measurements{1} == 0 | num_measurements{1}==1))=NaN;
bic{5}(find(eye(68)))=NaN;
bic_val{5} = nanmean(bic{5},'all')





%% Cross Validation
% LOOCV
% First model: f = a + B*s


cv_residuals{1} = 0;

for cv = 1:19
    % partition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition structural data
    str_test = str_data(cv);
    str_train = str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{1} = zeros(len); %len = 68
    cv_B{1} = zeros(len);
    
    for i = 1:len %68
        for j = 1:len %68
            func = [];
            structural = [];
            for t = 1:18
                if str_train{t}(i,j)~=0
                    if f_train{t}(i,j)~=0
                        func(end+1) = f_train{t}(i,j);
                        structural(end+1) = str_train{t}(i,j);
                    end
                end
            end
            
            if length(func)~=0
                veclen=length(func);
                design_mat = [ones(veclen,1),structural'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{1}(i,j) = X(1);
                cv_B{1}(i,j) = X(2);
            end
        end
    end
    
    % All str_test zero elements will be NaN
    f_test{1}(find(~str_test{1})) = NaN;
    store_res = f_test{1} - (cv_a{1} + cv_B{1}.*str_test{1});
    store_res = store_res.^2;
    cv_residuals{1} = cv_residuals{1} + nanmean(store_res,'all');

end

%% Cross Validation
% Second model: f = a + B*s + Y*s^2
 
cv_residuals{2} = 0;

for cv = 1:19
    % particition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition structural data
    str_test = str_data(cv); % not used
    str_train = str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{2} = zeros(len); %len = 68
    cv_B{2} = zeros(len);
    cv_Y{2} = zeros(len);
    
    for i = 1:len %68
        for j = 1:len %68
            func = [];
            structural = [];
            for t = 1:18
                if str_train{t}(i,j)~=0
                    if f_train{t}(i,j)~=0
                        func(end+1) = f_train{t}(i,j);
                        structural(end+1) = str_train{t}(i,j);
                    end
                end
            end
            
            if length(func)~=0
                veclen = length(func);
                design_mat = [ones(veclen,1),structural',(structural.^2)'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{2}(i,j) = X(1);
                cv_B{2}(i,j) = X(2);
                cv_Y{2}(i,j) = X(3);
            end
        end
    end
    
    f_test{1}(find(~str_test{1})) = NaN;
    
    store_res = f_test{1} - (cv_a{2} + cv_B{2}.*str_test{1} + cv_Y{2}.*(str_test{1}.^2));
    store_res = store_res.^2;
    cv_residuals{2} = cv_residuals{2}+ nanmean(store_res,'all');

end

%% Cross Validation
% third model: f = a + B*t


cv_residuals{3} = 0;

for cv = 1:19
    % particition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition indirect structural data
    indstr_test = ind_str_data(cv); % not used
    indstr_train = ind_str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{3} = zeros(len); %len = 68
    cv_B{3} = zeros(len);
    
    for i = 1:len %68
        for j = 1:len %68
            func = [];
            ind_structural = [];
            for t = 1:18
                if indstr_train{t}(i,j)~=0
                    if f_train{t}(i,j)~=0
                        func(end+1) = f_train{t}(i,j);
                        ind_structural(end+1) = indstr_train{t}(i,j);
                    end
                end
                
            end
            
            if length(func)~=0
                veclen = length(func);
                design_mat = [ones(veclen,1),ind_structural'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{3}(i,j) = X(1);
                cv_B{3}(i,j) = X(2);
            end
            
        end
    end
    
    f_test{1}(find(~indstr_test{1}))=NaN;
    
    store_res = f_test{1} - (cv_a{3} + cv_B{3}.*indstr_test{1});
    store_res = store_res.^2;
    
    cv_residuals{3} = cv_residuals{3} + nanmean(store_res,'all');
    
end




%% fourth model

% f = a + B*t + Y*t^2

cv_residuals{4} = 0;

for cv = 1:19
    % particition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition indirect structural data
    indstr_test = ind_str_data(cv); % not used
    indstr_train = ind_str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{4} = zeros(len); %len = 68
    cv_B{4} = zeros(len);
    cv_Y{4} = zeros(len);
    for i = 1:len %68
        for j = 1:len %68
            func = [];
            ind_structural = [];
            for t = 1:18
                if indstr_train{t}(i,j)~=0
                    if f_train{t}(i,j)~=0
                        func(end+1) = f_train{t}(i,j);
                        ind_structural(end+1) = indstr_train{t}(i,j);
                    end
                end
                
                
            end
            
            if length(func)~=0
                veclen = length(func);
                design_mat = [ones(veclen,1),ind_structural',(ind_structural.^2)'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{4}(i,j) = X(1);
                cv_B{4}(i,j) = X(2);
                cv_Y{4}(i,j) = X(3);
            end
            
        end
    end
    
    f_test{1}(find(~indstr_test{1}))=NaN;
    
    store_res = f_test{1} - (cv_a{4} + cv_B{4}.*indstr_test{1}+ cv_Y{4}.*(indstr_test{1}.^2));
    store_res = store_res.^2;
    
    cv_residuals{4} = cv_residuals{4} + nanmean(store_res,'all');
    
end





%% fifth model

% f = a + B*s + Y*t

cv_residuals{5} = 0;

for cv = 1:19
    % particition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition structural data
    str_test = str_data(cv); % not used
    str_train = str_data([1:cv-1,cv+1:19]);
    % particition indirect structural data    
    indstr_test = ind_str_data(cv); % not used
    indstr_train = ind_str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{5} = zeros(len); %len = 68
    cv_B{5} = zeros(len);
    cv_Y{5} = zeros(len);
    for i = 1:len %68
        for j = 1:len %68
            func = [];
            structural = [];
            ind_structural = [];
            for t = 1:18
                if f_train{t}(i,j)~=0 
                    if str_train{t}(i,j)~=0
                        if indstr_train{t}(i,j)~=0
                            func(end+1) = f_train{t}(i,j);
                            structural(end+1) = str_train{t}(i,j);
                            ind_structural(end+1) = indstr_train{t}(i,j);
                        end
                    end
                end
                
                
            end
            
            if length(func)~=0
                veclen = length(func);
                design_mat = [ones(veclen,1),structural',ind_structural'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{5}(i,j) = X(1);
                cv_B{5}(i,j) = X(2);
                cv_Y{5}(i,j) = X(3);
            end
            
        end
    end
    
    f_test{1}(find(~str_test{1})) = NaN;
    f_test{1}(find(~indstr_test{1}))=NaN;
    
    store_res = f_test{1} - (cv_a{5} + cv_B{5}.*str_test{1} + cv_Y{5}.*indstr_test{1});
    store_res = store_res.^2;
    cv_residuals{5} = cv_residuals{5} + nanmean(store_res,'all');

end


%%
cv_residuals





%%
subject = 9;
figure;
subplot(1,3,1)
imagesc(f_data{subject}); colorbar; colormap(jet(10));
title('Functional','FontSize',15)

subplot(1,3,2)
imagesc(str_data{subject}); colorbar
title('Structural','FontSize',15)

subplot(1,3,3)
imagesc(ind_str_data{subject}); colorbar
title('Indirect structural','FontSize',15)






%% negative correlation vs positive correlation












