

%% Repeat the previous subtask, but using a single set of coefficients 
%%for all edges, for each model.


%% First model: f = a + B*s

func = [];
structural = [];
indstructural = [];

for t = 1:19
    % make copy
    fdata = f_data{t};
    strdata = str_data{t};
    
    %exclude edges where str_data =0
    fdata(find(~strdata))=[];
    strdata(find(~strdata))=[];
    
    func=[func,fdata];
    structural=[structural,strdata];
end



veclen = length(func);
design_mat = [ones(veclen,1),structural'];

X = pinv(design_mat'*design_mat)*design_mat'*func';
all_a{1} = X(1);
all_B{1} = X(2);

% calculate residual
all_residual{1} = func - (all_a{1} + all_B{1}.*structural);
all_residual{1} = sum(all_residual{1}.^2);

%MSE
%all_residual{1} = all_residual{1}/veclen;

%AIC BIC
k=veclen;
all_aic{1} = 2*3 + k*log((1/k)*all_residual{1});
all_bic{1} = 3*log(k) + k*log((1/k)*all_residual{1});

%% Second model: f = a + B*s + Y*s^2


design_mat = [ones(veclen,1),structural',(structural.^2)'];
X = pinv(design_mat'*design_mat)*design_mat'*func';
all_a{2} = X(1);
all_B{2} = X(2);
all_Y{2} = X(3);

% calculate residual
all_residual{2} = func - (all_a{2} +   all_B{2}.*structural +   all_Y{2}.*(structural.^2)   );
all_residual{2} = sum(all_residual{2}.^2);
%MSE
%all_residual{2} = all_residual{2}/veclen;

%AIC BIC
k=veclen;
all_aic{2} = 2*4 + k*log((1/k)*all_residual{2});
all_bic{2} = 4*log(k) + k*log((1/k)*all_residual{2});


%% third model: f = a + B*t

func = [];
structural = [];
indstructural = [];

for t = 1:19
    % make copy
    fdata = f_data{t};
    indstrdata = ind_str_data{t};
    
    %exclude edges where ind_str_data =0
    fdata(find(~indstrdata))=[];
    indstrdata(find(~indstrdata))=[];

    func=[func,fdata];
    indstructural=[indstructural,indstrdata];
end

veclen = length(func);
design_mat = [ones(veclen,1),indstructural'];
X = pinv(design_mat'*design_mat)*design_mat'*func';
all_a{3} = X(1);
all_B{3} = X(2);


% calculate residual
all_residual{3} = func - (all_a{3} +   all_B{3}.*indstructural);
all_residual{3} = sum(all_residual{3}.^2);
%MSE
%all_residual{3} = all_residual{3}/veclen;

%AIC BIC
k=veclen;
all_aic{3} = 2*3 + k*log((1/k)*all_residual{3});
all_bic{3} = 3*log(k) + k*log((1/k)*all_residual{3});

%% Fourth model: f = a + B*t + Y*t^2

design_mat = [ones(veclen,1),indstructural',(indstructural.^2)'];

X = pinv(design_mat'*design_mat)*design_mat'*func';
all_a{4} = X(1);
all_B{4} = X(2);
all_Y{4} = X(3);

% calculate residual
all_residual{4} = func - (  all_a{4} +   all_B{4}.*indstructural    + all_Y{4}.*(indstructural.^2)  );
all_residual{4} = sum(all_residual{4}.^2);
%MSE
%all_residual{4} = all_residual{4} /veclen;

%AIC BIC
k=veclen;
all_aic{4} = 2*4 + k*log((1/k)*all_residual{4});
all_bic{4} = 4*log(k) + k*log((1/k)*all_residual{4});

%% Fifth model: f = a + B*s + Y*t



func = [];
structural = [];
indstructural = [];

for t = 1:19
    % make copy
    fdata = f_data{t};
    strdata = str_data{t};
    indstrdata = ind_str_data{t};
    
    %exclude edges where str_data =0
    fdata(find(~strdata))=[];
    indstrdata(find(~strdata))=[];
    strdata(find(~strdata))=[];

    func=[func,fdata];
    structural = [structural, strdata];
    indstructural=[indstructural,indstrdata];
end

veclen = length(func);

design_mat = [ones(veclen,1),structural',indstructural'];
X = pinv(design_mat'*design_mat)*design_mat'*func';
all_a{5} = X(1);
all_B{5} = X(2);
all_Y{5} = X(3);


% calculate residual
all_residual{5} = func - (  all_a{5} +   all_B{5}.*structural    + all_Y{5}.*indstructural  );
all_residual{5} = sum(all_residual{5}.^2);
%all_residual{5} = all_residual{5}/veclen;

k=veclen;
all_aic{5} = 2*4 + k*log((1/k)*all_residual{5});
all_bic{5} = 4*log(k) + k*log((1/k)*all_residual{5});

all_aic
all_bic










%% ADDITIONAL MODEL

figure;

for i = 1:8
    func=[];
    indstr=[];
    
    [row,col]=find(num_measurements{3}>15);
    n = randi([1 68]);

    for t = 1:19
        if ind_str_data{t}(row(n),col(n))~=0
            if f_data{t}(row(n),col(n))~=0
                func = [func, f_data{t}(row(n),col(n))];
                indstr = [indstr, ind_str_data{t}(row(n),col(n))];
            end
        end
    end
    subplot(4,4,i);
    
    x = linspace(-1,1);
    aval = a{4}(row(n),col(n));
    bval = B{4}(row(n),col(n));
    yval = Y{4}(row(n),col(n));
    % Model 4 formula
    plot(x, aval + bval.*x + yval.*x.^2)
    hold on
    
    scatter(indstr,func)
    title(sprintf('%.0f and %.0f', row(n),col(n)))
    ylim([-0.1 0.2])
    xlim([0.2 0.5])
    xlabel('Indirect Structural Connectivity') 


end
%sgtitle('Indirect structural')






for i = 9:16
    func=[];
    str=[];
    indstr=[];
    
    [row,col]=find(num_measurements{2}>15);
    n = randi([1 68]);

    for t = 1:19
        if str_data{t}(row(n),col(n))~=0
            if f_data{t}(row(n),col(n))~=0
                if ind_str_data{t}(row(n),col(n))~=0
                    func = [func, f_data{t}(row(n),col(n))];
                    str = [str, str_data{t}(row(n),col(n))];
                    indstr = [indstr, ind_str_data{t}(row(n),col(n))];
                end
              
            end
        end
    end
    subplot(4,4,i);
    
    x = linspace(-1,1);
    y = linspace(-1,1);
    
    aval = a{6}(row(n),col(n));
    bval = B{6}(row(n),col(n));
    yval = Y{6}(row(n),col(n));


    % Model formula
    plot3(x,y, aval + bval.*(x.^2) + yval.*(y.^2))
    hold on
    scatter3(str,indstr,func)
    title(sprintf('%.0f and %.0f', row(n),col(n)))
    zlim([-0.1 0.2])
    xlim([0.2 0.5])
    ylim([0.2 0.5])
     


end


%% Model 6

% f = a +B*

a{6} = zeros(len);
B{6} = zeros(len);
Y{6} = zeros(len);
num_measurements{6} = zeros(len);


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
            design_mat = [ones(veclen,1),(structural.^2)',(indstr.^2)'];
            X = pinv(design_mat)*func';
            a{6}(i,j) = X(1);
            B{6}(i,j) = X(2);
            Y{6}(i,j) = X(3);

            num_measurements{6}(i,j) = veclen;
            
        end
    end
end


% calculate residuals
mod6_sum = zeros(68);

for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{6};
    B_val = B{6};
    Y_val = Y{6};

    %FILTER
    f_store(find(~ind_str_data{t})) = 0;
    a_val(find(~ind_str_data{t})) = 0;
    B_val(find(~ind_str_data{t})) = 0;    
    Y_val(find(~ind_str_data{t})) = 0; 

    
    store_res = f_store - (a_val + B_val.*(str_data{t}.^2) + Y_val.*(ind_str_data{t}.^2));
    store_res = (store_res.^2);
    mod6_sum = mod6_sum +store_res;
end



aic{6} = 2*4 + (num_measurements{6} .* log(  (1./num_measurements{6}).*mod6_sum));
aic{6}(find(num_measurements{6} == 0 | num_measurements{6}==1))=NaN;
aic{6}(find(eye(68)))=NaN;
aic_val{6} = nanmean(aic{6},'all')



bic{6} = 4*log(num_measurements{6}) + num_measurements{6}.*log((1./num_measurements{6}).*mod6_sum);
bic{6}(find(num_measurements{6} == 0 | num_measurements{6}==1))=NaN;
bic{6}(find(eye(68)))=NaN;
bic_val{6} = nanmean(bic{6},'all')


%% Model 7

% f = a +B*t^3

a{7} = zeros(len);
B{7} = zeros(len);
num_measurements{7} = zeros(len);


for i = 1:len
    for j = 1:len
        
        func = [];
        indstr = [];
        
        for t = 1:size(f_data,2);
            if f_data{t}(i,j)~=0
                if ind_str_data{t}(i,j)~=0
                    func(end+1) = f_data{t}(i,j);
                    indstr(end+1) = ind_str_data{t}(i,j);
                end
            end
        end
        
        if length(func)~=0
            veclen = length(func);
            design_mat = [ones(veclen,1),(indstr.^3)'];
            X = pinv(design_mat'*design_mat)*design_mat'*func';
            a{7}(i,j) = X(1);
            B{7}(i,j) = X(2);

            num_measurements{7}(i,j) = veclen;
            
        end
    end
end


% calculate residuals
mod7_sum = zeros(68);

for t = 1:19
    %COPY
    f_store = f_data{t};
    a_val = a{7};
    B_val = B{7};

    %FILTER
    f_store(find(~ind_str_data{t})) = 0;
    a_val(find(~ind_str_data{t})) = 0;
    B_val(find(~ind_str_data{t})) = 0;    

    
    store_res = f_store - (a_val + B_val.*(ind_str_data{t}.^3));
    store_res = (store_res.^2);
    mod7_sum = mod7_sum +store_res;
end



aic{7} = 2*3 + (num_measurements{7} .* log(  (1./num_measurements{7}).*mod7_sum));
aic{7}(find(num_measurements{7} == 0 | num_measurements{7}==1))=NaN;
aic{7}(find(eye(68)))=NaN;
aic_val{7} = nanmean(aic{7},'all')



bic{7} = 3*log(num_measurements{7}) + num_measurements{7}.*log((1./num_measurements{7}).*mod7_sum);
bic{7}(find(num_measurements{7} == 0 | num_measurements{7}==1))=NaN;
bic{7}(find(eye(68)))=NaN;
bic_val{7} = nanmean(bic{7},'all')



%%


% 7th model: f = a + B*t^3


cv_residuals{7} = 0;

for cv = 1:19
    % particition functional data
    f_test=f_data(cv);
    f_train=f_data([1:cv-1,cv+1:19]);
    % particition indirect structural data
    indstr_test = ind_str_data(cv); % not used
    indstr_train = ind_str_data([1:cv-1,cv+1:19]);
    
    % create placeholders for model parameters
    cv_a{7} = zeros(len); %len = 68
    cv_B{7} = zeros(len);
    
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
                design_mat = [ones(veclen,1),(ind_structural.^3)'];
                X = pinv(design_mat'*design_mat)*design_mat'*func';
                cv_a{7}(i,j) = X(1);
                cv_B{7}(i,j) = X(2);
            end
            
        end
    end
    
    f_test{1}(find(~indstr_test{1}))=NaN;
    
    store_res = f_test{1} - (cv_a{7} + cv_B{7}.*(indstr_test{1}.^3));
    store_res = store_res.^2;
    
    cv_residuals{7} = cv_residuals{7} + nanmean(store_res,'all');
    
end

cv_residuals{7}





%% ADDITIONAL MODEL

figure;
[row,col] = find(num_measurements{7}>10);

for i = 1:16
    func=[];
    indstr=[];
    
    n = randi([1 4102])

    for t = 1:19
        if ind_str_data{t}(row(n),col(n))~=0
            if f_data{t}(row(n),col(n))~=0
                func = [func, f_data{t}(row(n),col(n))];
                indstr = [indstr, ind_str_data{t}(row(n),col(n))];
            end
        end
    end
    subplot(4,4,i);
    
    x = linspace(0,0.7);
    aval = a{7}(row(n),col(n));
    bval = B{7}(row(n),col(n));
    % Model 4 formula
    plot(x, aval + bval.*(x.^3))
    hold on
    
    scatter(indstr,func)
    title(sprintf('%.0f and %.0f', row(n),col(n)))
    %ylim([-0.1 0.2])
    xlim([0.2 0.7])


end
sgtitle({'f = a + βt^3','x-axis: indirect structural, y-axis: functional'})


%%
figure;
histogram(B{1}(:),10)

%%

