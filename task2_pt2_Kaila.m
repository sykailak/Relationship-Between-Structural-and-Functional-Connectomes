


%% Consider the variation in parameter estimates and fit quality
% across the different edges,and plot the relationship described by the 
% best-fitting models, for a small number of individual connections where 
% it holds particularly strongly or weakly.


%%


%Average SSE among number of measurements
mod4_sum_copy = mod4_sum./num_measurements{4};

% Exclude diagonal entries
mod4_sum_copy(logical(eye(68))) = NaN;

% Exclude entires where there are 1 or 0 measurements
mod4_sum_copy(find(num_measurements{4} ==1))=NaN;
mod4_sum_copy(find(num_measurements{4} ==0))=NaN;


min(mod4_sum_copy(:))
max(mod4_sum_copy(:))


[row,col] = find(mod4_sum_copy<0.0005 & mod4_sum_copy>0.0003)

func = [];
indstr = [];

% indices of 'good' fit quality
n = 17;
row(n)
col(n)

for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
                func = [func, f_data{t}(row(n),col(n))];
                indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

% f = a + B*t
figure;

subplot(2,2,1);
x = linspace(-1,1);
aval = a{4}(row(n),col(n));
bval = B{4}(row(n),col(n));
yval = Y{4}(row(n),col(n));
% Model 4 formula
plot(x, aval + bval.*x + yval.*x.^2)
hold on
scatter(indstr,func)
xlabel('Indirect Structural Connectivity') 
ylabel('Functional Connectivity')
title(sprintf('Fit of model 4 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','southwest')
ylim([-0.5 0.2])

func = [];
indstr = [];

% indices of 'good' fit quality
n = 16;
row(n)
col(n)

for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end


% f = a + B*t

subplot(2,2,3);
x = linspace(-1,1);
aval = a{4}(row(n),col(n));
bval = B{4}(row(n),col(n));
yval = Y{4}(row(n),col(n));
% Model 4 formula
plot(x, aval + bval.*x + yval.*x.^2)
hold on
scatter(indstr,func)
xlabel('Structural Connectivity') 
ylabel('Indirect Functional Connectivity')
title(sprintf('Fit of model 4 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','southwest')
ylim([-0.5 0.2])





row=[];
col=[];

[row,col] = find(mod4_sum_copy>0.006)


% indices of 'bad' fit quality
n = 1;
row(n)
col(n)

func=[];
indstr=[];

for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

subplot(2,2,2);

x = linspace(-1,1);
aval = a{4}(row(n),col(n));
bval = B{4}(row(n),col(n));
yval = Y{4}(row(n),col(n));
% Model 4 formula
plot(x, aval + bval.*x + yval.*x.^2)
hold on
scatter(indstr,func)
xlabel('Structural Connectivity') 
ylabel('Indirect Functional Connectivity')
title(sprintf('Fit of model 4 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','southwest')
ylim([-0.5 0.2])

% indices of 'bad' fit quality
n = 60;
row(n)
col(n)

func=[];
indstr=[];
for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

subplot(2,2,4);

x = linspace(-1,1);
aval = a{4}(row(n),col(n));
bval = B{4}(row(n),col(n));
yval = Y{4}(row(n),col(n));
% Model 4 formula
plot(x, aval + bval.*x + yval.*x.^2)
hold on
scatter(indstr,func)
xlabel('Structural Connectivity') 
ylabel('Indirect Functional Connectivity')
title(sprintf('Fit of model 4 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','southwest')
ylim([-0.5 0.2])


%%



%Average SSE among number of measurements
mod3_sum_copy = mod3_sum./num_measurements{3};

% Exclude diagonal entries
mod3_sum_copy(logical(eye(68))) = NaN;

% Exclude entires where there are 1 or 0 measurements
mod3_sum_copy(find(num_measurements{3} ==1))=NaN;
mod3_sum_copy(find(num_measurements{3} ==0))=NaN;


min(mod3_sum_copy(:))
max(mod3_sum_copy(:))



figure;
[row,col] = find(mod3_sum_copy<0.0005 & mod3_sum_copy>0.0003)

func = [];
indstr = [];

% indices of 'good' fit quality
n = 8;
row(n)
col(n)

for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

% f = a + B*t
figure;
subplot(2,1,1);
x = linspace(-1,1);
aval = a{3}(row(n),col(n));
bval = B{3}(row(n),col(n));
% Model 3 formula
plot(x, aval + bval.*x)

hold on
scatter(indstr,func)
xlabel('Indirect Structural Connectivity') 
ylabel('Functional Connectivity')
title(sprintf('Fit of model 3 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','northwest')
ylim([-0.5 0.2])


func = [];
indstr = [];

% indices of 'good' fit quality
n = 2;
row(n)
col(n)

for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

% f = a + B*t

subplot(2,1,2);
x = linspace(-1,1);
aval = a{3}(row(n),col(n));
bval = B{3}(row(n),col(n));
% Model 3 formula
plot(x, aval + bval.*x)
hold on
scatter(indstr,func)
xlabel('Indirect Structural Connectivity') 
ylabel('Functional Connectivity')
title(sprintf('Fit of model 3 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','northwest')
ylim([-0.5 0.2])


%%

figure;

row=[];
col=[];
[row,col] = find(mod3_sum_copy > 0.007)

% indices of 'bad' fit quality
n = 19;
row(n)
col(n)
func=[];
indstr=[];
for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

subplot(2,1,1);

x = linspace(-1,1);
aval = a{3}(row(n),col(n));
bval = B{3}(row(n),col(n));
% Model 3 formula
plot(x, aval + bval.*x)
hold on
scatter(indstr,func)
xlabel('Indirect Structural Connectivity') 
ylabel('Functional Connectivity')
title(sprintf('Fit of model 3 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','northwest')
ylim([-0.5 0.2])


% indices of 'bad' fit quality
n = 3;
row(n)
col(n)
func=[];
indstr=[];
for t = 1:19
    if ind_str_data{t}(row(n),col(n))~=0
        if f_data{t}(row(n),col(n))~=0
            func = [func, f_data{t}(row(n),col(n))];
            indstr = [indstr, ind_str_data{t}(row(n),col(n))];
        end
    end
end

subplot(2,1,2);

x = linspace(-1,1);
aval = a{3}(row(n),col(n));
bval = B{3}(row(n),col(n));
% Model 3 formula
plot(x, aval + bval.*x)
hold on
scatter(indstr,func)
xlabel('Indirect Structural Connectivity') 
ylabel('Functional Connectivity')
title(sprintf('Fit of model 3 to connectivity between region %.0f and %.0f', row(n),col(n)))
legend({'Fit','Data'},'Location','northwest')
ylim([-0.5 0.2])










