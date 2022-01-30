%% PROJECT 3: Car Application 
%
% AUTHOR: Anders Munk-Nielsen, University of Copenhagen 
% DATE: 2020-11-17
%
% SOURCES: Full credit to the amazing people that have collected the
% underlying datasets: 
% Frank Verboven's car dataset 
%   --> https://sites.google.com/site/frankverbo/data-and-software/data-set-on-the-european-car-market
% the World Inequality Database 
%   --> https://wid.world/

clc; clear all; 

% Setting: whether to drop any car types that were only available in too
% few (market, period)-pairs, since this can cause some problems with
% identifying a dummy for that car type or firm. 
DROPSINGLEYEARCARS = true; 
minyearobs = 3; % if above is true, this is the smallest required number of (market, periods)

% call the script that reads data 
read_data; 

%% Price variable 
% choose the price variable you wish to work with: 

%Cars.logpr = log(Cars.pr); 
Cars.logpr = log(Cars.eurpr); 
%Cars.logpr = log(Cars.princ); 

%% Car Attributes 
% Convenient to make a list of variables that we want to work with
% throughout the script
Dclass = dummyvar(Cars.cla);
Dclass = Dclass(:,2:end); %drop first class
Dclass = array2table(Dclass);
Cars = [Cars,Dclass ];
qvars = {'logpr', 'li', 'hp', 'wi', 'le', 'he', 'pl', 'home', 'Dclass1','Dclass2','Dclass3','Dclass4'}; 

%qvars = { 'li', 'cy'}; 
% Let's have a look at these variables for the 10 first rows
Cars(1:10, [qvars, 'type']);
%% Market share definition

% Compute sum of qu (quantity of cars sold) within (ma, ye)-pairs 
G = groupsummary(Cars, {'ma', 'ye'}, 'sum', 'qu'); % This creates an N*T Table with two variables, sum_su (sum of qu within groups) and 'GroupCount', which is the number of rows for that unit, which equals the size of the choiceset. 
G.Properties.VariableNames{'GroupCount'} = 'J_it'; 

% merge the new variables back to the Cars table 
Cars = join(Cars, G); 

% market share: quantity sold in percent of all cars sold in that (market,
% year). 
Cars.market_share = Cars.qu ./ Cars.sum_qu; 

%% Create x0, s0, etc. variables: 
% characteristics for the first car within a (ma, ye)-pair. 
% This gives us a "comparison" choice, j=0, for use in logit. 
% we simply pick the first car within each (ma, ye)

% the variables to take the first for 
these_vars = ['market_share', qvars, 'co', 'class', 'firm_country']; 

% find the first within each (ye, ma) 
first = @(x) x(1, :); 
C2 = sortrows(Cars, {'ma', 'ye', 'qu'});  % then the first has lowest sales 
G = groupsummary(C2, {'ma', 'ye'}, first, these_vars); % extract attributes for lowest selling car

% rename the variables to something sensible 
for i=1:numel(these_vars) 
    v = these_vars{i}; 
    oldname = sprintf('fun1_%s', v); 
    newname = sprintf('first_%s', v); 
    G.Properties.VariableNames{oldname} = newname; 
end

% delete the GroupCount variable: we already have it from previously 
G.GroupCount = []; 

% merge back to original dataset 
Cars = join(Cars, G); 

%% OLS for quantity

Cars.logqu = log(Cars.qu); 
oo = ones(size(Cars,1), 1);

[b,se] = ols.estimate(Cars.logqu, [oo, Cars.logpr, Cars.home])
print_par(b,se,{'const', 'logpr', 'home', 'sigma'}, 'OLS')

D1 = dummyvar(Cars.i); 
D2 = dummyvar(Cars.t); 
D3 = dummyvar(Cars.j); 
D4 = dummyvar(Cars.brand); 
D5 = dummyvar(Cars.frm); 

% ... FILL IN ... 
%% QUESTION 3
%% CONDITIONAL LOGIT @ OLS 
% Requires us to construct differences between each given car, j, and some
% baseline that was also available in the same choiceset (market i at time
% t). 
Cars.logdiff_s = log(Cars.market_share) - log(Cars.first_market_share); 

Cars.d_logpr = Cars.logpr - Cars.first_logpr; 
Cars.d_li = Cars.li - Cars.first_li;
Cars.d_hp = Cars.hp - Cars.first_hp;
Cars.d_wi = Cars.wi - Cars.first_wi;
Cars.d_le = Cars.le - Cars.first_le;
Cars.d_he = Cars.he - Cars.first_he;
Cars.d_pl = Cars.pl - Cars.first_pl;
Cars.d_home = Cars.home - Cars.first_home; 
Cars.d_Dclass1 = Cars.Dclass1 - Cars.first_Dclass1; 
Cars.d_Dclass2 = Cars.Dclass2 - Cars.first_Dclass2; 
Cars.d_Dclass3 = Cars.Dclass3 - Cars.first_Dclass3; 
Cars.d_Dclass4 = Cars.Dclass4 - Cars.first_Dclass4; 

%inlcude raw dummies (with notmalization) 
x_l = [ Cars.d_logpr, Cars.d_li, Cars.d_hp, Cars.d_wi,Cars.d_le, Cars.d_he, Cars.d_pl,Cars.d_home,Cars.d_Dclass1,Cars.d_Dclass2,Cars.d_Dclass3,Cars.d_Dclass4];

% x matrix without any dummies
%x_l = [ Cars.d_logpr, Cars.d_li, Cars.d_hp, Cars.d_wi,Cars.d_le, Cars.d_he, Cars.d_pl,Cars.d_home];

y_l = Cars.logdiff_s;
[b_l,se] = ols.estimate(y_l,x_l);
%[b_l,se] = ols.estimate(Cars.logdiff_s, [ Cars.d_logpr, Cars.d_home,Cars.d_we,Cars.d_li,Cars.d_hp,Cars.d_wi,Cars.d_le,Cars.d_he,Cars.d_cy])

%print_par(b_l,se,{'logpr', 'liters', 'horsepower','width','length','height','places','home', 'sigma'}, 'Logdiff')
print_par(b_l,se,{'logpr', 'liters', 'horsepower','width','length','height','places','home','class2','class3','class4','class5', 'sigma'}, 'Logdiff')

% perhaps this code can become useful 
% for i=1:numel(qvars) 
%     v = qvars{i}; 
%     firstname = sprintf('first_%s', v); 
%     Cars.(varname) = Cars.(v) - Cars.(firstname); 
% end

k = numel(x_l(1,:));
se_l = se;
s2_l = b_l(end)^2;              
res_l = y_l-x_l*b_l(1:k);   
ybar = y_l-mean(y_l);        
R2_l = 1-(res_l'*res_l)/(ybar'*ybar);  

%% From unbalanced to balanced (with missings) 
% Cars is balanced in i,t but not in j (since not all cars were available
% in all years/markets). In the logit model, it is simpler to work with a
% balanced dataset. 

tt = Cars.t; 
ii = Cars.i; 
jj = Cars.j; 

N = numel(unique_ma); 
J = numel(unique_co); 
T = numel(unique_ye); 
K_q = numel(qvars); 

% convenient conversion vector: converts from the unbalanced Cars table 
% format to a fully N*T*J-balanced array 
I = sub2ind([N,T,J],ii,tt,jj);

% initialize scalar variables
brand               = nan(N*T*J,1); 
firm                = nan(N*T*J,1); 
class               = nan(N*T*J,1); 
model               = nan(N*T*J,1); 
type                = nan(N*T*J,1); 
q                   = nan(N*T*J,1); 
market_share        = nan(N*T*J,1); 
logdiff_s           = nan(N*T*J,1); 
first_market_share  = nan(N*T*J,1); 

% convert 
q(I,1)                  = Cars.qu; 
brand(I,1)              = Cars.brand; 
firm(I,1)               = Cars.firm; 
class(I,1)              = Cars.class; 
model(I,1)              = Cars.model; 
type(I,1)               = Cars.type; 
market_share(I,1)       = Cars.market_share; 
logdiff_s(I,1)          = Cars.logdiff_s; 
first_market_share(I,1) = Cars.first_market_share; 

% extract a whole vector of variables at once 
X          = nan(N*T*J,K_q); 
X(I,:)     = Cars(:, qvars).Variables; 

% useful: extracting the corresponding names for categorical variables 
firm_names  = categories(Cars.firm); 
brand_names = categories(Cars.brand); 
model_names = categories(Cars.model); 
type_names  = categories(Cars.type); 

testX = reshape(X,N*T,J,K_q) % in exercise 8 it was NxJxK in this case N = N*T

 

%% Market shares 

q = reshape(q,N,T,J); 

tot_sales = nansum(q, 3); % sum over 3rd dimension (= cars); nansum skips over NaN (missing) cells 
s1 = q ./ tot_sales; % inside market share 
[s1_min, j_min] = min(s1, [], 3); %j_min will be N*T

% Find car attributes for the car with the lowest market share
X_min    = nan(N, T, 1, K_q); 
X = reshape(X, N, T, J, K_q); 
firm_min = nan(N,T); 
class_min = nan(N,T); 
brand_min = nan(N,T); 
for i=1:N
    for t=1:T
        j = j_min(i,t); 
        X_min(i,t,1,:) = X(i,t,j_min(i,t),:); 
        I_ = sub2ind([N,T,J],i,t,j_min(i,t)); 
        firm_min(i,t) = firm(I_); 
        class_min(i,t) = class(I_); 
        brand_min(i,t) = brand(I_); 
    end
end




%% QUESTION 1
%Select Variables
D1 = dummyvar(Cars.i); 
D2 = dummyvar(Cars.t); 
D3 = dummyvar(Cars.j); 
D4 = dummyvar(Cars.brand); 
D5 = dummyvar(Cars.cla); 

y_ols = log(Cars.qu);
const = ones(size(Cars,1), 1);
x_ols = [const,Cars.logpr,Cars.home];
zvars = [Cars.pop/1000000,Cars.tax,Cars.IncShare_p50p90,Cars.IncShare_p90p100,Cars.avexr./100,Cars.avdexr./100,Cars.engdp/10^12];
% we do not use doors, acceleration and cylinder due to NaNs and time
% invariance - weight not significant and thus dropped
qvars = [Cars.li,Cars.hp,Cars.wi,Cars.le/100,Cars.he/100,Cars.pl,];

xzq_ols = [const,Cars.logpr,Cars.home,zvars,qvars];

%POLS
%(1) WITHOUT CONTROLS AND DUMMIES
n= size(y_ols,1);
k = numel(x_ols(1,:));

lbly = {'SALES'};
lblx = {'CONSTANT'; 'PRICE';'HOME'};

[b,se] = ols.estimate(y_ols,x_ols)

b_ols1 = b;
se_ols1 = se;
s2_ols1 = b_ols1(end)^2;              
res_ols1 = y_ols-x_ols*b_ols1(1:k);   
ybar1 = y_ols-mean(y_ols);        
R2_ols1 = 1-(res_ols1'*res_ols1)/(ybar1'*ybar1);  


disp('OLS estimates (1) - NO control/FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols1(r,1),se_ols1(r,1),b_ols1(r,1)./se_ols1(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols1)]);
disp(['Sigma         ', num2str(sqrt(s2_ols1))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');


%POLS with FE
%(2) COUNTRY DUMMIES
xd1_ols2 = [x_ols,D1(:,2:5)];

n= size(y_ols,1);
k = numel(xd1_ols2(1,:));

[b,se] = ols.estimate(y_ols,xd1_ols2);

b_ols2 = b;
se_ols2 = se;
s2_ols2 = b_ols2(end)^2;                            
res_ols2 = y_ols-xd1_ols2*b_ols2(1:k);                
ybar2 = y_ols-mean(y_ols);       
R2_ols2 = 1-(res_ols2'*res_ols2)/(ybar2'*ybar2);  

disp('OLS estimates (2) - Country FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols2(r,1),se_ols2(r,1),b_ols2(r,1)./se_ols2(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols2)]);
disp(['Sigma         ', num2str(sqrt(s2_ols2))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');
%POLS
%(3) YEAR DUMMIES
xd2_ols3 = [x_ols,D2(:,2:30)];

n= size(y_ols,1);
k = numel(xd2_ols3(1,:));

[b,se] = ols.estimate(y_ols,xd2_ols3);

b_ols3 = b;
se_ols3 = se;
s2_ols3 = b_ols3(end)^2;                            
res_ols3 = y_ols-xd2_ols3*b_ols3(1:k);                
ybar3 = y_ols-mean(y_ols);       
R2_ols3 = 1-(res_ols3'*res_ols3)/(ybar3'*ybar3); 

% Print results
disp('OLS estimates (3) - Year FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols3(r,1),se_ols3(r,1),b_ols3(r,1)./se_ols3(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols3)]);
disp(['Sigma         ', num2str(sqrt(s2_ols3))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');

%POLS
%(4) CAR DUMMIES
xd3_ols5 = [x_ols,D3(:,2:334)];

n= size(y_ols,1);
k = numel(xd3_ols5(1,:));

[b,se] = ols.estimate(y_ols,xd3_ols5);

b_ols5 = b;
se_ols5 = se;
s2_ols5 = b_ols5(end)^2;                            
res_ols5 = y_ols-xd3_ols5*b_ols5(1:k);                
ybar5 = y_ols-mean(y_ols);       
R2_ols5 = 1-(res_ols5'*res_ols5)/(ybar5'*ybar5);  

% Print results
disp('OLS estimates (4) - Car FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols5(r,1),se_ols5(r,1),b_ols5(r,1)./se_ols5(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols5)]);
disp(['Sigma         ', num2str(sqrt(s2_ols5))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');

%POLS
%(6) COUNTRY,YEAR & CAR DUMMIES
xd1d2d3_ols6 = [x_ols,D1(:,2:5),D2(:,2:30),D3(:,2:334)];

n= size(y_ols,1);
k = numel(xd1d2d3_ols6(1,:));

[b,se] = ols.estimate(y_ols,xd1d2d3_ols6);

b_ols6 = b;
se_ols6 = se;
s2_ols6 = b_ols6(end)^2;                            
res_ols6 = y_ols-xd1d2d3_ols6*b_ols6(1:k);                
ybar6 = y_ols-mean(y_ols);       
R2_ols6 = 1-(res_ols6'*res_ols6)/(ybar6'*ybar6);  

% Print results
disp('OLS estimates (5) - Country,Year & Car FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols6(r,1),se_ols6(r,1),b_ols6(r,1)./se_ols6(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols6)]);
disp(['Sigma         ', num2str(sqrt(s2_ols6))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');
%POLS
%(5) Segment DUMMIES
xd5_ols7 = [x_ols,D5(:,2:end)];

n= size(y_ols,1);
k = numel(xd5_ols7(1,:));

[b,se] = ols.estimate(y_ols,xd5_ols7);

b_ols7 = b;
se_ols7 = se;
s2_ols7 = b_ols7(end)^2;                            
res_ols7 = y_ols-xd5_ols7*b_ols7(1:k);                
ybar7 = y_ols-mean(y_ols);       
R2_ols7 = 1-(res_ols7'*res_ols7)/(ybar7'*ybar7);  

latex_table = latex(sym(b_ols7));

disp('OLS estimates (6) - Segment FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols7(r,1),se_ols7(r,1),b_ols7(r,1)./se_ols7(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols7)]);
disp(['Sigma         ', num2str(sqrt(s2_ols7))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');

%POLS
%(7) COUNTRY,YEAR & Segment DUMMIES
xd1d2d5_ols8 = [x_ols,D1(:,2:5),D2(:,2:30),D5(:,2:end)];

n= size(y_ols,1);
k = numel(xd1d2d5_ols8(1,:));

[b,se] = ols.estimate(y_ols,xd1d2d5_ols8);

b_ols8 = b;
se_ols8 = se;
s2_ols8 = b_ols8(end)^2;                            
res_ols8 = y_ols-xd1d2d5_ols8*b_ols8(1:k);                
ybar8 = y_ols-mean(y_ols);       
R2_ols8 = 1-(res_ols8'*res_ols8)/(ybar8'*ybar8); 

disp('OLS estimates (7) - Country, Year & Segment FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols8(r,1),se_ols8(r,1),b_ols8(r,1)./se_ols8(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols8)]);
disp(['Sigma         ', num2str(sqrt(s2_ols8))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');


%POLS
%(8) ONLY CONTROLS
y_ols = log(Cars.qu);
const = ones(size(Cars,1), 1);
x_ols = [const,Cars.logpr,Cars.home];

xzq_ols = [const,Cars.logpr,Cars.home,zvars,qvars];

data = [y_ols,xzq_ols];
data(any(isnan(data), 2), :) = [];
 
 y_ols9  =  data(:,1);
 xzq_ols9 = [data(:,2:end)];
 
 n= size(y_ols9,1);
 k = numel(xzq_ols9(1,:));
 
 
 lbly = {'SALES'};
 
lblx = {'Constant';'Price';'Home';'pop';'tax';'IncShare_p50p90';'IncShare_p90p100';'avexr';'avdexr';'engdp';'li';'hp';'wi';'le';'he';'pl'};

 [b,se] = ols.estimate(y_ols9,xzq_ols9)
 
 b_ols9 = b;
 se_ols9 = se;
 s2_ols9 = b_ols9(end)^2;              
 res_ols9 = y_ols9-xzq_ols9*b_ols9(1:k);    
 ybar9 = y_ols9-mean(y_ols9);        
 R2_ols9 = 1-(res_ols9'*res_ols9)/(ybar9'*ybar9);  
 
  disp('OLS estimates (8) - Controls');
  disp('================================================');
  disp('             Param. est.  s.e.         t-value');
  for r=1:length(lblx)
    fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols9(r,1),se_ols9(r,1),b_ols9(r,1)./se_ols9(r,1)]);
  end
  disp('------------------------------------------------');
  disp(['R-squared     ', num2str(R2_ols9)]);
  disp(['Sigma         ', num2str(sqrt(s2_ols9))]);
  disp(['Observations         ', num2str(n)]);
  disp('================================================');
  
%POLS
%(9) CONTROLS + COUNTRY,YEAR & CAR DUMMIES
xzqd1d2d3_ols10 = [const,Cars.logpr,Cars.home,zvars,qvars,D1(:,2:5),D2(:,2:30),D3(:,2:334)];

data = [y_ols,xzqd1d2d3_ols10];
data(any(isnan(data), 2), :) = [];

y_ols10  =  data(:,1);
xzqd1d2d3_ols10 = [data(:,2:end)];

n= size(y_ols10,1);
k = numel(xzqd1d2d3_ols10(1,:));

[b,se] = ols.estimate(y_ols10,xzqd1d2d3_ols10)
 
b_ols10 = b;
se_ols10 = se;
s2_ols10 = b_ols10(end)^2;              
res_ols10 = y_ols10-xzqd1d2d3_ols10*b_ols10(1:k);    
ybar10 = y_ols10-mean(y_ols10);        
R2_ols10 = 1-(res_ols10'*res_ols10)/(ybar10'*ybar10);  
 
disp('OLS estimates (9) - Controls and all FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
 for r=1:length(lblx)
   fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
         lblx{r,1}, [b_ols10(r,1),se_ols10(r,1),b_ols10(r,1)./se_ols10(r,1)]);
 end
 disp('------------------------------------------------');
 disp(['R-squared     ', num2str(R2_ols10)]);
 disp(['Sigma         ', num2str(sqrt(s2_ols10))]);
 disp(['Observations         ', num2str(n)]);
 disp('================================================');

%POLS
%(10) CONTROLS + COUNTRY,YEAR & Segment DUMMIES
xzqd1d2d5_ols11 = [const,Cars.logpr,Cars.home,zvars,qvars,D1(:,2:5),D2(:,2:30),D5(:,2:end)];

data = [y_ols,xzqd1d2d5_ols11];
data(any(isnan(data), 2), :) = [];

y_ols11  =  data(:,1);
xzqd1d2d5_ols11 = [data(:,2:end)];

n= size(y_ols11,1);
k = numel(xzqd1d2d5_ols11(1,:));

[b,se] = ols.estimate(y_ols11,xzqd1d2d5_ols11)
 
b_ols11 = b;
se_ols11 = se;
s2_ols11 = b_ols11(end)^2;              
res_ols11 = y_ols11-xzqd1d2d5_ols11*b_ols11(1:k);    
ybar11 = y_ols11-mean(y_ols11);        
R2_ols11 = 1-(res_ols11'*res_ols11)/(ybar11'*ybar11);  
 
disp('OLS estimates (10) - Controls and FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
 for r=1:length(lblx)
   fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
         lblx{r,1}, [b_ols11(r,1),se_ols11(r,1),b_ols11(r,1)./se_ols11(r,1)]);
 end
 disp('------------------------------------------------');
 disp(['R-squared     ', num2str(R2_ols11)]);
 disp(['Sigma         ', num2str(sqrt(s2_ols11))]);
 disp(['Observations         ', num2str(n)]);
 disp('================================================');
%% HERE STARTS MALTE'S PART
%% QUESTION 4
%% %% Estimation preparation
%y = logdiff_s(I,1); % 50100x1 vector
%x = X; % NxTxJ 3 dimensional matrix

%SELECT THETA0
%--------------
s1 = reshape(s1,N*T*J,1);
theta0 = clogit_w8.starting_values(s1,testX); 
%theta_bl = b_l(1:end-1);
%theta1 = ones(10,1);
%theta1 = theta1*0.00001;

q = @(theta) clogit_w8.criterion(s1,testX,theta); %update to chosen theta0

% ----------------------------------------------------------------
% Optimizer settings
opt = optimoptions(@fminunc,...
    'Display','iter',...                 % 'iter' (for every iteration), 'final' or 'off'.
    'Algorithm','quasi-newton',...       % Alternatively, we can use 'trust-region'
    'StepTolerance',1e-08,...            % Optimizer will terminate if the norm of (thetaNew - thetaOld) is smaller than this
    'FiniteDifferenceType','central',... % These are more precise than 'forward', but require twice as many computations
    'MaxFunEvals',5000,...
    'HessUpdate','bfgs'...               % 'iter' (for every iteration), 'final' or 'off'.
    );

opt_search =  optimset('Display','iter',...
    'TolX',1e-08,...
    'TolFun',1e-08,...
    'MaxIter',5000,...
    'MaxFunEvals',50000,'HessUpdate','bfgs');

%% Estimation using BFGS (fminunc) 
N=150;
[thetahat,se_thetahat,out] = estimation_clog.estimate_m(q,theta0,N,opt,'Hessian');%fminunc
%[thetahat,se,out] = estimation_clog_search.estimate_m(q,theta0,N,opt_search,'Outer Product');%fminsearch
%[thetahat,se,out,hess] =estimation_clog.estimate_m(q,theta_bl,opt,'Outer Product');%fminsearch

g = estimation_clog.centered_grad(q,thetahat); 
fprintf('Norm of gradient: %g.\n',norm(sum(g,1))); 

n= size(s1,1);
k = numel(testX(1,1,:));
X_r = reshape(testX,N*J,k);
res_opt = s1-X_r*thetahat;    
yopt = s1-mean(s1);        
R2_opt = 1-(res_opt'*res_opt)/(yopt'*yopt); 
%% Print
%for i=1:numel(thetahat);
%    fprintf('%24s: %10f (t = %6.2f) \n',labels{i},thetahat(i),thetahat(i)/se(i));
%end; 


%% %%QUESTION 5
%% Now: compute the elasticity of demand wrt. "price/log(income)" 
k_price = 1; 

% 1. baseline probabilities 
ccp = clogit_w8.choice_probabilities(testX,thetahat); 

% 2. coutnerfactual probabilities 
E_own = nan(N,J); 
E_cross = nan(N,J); 
dpdx = nan(N,J); 
ii=1:J; 
for j=1:J
    %keyboard
    x2 = testX; 
    x2(:,j,k_price) = x2(:,j,k_price)*1.00001; 
    ccp2 = clogit_w8.choice_probabilities(x2,thetahat); 
    dccp = ccp2-ccp; 
    dx = x2(:,j,k_price)-testX(:,j,k_price); 
    x0 = testX(:,j,k_price); 
    
    % own-price elasticity 
    dy = dccp(:,j); 
    y0 = ccp(:,j); 
    
    E_own(:,j) = dy./y0 .* x0./dx; 
    
    % cross-price elasticity
    dy = dccp(:,ii(ii~=j));  % delete jth column => N*(J-1)
    y0 = ccp(:,ii(ii~=j));   % N*(J-1)
    E_cross(:,j) = mean(dy./y0,2,'omitnan').* x0./dx; 
end
fprintf('Avg. own-price elasticities:   %.2f%%\n',mean(E_own(:),'omitnan')); 
fprintf('Avg. cross-price elasticities: %.2f%%\n',mean(E_cross(:),'omitnan')); 


%% Comparing own-price elasticities 
% Comparing the avg. elasticity for cars with x_k =0 vs. =1. 
k = 8; 
idx1 = testX(:,:,k)==1; 
idx0 = testX(:,:,k)==0; 
%fprintf('Mean own elast.: %f (with %12s=1)\n',mean(E_own(idx1),'omitnan'),labels{k}); 
%fprintf('Mean own elast.: %f (with %12s=0)\n',mean(E_own(idx0),'omitnan'),labels{k}); 
disp('Mean own elast. home = 1:');
disp(mean(E_own(idx1)));

disp('Mean own elast. home = 0:');
disp(mean(E_own(idx0)));

disp('Mean cross elast. home = 1:');
disp(mean(E_cross(idx1)));

disp('Mean cross elast. home = 0:');
disp(mean(E_cross(idx0)));


%% 
% to 5. heterogenious price effects
Cars.logshare = log(Cars.market_share); 
Interaction = Cars.logpr.*Cars.home
oo = ones(size(Cars,1), 1);

[b_heterog,se] = ols.estimate(Cars.logshare, [oo, Cars.logpr, Interaction])
print_par(b_heterog,se,{'const','logpr','Price*Home','sigma' }, 'Heterogenious OLS')

% repeat with dummies

[b_heterog2,se] = ols.estimate(Cars.logshare, [oo, Cars.logpr, Interaction, Cars.Dclass1, Cars.Dclass2, Cars.Dclass3,Cars.Dclass4 ])
print_par(b_heterog2,se,{'const','logpr','Price*Home', 'class2','class3','class4','class5','sigma' }, 'Heterogenious OLS')


lbly = {'PRICE'};
lblx = {'const','logpr','Price*Home', 'class2','class3','class4','class5'};
y_ols = Cars.logshare;
x_ols = [oo, Cars.logpr, Interaction];
n= size(y_ols,1);
k = numel(x_ols(1,:));

b_ols = b_heterog2;
se_ols = se;
s2_ols = b_ols(end)^2;                            
res_ols = y_ols-x_ols*b_ols(1:k);                
ybar = y_ols-mean(y_ols);       
R2_ols = 1-(res_ols'*res_ols)/(ybar'*ybar);  

disp(['R-squared     ', num2str(R2_ols)]);

%% QUESTION 6
%%POLS ONLY CONTROLS price effects on home
%Cars.logshare = log(Cars.market_share); 
oo = ones(size(Cars,1), 1);

[b_home,se_home] = ols.estimate(Cars.logpr, [oo, Cars.home])

print_par(b_home,se_home,{'const','home','sigma' }, 'Price on home OLS')

y_ols = Cars.logpr;
x_ols = [oo,Cars.home];
n= size(y_ols,1);
k = numel(x_ols(1,:));

lbly = {'PRICE'};
lblx = {'CONSTANT';'HOME'};

b_ols = b_home;
se_ols = se_home;
s2_ols = b_ols(end)^2;                            
res_ols = y_ols-x_ols*b_ols(1:k);                
ybar = y_ols-mean(y_ols);       
R2_ols = 1-(res_ols'*res_ols)/(ybar'*ybar);  

% Print results
disp('OLS estimates (9) - price on home');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols(r,1),se_ols(r,1),b_ols(r,1)./se_ols(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols)]);
disp(['Sigma         ', num2str(sqrt(s2_ols))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');

%% 

%%POLS ONLY CONTROLS AND CAR FIXED EFFECTS
%Cars.logshare = log(Cars.market_share); 
oo = ones(size(Cars,1), 1);

[b_homefe,se_homefe] = ols.estimate(Cars.logpr, [oo, Cars.home, D2])

y_ols_ = Cars.logpr;
x_ols_ = [oo,Cars.home, D2];
n= size(y_ols_,1);
k = numel(x_ols_(1,:));

lbly = {'PRICE'};
lblx = {'CONSTANT';'HOME';};

b_ols_ = b_homefe;
se_ols_ = se_homefe;
s2_ols_ = b_ols_(end)^2;                            
res_ols_ = y_ols_-x_ols_*b_ols_(1:k);                
ybar_ = y_ols_-mean(y_ols_);       
R2_ols_ = 1-(res_ols_'*res_ols_)/(ybar_'*ybar_);  

% Print results
disp('OLS estimates (10) - price on home + type FE');
disp('================================================');
disp('             Param. est.  s.e.         t-value');
for r=1:length(lblx)
  fprintf('%-10s %11.4f %11.4f %11.4f\n', ...
          lblx{r,1}, [b_ols_(r,1),se_ols_(r,1),b_ols_(r,1)./se_ols_(r,1)]);
end
disp('------------------------------------------------');
disp(['R-squared     ', num2str(R2_ols_)]);
disp(['Sigma         ', num2str(sqrt(s2_ols_))]);
disp(['Observations         ', num2str(n)]);
disp('================================================');

%% Test for graph

p_nohome = nan(J,1);
p_home = nan(J,1);
ii=1:J; 
for j=1:J
    p_home(j,1) = mean(Cars.logpr(idx1(:,j)));
    p_nohome(j,1) = mean(Cars.logpr(idx0(:,j)),'omitnan');
end
%fprintf('Avg. own-price elasticities:   %.2f%%\n',mean(p_nohome(:),'omitnan')); 

%hist(p_home,100, 'facecolor',  'blue', 'facealpha',.5,'EdgeColor','none');
%hold on
%hist(p_nohome,100, 'facecolor',  'red', 'facealpha',.5, 'EdgeColor','none');
%title('Comparision between home and away prices per car type j');
%legend('home market price distr.','foreign market price distr.');
%xlabel('Car type j');
%ylabel('log(price)'); 


x = p_home;
y = p_nohome;
h1 = histogram(x,100);
hold on
h2 = histogram(y,100);
xlabel('Log price');
legend('home market price distr.','foreign market price distr.');


%% graph 2
home = Cars.home;
EOWN = reshape(E_own, 150*334,1);
EOWN = EOWN(~isnan(EOWN));
el_nohome = [EOWN, home];
el_home = [EOWN, home];
el_nohome(el_nohome(:, 2)== 1, :)= [];
el_home(el_home(:, 2)== 0, :)= [];

%fprintf('Avg. own-price elasticities:   %.2f%%\n',mean(p_nohome(:),'omitnan')); 
x = el_home(:,1);
y = el_nohome(:,1);
h1 = histogram(x,100, 'normalization','probability');
hold on
h2 = histogram(y,100, 'normalization','probability');
xlabel('Own price elast.');
ylabel('Probability');
legend('home market own elast. distr.','foreign market own elast. distr.');


%% graph 3
price = Cars.logpr;
EOWN = reshape(E_own, 150*334,1);
EOWN = EOWN(~isnan(EOWN));
p_nohome = [price, home];
p_home = [price, home];
p_nohome(p_nohome(:, 2)== 1, :)= [];
p_home(p_home(:, 2)== 0, :)= [];

%fprintf('Avg. own-price elasticities:   %.2f%%\n',mean(p_nohome(:),'omitnan')); 
x = p_home(:,1);
y = p_nohome(:,1);
h1 = histogram(x,100, 'normalization','probability');
hold on
h2 = histogram(y,100, 'normalization','probability');
xlabel('Log price');
ylabel('Probability');

legend('home market price distr.','foreign market price distr.');
