%% Read data 
% Reads the car dataset and merges on income inequality data 
% AUTHOR: Anders Munk-Nielsen, U. of Copenhagen 
% DATE: 2020-11-17
%
% SOURCES: Full credit to the amazing people that have collected the
% underlying datasets: 
% Frank Verboven's car dataset 
%   --> https://sites.google.com/site/frankverbo/data-and-software/data-set-on-the-european-car-market
% the World Inequality Database 
%   --> https://wid.world/

Cars = readtable('cars.xlsx'); 
Cars = convertvars(Cars, {'type', 'brand', 'model'}, 'categorical'); 


%% Drop observations of cars (j) that were sold in only one market-period 

if DROPSINGLEYEARCARS
    % count number of (market, year)-observations for each car (brand,
    % type)-pair. 
    G = groupcounts(Cars, {'brand', 'type'});
    G.Properties.VariableNames{'GroupCount'} = 'num_yearobs_this_car'; 
    
    % merge back 
    Cars = join(Cars, G); 
    I = Cars.num_yearobs_this_car < minyearobs; 
    fprintf('Dropping %d rows of cars sold in too few (market, country)s.\n', sum(I)); 
    Cars(I, :) = []; 
    
    % clean up
    Cars.num_yearobs_this_car = []; 
end

%% 

Cars = sortrows(Cars, {'ma', 'ye', 'org'}); 

Cars.we = Cars.we / 1000; % rescale to weight in tonnes (makes coefficients nicer to look at) 

market_name        = {'Belgium','France','Germany','Italy','UK'}; 
class_name         = {'subcompact','compact','intermediate','standard','luxury'}; 
firm_code_name     = {'AlfaRomeo','Audi','BMW','Citroen','Daihatsu','Ferrari','Fiat','Ford','Honda','Hyundai','Innocenti','Jaguar','Lada','Lancia','Mazda','Mercedes','Mitsubishi','NissanDatsun','OpelVauxhall','Peugeot','Porsche','Renault','RoverTriumph','Saab','Seat','Skoda','Subaru','Suzuki','Toyota','Volkswagen','Volvo','Yugo','Talbot','Kia','Daewoo','Rover','MCC','Lexus','DAF','Dino','Autobianchi','Triumph','Princess','TalbotHillmanChrysler','TalbotSimca','TalbotMatra   1   AlfaRomeo','Audi','BMW','Citroen','Daihatsu','Ferrari','Fiat','Ford','Honda','Hyundai','Innocenti','Jaguar','Lada','Lancia','Mazda','Mercedes','Mitsubishi','NissanDatsun','OpelVauxhall','Peugeot','Porsche','Renault','RoverTriumph','Saab','Seat','Skoda','Subaru','Suzuki','Toyota','Volkswagen','Volvo','Yugo','Talbot','Kia','Daewoo','Rover','MCC','Lexus','DAF','Dino','Autobianchi','Triumph','Princess','TalbotHillmanChrysler','TalbotSimca','TalbotMatra'}; 
origin_code_name   = {'France','Germany','Italy','JapanKorea','Spain','Sweden','UK','EasternEurope','US','Brazil'}; 
location_code_name = {'Belgium','CheckRepublic','France','Germany','Italy','Japan','Korea','Netherlands','Romania','Spain','Sweden','UK','Russia','Yugoslavia','Poland','US','Finland','Australia','Hungary','Portugal','India','Mexico'}; 

[unique_ye,~,Cars.t] = unique(Cars.ye); 
[unique_co,~,Cars.j] = unique(Cars.co); 
[unique_ma,~,Cars.i] = unique(Cars.ma); 

Cars.market = categorical(market_name(Cars.ma)'); 
Cars.class = categorical(class_name(Cars.cla)'); 
Cars.firm = categorical(firm_code_name(Cars.frm)'); 
Cars.firm_country = categorical(origin_code_name(Cars.org)'); 
Cars.firm_production = categorical(location_code_name(Cars.loc)'); 


%% Income shares
% Variables that measure income inequality from the World Inequality
% Database. 
% 
% Creates three variables added to Cars: 
%   IncShare_p0p50: share of national income going to percentiles 0 through 50
%   IncShare_p50p90: same, but for 50-90. 
%   IncShare_p90p100: same, but for the top decile. 

assert(isfile('income_shares.csv'), 'Unable to find required file, "income_shares.csv". Sure you put it in the folder I am currently standing in?')

is = readtable('income_shares.csv'); 
is = convertvars(is, 'Percentile', 'categorical'); 

% fill missings: variables are country-specific 
is = grouptransform(is, {'Percentile'}, 'linearfill'); 

is = stack(is, {'Belgium', 'France', 'Germany', 'Italy', 'UK'}, ...
    'NewDataVariableName', 'IncomeShare', 'IndexVariableName', 'market'); 

% create year variable with same form as that in Cars (YY, not YYYY)
is.ye = is.Year - 1900; 
is.Year = []; 

is = sortrows(is, {'market', 'Percentile', 'ye'});

% tall->wide
is = unstack(is, 'IncomeShare', 'Percentile'); 

% rename variables 
is.Properties.VariableNames{'p0p50'}   = 'IncShare_p0p50'; 
is.Properties.VariableNames{'p50p90'}  = 'IncShare_p50p90'; 
is.Properties.VariableNames{'p90p100'} = 'IncShare_p90p100'; 

% merge back 
Cars = join(Cars, is) ; 

%% National income
% Creates the variable, NationalInc, which is avg. national income per
% capita in the country in that year. 

assert(isfile('income.csv'), 'Unable to find required file, "income.csv". Sure you put it in the folder I am currently standing in?')

inc = readtable('income.csv'); 
inc = convertvars(inc, 'Percentile', 'categorical'); 
assert(not(any(ismissing(inc), 'all')))

inc = stack(inc, {'Belgium', 'France', 'Germany', 'Italy', 'UK'}, ...
    'NewDataVariableName', 'NationalInc', 'IndexVariableName', 'market'); 

% the year variable must conform to "ye" in Cars, which is YY. 
inc.ye = inc.Year - 1900; % create the new 
inc.Year = [];            % delete the old 

inc = sortrows(inc, {'market', 'Percentile', 'ye'});
inc.Percentile = []; % this variable is superfluous, as there is only one group

% merge back 
Cars = join(Cars, inc); 

% remove superfluous values from categorical variables: otherwise, we risk
% creating dummies for cars that have been deleted from the dataset. 
cat_vars = {'brand', 'type', 'model', 'firm', 'market', 'firm_country', 'firm_production'}; 
for i=1:numel(cat_vars)
    v = cat_vars{i}; 
    Cars.(v) = removecats(Cars.(v)); 
end
