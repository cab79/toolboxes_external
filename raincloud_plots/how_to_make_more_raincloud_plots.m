% download socio-economic status dataset from kaggle...
% https://mwww.kaggle.com/sdorius/countryses

% insert your path to the plotting scripts here
addpath('C:\Data\Matlab\raincloud_plots');
addpath('C:\Data\Matlab\cbrewer');

% insert your path to the data here
datadir     = '/Users/marshall/Desktop';
fn          = 'GLOBCSES.Final20170714.csv';
format_spec = '%f%C%C%f%f%C%f%f%C%C';

d           = readtable(fullfile(datadir, fn), 'Delimiter', ',' , ...
                'Format', format_spec);

%% get nice colours from color brewer
% (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

[cb] = cbrewer('qual', 'Set1', 10, 'pchip');

%% pull relevant variables out of table

% column vectors from the table columns
yr          = table2array(d(:,4));
ses         = table2array(d(:,5));
gdp_ppc     = table2array(d(:,7));
yrs_educ    = table2array(d(:,8));
region      = table2array(d(:,10));

% logicals to determine region
is_eur      = (region == 'East Europe' | region == 'West Europe' | ...
                region == 'South Europe' | region == 'North Europe' );

is_afr      = (region == 'East Africa' | region == 'West Africa' | ...
                region == 'South Africa' | region == 'North Africa' | ...
                region == 'Middle Africa');

is_amr      = (region == 'South America' | region == 'Central America' | ...
                region == 'Caribbean');

is_asi      = (region == 'East Asia' | region == 'West Asia' | ...
                region == 'South Asia' | region == 'Southeast Asia');

%% example of a two-d raincloud

yrs_to_select   = [1950 2010];
nyrs            = length(yrs_to_select);

data = [];
for i = 1:nyrs
    
    target_yr = yrs_to_select(i);
    
    amrdata = is_amr & yr == target_yr;
    
    data{i,1}(2,:) = ses(amrdata);
    data{i,1}(1,:) = log10(gdp_ppc(amrdata));
    
end

clz = cb([2 8],:);

% plot
figure('units', 'normalized', 'outerposition', [0 0 1 1])
[h, ax] = twod_raincloud(data, clz, 0);

h.s{1}.SizeData = 180;
h.s{2}.SizeData = 180;
h.s{3}.SizeData = 180;
h.s{4}.SizeData = 180;

xlabel(ax(1), 'log10 per-capita GDP ($)');
ylabel(ax(1), 'socio-economic status');
legend(ax(1), {'1950';'2010'}, 'Location', 'SouthEast');
set(ax(1), 'FontSize', 14);

%% example of n-rainclouds

clear h

yrs_to_select = [1950 1980 2010];
nyrs = length(yrs_to_select);

clz = cb([3 5],:);

figure('units', 'normalized', 'outerposition', [0.05 0.05 0.3 0.6])

for i = 1:nyrs
    
    target_yr = yrs_to_select(i);
    
    data = [];
    data{1} = yrs_educ(is_eur & yr == target_yr);
    data{2} = yrs_educ(is_asi & yr == target_yr);
    
    subplot(1, nyrs, i)
    h(i) = n_rainclouds(data, clz);
    set(gca, 'XLim', [0 14]);
    
    % flip plot so 'x' axis is vertical and 'y' is horizontal
    view([-90 -90]);
    set(gca, 'Xdir', 'reverse');
    
    % labels
    title(num2str(target_yr));
    if i == 1
        xlabel('education (years)');
    end
end

legend({'europe';'asia'}, 'Location', 'SouthEast');

%% another example of nrainclouds

yrs_to_select = [1970 1990 2010];
nyrs = length(yrs_to_select);
data = [];

for i = 1:nyrs
    target_yr = yrs_to_select(i);
    data{i} = yrs_educ(yr == target_yr);
end

clz = cbrewer('qual', 'Set1', 3);
figure('units', 'normalized', 'outerposition', [0.55 0.05 0.6 0.4])
n_rainclouds(data, clz);
set(gca, 'XLim', [-1 14]);
legend({'1970';'1990';'2010'});
xlabel('education (years)');