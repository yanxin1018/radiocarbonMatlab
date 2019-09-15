function C14Time(Year, Err, h)
data = xlsread('IntCal13.xlsx');
H = min(abs(diff(data(:,1)))); % step size of the years in the given data. It is 5
if nargin < 3
    h = 0.1; % step size of the years
end

CalendarAge = (data(1,1):h:data(end,1))';
Interpolation_Fun = @pchip; % other options: spline, interp1
RadiocarbonAge = Interpolation_Fun(data(:,1), data(:,2), CalendarAge);
RadiocarbonErr = Interpolation_Fun(data(:,1), data(:,3), CalendarAge);

% find the radiocarbon err corresponding to the given "Year"
idx = find(RadiocarbonAge > Year - H/2 & RadiocarbonAge < Year + H/2);
err_rd = sum(RadiocarbonErr(idx)) / numel(idx);
sigma = sqrt(Err^2 + err_rd^2);

% select sufficient large range
k = 8;
range =find(RadiocarbonAge <=Year+k*Err & RadiocarbonAge>=Year-k*Err);
CalIn3Range=CalendarAge(range);
% locate on the y axis
RCIn3Range=RadiocarbonAge(range);
% calculate the probality densities 3 sigma range
PDensity3Sigma=exp(-0.5*(RCIn3Range-Year).^2/sigma^2);
% calculate the probabilities for 3 sigma range
Probability3Sigma=PDensity3Sigma/sum(PDensity3Sigma);
% sort from the highest to lowest probability
[probability, index]=sort(Probability3Sigma,'Descend');
% locate the 1 sigma range and 2 sigma range by summing up probabilities
P = 0; i=0;
while P<=0.682
    i=i+1;
    P=P+probability(i);
end
Sigma1thresh=i;
while P<=0.954
    i=i+1;
    P=P+probability(i);
end
Sigma2thresh=i;
Sigma1Range=range(index(1:Sigma1thresh));
Sigma2Range=range(index(1:Sigma2thresh));
% locate the calendar ranges for 1 sigma and 2 sigma
Cal1Sigma=CalendarAge(Sigma1Range);
Cal2Sigma=CalendarAge(Sigma2Range);
% find the probabilities for each calendar age
Probability1Sigma=Probability3Sigma(index(1:Sigma1thresh));
Probability2Sigma=Probability3Sigma(index(1:Sigma2thresh));
% plot the 1, 2, 3 sigma range
figure;
idx = Probability3Sigma > max(Probability3Sigma(:)) * 1e-5;
year_left = floor(min(CalIn3Range(idx))/50)*50;
year_right = ceil(max(CalIn3Range(idx))/50)*50;
idx2 = CalendarAge >= year_left & CalendarAge <= year_right;
cal_age = CalendarAge(idx2);
age_low = RadiocarbonAge(idx2) - RadiocarbonErr(idx2);
age_up = RadiocarbonAge(idx2) + RadiocarbonErr(idx2);

hold on
[AX H1 H2] = plotyy( cal_age, age_up, CalIn3Range(idx), Probability3Sigma(idx),  'plot','bar' ,'c');
hold(AX(2));
%bar(AX(2),CalIn3Range(idx),Probability3Sigma(idx), 'c','EdgeColor','c');
bar(AX(2), Cal2Sigma,Probability2Sigma, 'c','EdgeColor','c');
bar(AX(2), Cal1Sigma,Probability1Sigma, 'm','EdgeColor','m');
legend(AX(2),'99.7% ', '95.4%','68.2%');
fill([cal_age; flipud(cal_age)], [age_up; flipud(age_low)], 'g');
ylabel(AX(1),'Radiocarbon Age (BP)');
xlabel('Calendar Age (Cal BC/AD)');
set(AX(2),'YTick',[]);
set(AX(2),'YTickLabel','');
ylim(AX(1),[2*min(age_low)-max(age_low),1.1*max(age_up)-0.1*min(age_up)]);
ylim(AX(2),[-0.5*max(Probability3Sigma(idx)), 2*max(Probability3Sigma(idx))]);
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cal1, idx1] = sort(Cal1Sigma, 'ascend');
prob1 = Probability1Sigma(idx1);
idx_split = find(diff(cal1) > 1*H);
if numel(idx_split) > 0
    section1 = zeros(numel(idx_split)+1, 3);
    section1(1, :) = [cal1(1), cal1(idx_split(1)), sum(prob1(1:idx_split(1)))];
    for k = 1:numel(idx_split)-1
        section1(k+1, :) = [cal1(idx_split(k)+1), cal1(idx_split(k+1)),...
            sum(prob1(idx_split(k)+1:idx_split(k+1)))];
    end
    section1(end, :) = [cal1(idx_split(end)+1), cal1(end),...
        sum(prob1(idx_split(end)+1:end))];
else
    section1 = [cal1(1), cal1(end), sum(prob1)];
end
fprintf('68.2%% Probability:\n');
for k = 1:numel(idx_split)+1
    fprintf('%s (%.2f%%) %s \n', convertAge(floor(section1(k, 1))),...
        100*section1(k, 3), convertAge(ceil(section1(k, 2))));
end

[cal2, idx2] = sort(Cal2Sigma, 'ascend');
prob2 = Probability2Sigma(idx2);
idx_split = find(diff(cal2) > 1*H);
if numel(idx_split) > 0
    section2 = zeros(numel(idx_split)+1, 3);
    section2(1, :) = [cal2(1), cal2(idx_split(1)), sum(prob2(1:idx_split(1)))];
    for k = 1:numel(idx_split)-1
        section2(k+1, :) = [cal2(idx_split(k)+1), cal2(idx_split(k+1)), ...
            sum(prob2(idx_split(k)+1:idx_split(k+1)))];
    end
    section2(end, :) = [cal2(idx_split(end)+1), cal2(end),...
        sum(prob2(idx_split(end)+1:end))];
else
    section2 = [cal2(1), cal2(end), sum(prob2)];
end
fprintf('95.4%% Probability:\n');
for k = 1:numel(idx_split)+1
    fprintf('%s (%.2f%%) %s \n', convertAge(floor(section2(k, 1))),...
        100*section2(k, 3), convertAge(ceil(section2(k, 2))));
end
