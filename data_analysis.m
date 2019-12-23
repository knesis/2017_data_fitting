% Data Fitting and Analysis

% Assesses diagnostic significance of parametric fits on two-class data
% under various distributions, and estimates parameters.

% Computes statistical transformations of the estimated distributions, and
% evaluates system improvement.

%% Import data

file = 'probdata.xlsx'; % selects the data file

probdata0 = xlsread(file,'A:A'); % Data with target absent
probdata1 = xlsread(file,'B:B'); % Data with target present

%% ROC Curve

numN = length(probdata0); % Number absent
numA = length(probdata1); % Number present
[PPV, Az, PerfInd, ROC_threshold] = ROC(probdata0, probdata1);


%% Estimated Densities
   
% Determine regions of misclassification
Xi = linspace(0,13.5,200);
Xi_absent_shaded = Xi(Xi > ROC_threshold);
Xi_present_shaded = Xi(1:(200 - length(Xi_absent_shaded))+1);

% Fit distributions to classes
F_absent = ksdensity(probdata0, Xi);
F_present = ksdensity(probdata1, Xi);
F_absent_shaded = F_absent(length(Xi_present_shaded):200);
F_present_shaded = F_present(1:length(Xi_present_shaded));


%% ROC Sigma Comparison

nA = numA:1:numN;
A1 = Az./(2-Az); % False positive
A2 = (2*Az.^2)./(1+Az); % False negative

% Change in area under ROC as function of threshold
sigma = sqrt(((Az.*(1-Az)) + (nA-1).*(A1 - Az.^2) + ((1-nA)-1).*(A2 - Az.^2))./(nA.*(1-nA)));

%% Figure Display

% Estimated Densities
figure;
plot(Xi, F_absent, 'r-', Xi, F_present, 'b-'); hold on;
area(Xi_present_shaded, F_present_shaded, 'FaceAlpha', 0.4, 'FaceColor', 'blue');
area(Xi_absent_shaded, F_absent_shaded, 'FaceAlpha', 0.4, 'FaceColor', 'red');
xlabel('Input data'); ylabel('Estimated pdf');
title('Estimated Densities');
text(6,0.31, ['Threshold (optimal) = ', num2str(ROC_threshold)]);
text(6,0.285, ['PPV at Optimal Threshold = ', num2str(PPV)]);
text(6,0.25, 'BLUE - TYPE II error (False negative)');
text(6,0.225, 'RED - TYPE I error (False positive)');
legend('Normal (target absent, no disease etc.)', 'Abnormal (target present, disease present, etc.)');


% ROC Standard Deviation
figure;
plot(nA,sigma, 'r', 'LineWidth', 2);
xlabel('Number of samples with target');
ylabel('Standard Deviation of ROC area');
text(50,0.029, {['Number of normals: ', num2str(numN)];
    ['Number of abnormals: ', num2str(numA)];
    ['A_z = ', num2str(Az)];['\sigma(A_z) = ', num2str(sigma(1))]});
title({'Performance Increase, Decline in \sigma(A_z)';
    'Hanley and McNeill, Radiology, Vol. 143, No. 1, pp.29-36, Apr 1982'});

%% PDF Fits and Chi-Squared Tests

% Density fits - No Target
ray0 = dist_test(probdata0, 'Rayleigh');
nak0 = dist_test(probdata0, 'Nakagami');
gam0 = dist_test(probdata0, 'Gamma');
ric0 = dist_test(probdata0, 'Rician');

% Density fits - Target Present
ray1 = dist_test(probdata1, 'Rayleigh');
nak1 = dist_test(probdata1, 'Nakagami');
gam1 = dist_test(probdata1, 'Gamma');
ric1 = dist_test(probdata1, 'Rician');


%% Best Fit Density & Summary of All Chi-Squared Tests

% No Target
pdf0 = [ray0.pdf, nak0.pdf, gam0.pdf, ric0.pdf];
chi0 = [ray0.st, nak0.st, gam0.st, ric0.st];

% Retrieve index of best (minimum) chi-squared fit
ind0 = find([chi0.chi2stat] == min([chi0.chi2stat]));

figure;
set(gca,'Visible','off')
text(0.3, 1, "Hypothesis testing");
text(-0.05,0.75, {'        Target Absent (70)       ';
    '';
    ['Best fit density: ',num2str(pdf0(ind0).DistributionName)];
    ['Chi-squared value: ', num2str(chi0(ind0).chi2stat)]
    ['Degrees of Freedom: ', num2str(chi0(ind0).df)];
    'Gamma Parameters: ';
    ['A: ', num2str(gam0.pdf.a)];
    ['B: ', num2str(gam0.pdf.b)];
    });

% Target Present
pdf1 = [ray1.pdf, nak1.pdf, gam1.pdf, ric1.pdf];
chi1 = [ray1.st, nak1.st, gam1.st, ric1.st];

% Find index of smallest chi-squared value
ind1 = find([chi1.chi2stat] == min([chi1.chi2stat]));

text(0.55,0.75, {'        Target Present (30)       ';
    '';
    ['Best fit density: ',num2str(pdf1(ind1).DistributionName)];
    ['Chi-squared value: ', num2str(chi1(ind1).chi2stat)]
    ['Degrees of Freedom: ', num2str(chi1(ind1).df)];
    'Rician Parameters: ';
    ['S: ', num2str(ric1.pdf.s)];
    ['Sigma: ', num2str(ric1.pdf.sigma)];
    });

names = ["Rayleigh", "Nakagami", "Gamma", "Rician"]';
h0 = [ray0.h, nak0.h, gam0.h, ric0.h]';
df0 = [chi0.df]';
st0 = [chi0.chi2stat]';

text(0.25, 0.45, "Summary of all Chi-Squared Tests");
text(-0.05, 0.25, {"         Target Absent         ";("Density        h  DoF   Chi^2 ");""});
text(-0.05, 0.1, names);
text(0.125, 0.1, num2str(h0));
text(0.175, 0.1, num2str(df0));
text(0.225, 0.1, num2str(st0));

h1 = [ray1.h, nak1.h, gam1.h, ric1.h]';
df1 = [chi1.df]';
st1 = [chi1.chi2stat]';

text(0.55, 0.25, {"         Target Present         ";("Density        h  DoF   Chi^2 ");""});
text(0.55, 0.1, names);
text(0.725, 0.1, num2str(h1));
text(0.775, 0.1, num2str(df1));
text(0.825, 0.1, num2str(st1));

%% Random Data

newprob0 = random(pdf0(ind0).DistributionName, gam0.pdf.a, gam0.pdf.b, 70, 2);
newprob1 = random(pdf1(ind1).DistributionName, ric1.pdf.s, ric1.pdf.sigma, 30, 2);

% Transformations - New Data (Absent)
AM_prob0 = 0.5*(sum(newprob0,2));
MAX_prob0 = max(newprob0,[],2);
GM_prob0 = sqrt(prod(newprob0,2));

% Transformations - New Data (Present)
AM_prob1 = 0.5*(sum(newprob1,2));
MAX_prob1 = max(newprob1,[],2);
GM_prob1 = sqrt(prod(newprob1,2));

trans_prob = {newprob0(:,1), AM_prob0, GM_prob0, MAX_prob0;
    newprob1(:,1), AM_prob1, GM_prob1, MAX_prob1};

titles = ["Data Supplied", "Arithmetic Mean", "Maximum", "Geometric Mean"];

figure;
for i = 1:size(trans_prob,2)
    subplot(2,2,i);
    for j = 1:size(trans_prob,1)
        plot(Xi,ksdensity(trans_prob{j,i},Xi), 'LineWidth', 2);
        hold on;
    end
    title(titles(i));
    xlabel('Data'); ylabel('Estimated PDF');
    L = legend('Absent', 'Present', 'Location', 'ne');
    L.set('FontSize', 7);
    hold off;
end

%% ROC Curves (Transformations)

figure;
[PPV_MAX, Az_MAX, PerfInd_MAX] = ROC(MAX_prob0,MAX_prob1,false);
[PPV_GM, Az_GM, PerfInd_GM] = ROC(GM_prob0,GM_prob1,false); 
[PPV_AM, Az_AM, PerfInd_AM] = ROC(AM_prob0,AM_prob1,false);

plot(linspace(0,1,2),linspace(0,1,2),'--g');
legend(['ROC Curve MAX, A_z = ', num2str(Az_MAX)], ...
    ['ROC Curve GM, A_z = ', num2str(Az_GM)], ...
    ['ROC Curve AM, A_z = ', num2str(Az_AM)], 'Location', 'se');

%% Summary

figure;
set(gca,'Visible','off')
text(0.4, 1.05, "Summary");
text(-0.1,0.8, {'Hypothesis testing- Target Absent (70)';
    ['Best fit density: ',num2str(pdf0(ind0).DistributionName)];
    ['Chi-squared value: ', num2str(chi0(ind0).chi2stat)]
    ['Degrees of Freedom: ', num2str(chi0(ind0).df)];
    'Gamma Parameters: ';
    ['A: ', num2str(gam0.pdf.a)];
    ['B: ', num2str(gam0.pdf.b)];
    });

text(0.55,0.8, {'Hypothesis testing- Target Present (30)';
    ['Best fit density: ',num2str(pdf1(ind1).DistributionName)];
    ['Chi-squared value: ', num2str(chi1(ind1).chi2stat)]
    ['Degrees of Freedom: ', num2str(chi1(ind1).df)];
    'Rician Parameters: ';
    ['S: ', num2str(ric1.pdf.s)];
    ['Sigma: ', num2str(ric1.pdf.sigma)];
    });

text(-0.1, 0.5, {'Performance Measures';
    ['PPV at optimal threshold = ', num2str(PPV), ',  Performance Index = ', num2str(PerfInd), ',  Az = ', num2str(Az)]});

text(-0.1, 0.2, {'Post-Processing Results';
    'PPV at optimal threshold';
    [num2str(PPV_AM),' (Arithmetic Mean), ',num2str(PPV_GM), ' (Geometric Mean), ',num2str(PPV_MAX) ' (Maximum)'];
    'Performance Index';
    [num2str(PerfInd_AM),' (Arithmetic Mean), ',num2str(PerfInd_GM), ' (Geometric Mean), ',num2str(PerfInd_MAX) ' (Maximum)'];
    'Area under the ROC Curve';
    [num2str(Az_AM),' (Arithmetic Mean), ',num2str(Az_GM), ' (Geometric Mean), ',num2str(Az_MAX) ' (Maximum)']});

