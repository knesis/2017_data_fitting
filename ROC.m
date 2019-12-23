function [PPV, Az, PerfInd, ROC_threshold] = ROC(probdata0, probdata1, plotv)

if nargin<3
    plotv=true;
end

numN = length(probdata0);
numA = length(probdata1);

no_target = [zeros(1,numN)', probdata0];
target = [ones(1,numA)', probdata1];

alldata = [no_target; target];
alldata = sortrows(alldata, -2);

Nc = zeros(1,size(alldata,1))';
Nf = zeros(1,size(alldata,1))';

for i = 1:size(alldata,1)
    for j = 1:size(alldata,1)
        if alldata(j,2) > alldata(i,2)
            if alldata(j,1) == 0
               Nf(i) = Nf(i) + 1;
            else
               Nc(i) = Nc(i) + 1;
            end
        end
    end
end

Pd = Nc./numA;     % Prob of detection
Pf = Nf./numN;     % Prob of failure
Az = 0;
min_op = sqrt(max(Pf).^2 + max(Pd).^2);
OOP = [0,0];
Ptt = 0; Pth = 0;

for i = 2:size(Pf,1)
    Az = Az + Pd(i)*(Pf(i) - Pf(i-1));
end

for i = 2:size(Pf,1)
    for j = 2:size(Pd,1)
        dist = sqrt(Pf(i).^2 + (1-Pd(i)).^2);
        if dist < min_op
            min_op = dist;
            index = i;
            OOP = [Pf(i),Pd(i)];
        end
    end
end

ROC_threshold = alldata(index,2);

for i = 1:size(alldata,1)
    if (alldata(i,1) == 1) && (alldata(i,2) > ROC_threshold)
        Ptt = Ptt + 1;
    end
    if (alldata(i,2) > ROC_threshold)
        Pth = Pth + 1;
    end
end

PPV = (Ptt/Pth);

nA = numA:1:numN;
A1 = Az./(2-Az);
A2 = (2*Az.^2)./(1+Az);

avg = abs(mean(probdata0) - mean(probdata1));
std = sqrt(var(probdata0) + var(probdata1));
PerfInd = avg/std;

sigma = sqrt(((Az.*(1-Az)) + (nA-1).*(A1 - Az.^2) + ((1-nA)-1).*(A2 - Az.^2))./(nA.*(1-nA)));

% ROC Curve
if plotv==true
    figure;
    plot(Pf, Pd, 'r', 'LineWidth',2); hold on;
    plot(linspace(0,1,2),linspace(0,1,2),'--g');
    xlabel('1 - Specificity'); ylabel('Sensitivity'); title('ROC Curve');
    text(0.4,0.3,['Threshold = ',num2str(alldata(index,2))]);
    text(0.4,0.25,['PPV at Opt. Threshold is ',num2str(PPV)]);
    plot(OOP(1),OOP(2), 'b*', 'LineWidth', 2);
    legend(['ROC Curve,   ', 'A_z = ', num2str(Az)],'AUC = 0.5',...
        ['Optimal Operating Point = [',num2str(OOP(1)),', ',num2str(OOP(2)), ']'], 'Location', 'se');
else
   plot(Pf, Pd, 'LineWidth',2); hold on; 
end

end