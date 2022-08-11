% Example 4: Comparisons with one PTV and one OAR

clear all; close all; clc;

% Add data and functions to path
currentFolder = pwd;
cd ..
addpath(genpath(pwd));
cd(currentFolder);

% PTV - prostate
prostate.name = 'PTV_68';
prostate.terms = {struct('type','unif','dose',81,'weight',1),...
    struct('type','ldvc','dose',81,'percent',5,'weight',1)};

% OAR - rectum
rectum.name = 'Rectum';
rectum.terms = {struct('type','udvc','dose',30,'percent',30,'weight',1)};

% Create problem instance
fprintf('\nExample 4\n');
structs = {prostate,rectum};
prob = FluenceMapOpt(structs,'tol',1e-2);

% Calculate approximate doses
fprintf('\nCalculating approximate doses\n');
labels = ['a' 'b' 'c'];
for ii = 1:length(labels)
    fprintf('\nExample 4%s\n\n',labels(ii));
    if ii == 1
        prob.calcBeams();
    elseif ii == 2
        prob.calcBeamsIter();
    else
        prob.calcBeamsSlack();
    end
    fprintf('\nIterations: %d, Time: %.2f\n',prob.nIter,prob.time);
    prob.saveResults(['ex4' labels(ii) 'Approx.mat']);
end

% Calculate polished doses
fprintf('\nCalculating polished doses\n');
for ii = 1:length(labels)
    fprintf('\nExample 4%s\n',labels(ii));
    if ii == 1
        load('ex4aApprox.mat');
        x = results.x;
        prob.calcBeamsPolish(x);
        t1 = prob.time;
    elseif ii == 2
        x0 = prob.x0;
        prob.calcBeamsPolish(x0);
        t1 = prob.time;
    else
        prob.calcBeamsConvex();
        t1 = prob.time;
        prob.calcBeamsPolish(prob.x);
        t1 = t1 + prob.time;
    end
    fprintf('\nObjective: %.4e, Time: %.2f\n',prob.getObj('unif'),t1);
    prob.saveResults(['ex4' labels(ii) 'Polish.mat']);
end

% Calculating approximate dose with continuation
fprintf('\nCalculating approximate dose with continuation\n')
prob = calcBeamsContinue(prob,structs,0.99,0.01,100,true,0,false);
fprintf('Iterations: %d, Time: %.2f\n',prob.nIter,prob.time);
prob.saveResults('ex4Continue.mat');