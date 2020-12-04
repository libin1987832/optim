matSize = [1000:floor(66146/5):76146];
fprintf(['Starting benchmarks with %d different matrix sizes ranging\n' ...
         'from %d-by-%d to %d-by-%d.\n'], ...
        length(matSize), matSize(1), matSize(1), matSize(end), ...
        matSize(end));
gflops = zeros(size(matSize));
for i = 1:length(matSize)
    gflops(i) = benchFcn(matSize(i));
    fprintf('Gigaflops: %f\n\n', gflops(i));
end
% results.matSize = matSize;
% results.gflops = gflops;