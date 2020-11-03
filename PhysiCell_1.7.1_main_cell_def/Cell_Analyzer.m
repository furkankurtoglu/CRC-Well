cd output

s=what;
MatFiles = s.mat;
OutMatFiles = MatFiles(contains(MatFiles,'physicell'));
OutMatFiles(1) = [];
OutMatFiles(1) = [];
% for i = 1:length(OutMatFiles)
%     OutMatFiles{i}=OutMatFiles{i}(1:14);
% end
cell_per_organoid=zeros(length(OutMatFiles),1);

for i=1:length(OutMatFiles)
    load(OutMatFiles{i})
    total_cell_count = size(cells,2);
    cell_per_organoid(i) = total_cell_count/250;
end
plot(cell_per_organoid)
