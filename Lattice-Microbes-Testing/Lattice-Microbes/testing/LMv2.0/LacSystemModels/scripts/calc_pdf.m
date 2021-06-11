function calc_pdf(outputFilename, inputFilename, species, numberReplicates, numberTimeSteps, i, j, k, l, m)

if nargin < 8 || nargin > 10
    error('Usage: calc_pdf(outputDirname, inputFilename, species, numberReplicates, numberTimeSteps, i, j, k, [l, m])');
end

d=nargin-5;
if d == 3
    pdfs=cell(i,j,k);
elseif d == 4
    pdfs=cell(i,j,k,l);
elseif d == 5
    pdfs=cell(i,j,k,l,m);
else
    error('Unsupported number of dimensions');
end

Rs=[1:numberReplicates];

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loading %s into position %d,%d,%d',clock,inputFilename,i,j,k));

counts=zeros(length(Rs)*numberTimeSteps,1);
countsUsed=0;
shortRuns=0;
for ri=[1:length(Rs)]
    R=Rs(ri);
    if mod(R,1000) == 0, disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Processing replicate %d',clock,R));, end
    data=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCounts',R)),[2,1]),'double');
    if length(data) < numberTimeSteps, shortRuns=shortRuns+1;, end
    first=6;
    last=length(data)-2;
    counts((countsUsed+1):(countsUsed+1+last-first),1)=data(first:last,species);
    countsUsed=countsUsed+last-first+1;
end
disp(sprintf('Short runs: %d',shortRuns));
[C,X]=hist(counts(1:countsUsed,1),[0:max(counts(1:countsUsed,1))]);

if d == 3
    pdfs{i,j,k}=zeros(2,length(X));
    pdfs{i,j,k}(1,:)=X;
    pdfs{i,j,k}(2,:)=[C./(sum(C)*(X(2)-X(1)))];
elseif d == 4
    pdfs{i,j,k,l}=zeros(2,length(X));
    pdfs{i,j,k,l}(1,:)=X;
    pdfs{i,j,k,l}(2,:)=[C./(sum(C)*(X(2)-X(1)))];
elseif d == 5
    pdfs{i,j,k,l,m}=zeros(2,length(X));
    pdfs{i,j,k,l,m}(1,:)=X;
    pdfs{i,j,k,l,m}(2,:)=[C./(sum(C)*(X(2)-X(1)))];
end
cellsave(outputFilename, pdfs, d);
disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished.',clock));

