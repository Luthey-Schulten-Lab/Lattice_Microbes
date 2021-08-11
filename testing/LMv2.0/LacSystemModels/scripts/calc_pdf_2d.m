function calc_pdf(outputFilename, inputFilename, species, i, j, k, l, m)

if nargin < 6 || nargin > 8
    error('Usage: calc_pdf(outputDirname, inputFilename, species, i, j, k, [l, m])');
end

d=nargin-3;
if d == 3
    pdfs=cell(i,j,k);
    l=1;
    m=1;
elseif d == 4
    pdfs=cell(i,j,k,l);
    m=1;
elseif d == 5
    pdfs=cell(i,j,k,l,m);
else
    error('Unsupported number of dimensions');
end

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loading %s into position %d,%d,%d',clock,inputFilename,i,j,k));

counts=[];
R=0;
try
    while 1
    
        % Try to load the next replicate.
        data=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCounts',R+1)),[2,1]),'double');
    
        % If we loaded a data set, process it.
        R=R+1;
        if mod(R,1000) == 0, disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Processing replicate %d',clock,R));, end        
        first=6;
        last=length(data)-2;
        counts(end+1:end+1+last-first,1)=data(first:last,species);
    end
catch
end
if R > 0
    [C,X]=hist(counts,[0:max(counts)]);
    pdfs{i,j,k,l,m}=zeros(2,length(X));
    pdfs{i,j,k,l,m}(1,:)=X;
    pdfs{i,j,k,l,m}(2,:)=[C./(sum(C)*(X(2)-X(1)))];
    cellsave(outputFilename, pdfs, d);
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished %d replicates.',clock, R));
else
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): No replicates in file.',clock));
end

