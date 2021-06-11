function calc_pdf_2d(outputFilename, inputFilename, species, parameter, i, j, k, l, m)

if nargin < 7 || nargin > 9
    error('Usage: calc_pdf(outputDirname, inputFilename, species, i, j, k, [l, m])');
end

d=nargin-4;
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

valuesUsed=0;
values=zeros(1000000,2);
R=0;
try
    while valuesUsed < 50000000
    
        % Try to load the next replicate.
        data=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCounts',R+1)),[2,1]),'double');
        data2=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/ParameterValues/%07d',R+1,parameter)),[2,1]),'double');
        if size(data,1) ~= size(data2,1)
            error('Species and parameter values sizes do not match.');
        end
        
        % If we loaded a data set, process it.
        R=R+1;
        if mod(R,1000) == 0, disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Processing replicate %d',clock,R));, end        
        first=6;
        last=length(data)-2;
        nextValuesUsed=valuesUsed+last-first+1;
        while size(values,1) < nextValuesUsed
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanding data size.',clock));
            newSize=size(values,1)*2;
            values(end+1:newSize,:)=zeros(newSize-size(values,1),2);
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanded data size to %d.',clock,size(values,1)));
            drawnow('update');
        end
        values(valuesUsed+1:nextValuesUsed,1)=data(first:last,species);
        values(valuesUsed+1:nextValuesUsed,2)=data2(first:last,2);
        valuesUsed=nextValuesUsed;
    end
catch err
    if (strcmp(err.identifier,'MATLAB:hdf5readc:notAttributeOrDataset') == 0)
        rethrow(err);
    end
end
if R > 0
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Binning data.',clock));
    [N,X,Y]=hist2d(values(1:valuesUsed,1),values(1:valuesUsed,2),[0:max(values(:,1))],50);
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Done binning data.',clock));
    pdfs{i,j,k,l,m}=cell(3,1);
    pdfs{i,j,k,l,m}{1}=X;
    pdfs{i,j,k,l,m}{2}=Y;
    pdfs{i,j,k,l,m}{3}=N./(sum(sum(N))*(X(2)-X(1))*(Y(2)-Y(1)));
    cellsave(outputFilename, pdfs, d);
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished %d replicates.',clock, R));
else
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): No replicates in file.',clock));
end

