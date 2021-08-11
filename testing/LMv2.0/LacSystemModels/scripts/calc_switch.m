function calc_fpt_ts(outputDirname, inputFilename, species, lowerbound, upperbound, i, j, k, l, m)

if nargin < 8 || nargin > 10
    error('Usage: calc_switch_pdf(outputDirname, inputFilename, species, lowerbound, upperbound, i, j, k, [l, m])');
end

d=nargin-5;
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

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loading %s into position %d,%d,%d,%d,%d',clock,inputFilename,i,j,k,l,m));

% Load the data.
ts=cast(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCountTimes',1)),'double');
data=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCounts',1)),[2,1]),'double');

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loaded data.',clock));

% Find all of the transitions.
sfpCrossings=[];
sfpCrossingIndex=[];
for x=[100:1:size(data,1)-1]
    if data(x,species) < lowerbound && data(x+1,species) >= lowerbound
        sfpCrossings(end+1) = 0;
        sfpCrossingIndex(end+1) = x+1;
    elseif data(x,species) >= lowerbound && data(x+1,species) < lowerbound
        sfpCrossings(end+1) = 0;
        sfpCrossingIndex(end+1) = x;
    elseif data(x,species) <= upperbound && data(x+1,species) > upperbound
        sfpCrossings(end+1) = 1;
        sfpCrossingIndex(end+1) = x;
    elseif data(x,species) > upperbound && data(x+1,species) <= upperbound
        sfpCrossings(end+1) = 1;
        sfpCrossingIndex(end+1) = x+1;
    end
end
lowToHigh=[];
highToLow=[];
for x=[1:size(sfpCrossings,2)-1]
    if sfpCrossings(x) == 0 && sfpCrossings(x+1) == 1
        lowToHigh(end+1,1:2)=[sfpCrossingIndex(x) sfpCrossingIndex(x+1)];
    elseif sfpCrossings(x) == 1 && sfpCrossings(x+1) == 0
        highToLow(end+1,1:2)=[sfpCrossingIndex(x) sfpCrossingIndex(x+1)];
    end
end

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Found %d transitions.',clock, length(lowToHigh)));

% Bin the data.
pdfs{i,j,k,l,m}=cell(4,2);
values=zeros(size(data,1),2);
valuesUsed=0;
for x=[1:size(lowToHigh,1)]
    for y=[lowToHigh(x,1):lowToHigh(x,2)]
        valuesUsed=valuesUsed+1;
        values(valuesUsed,1)=data(y,species);
        values(valuesUsed,2)=pdata(y,2);
    end
end
values=values(1:valuesUsed,:);
[N,X,Y]=hist2d(values(:,1),values(:,2),[0:max(values(:,1))],50);
pdfs{i,j,k,l,m}{1,1}=X;
pdfs{i,j,k,l,m}{2,1}=Y;
pdfs{i,j,k,l,m}{3,1}=N./(sum(sum(N))*(X(2)-X(1))*(Y(2)-Y(1)));
pdfs{i,j,k,l,m}{4,1}=lowToHigh;

% Bin the data.
values=zeros(size(data,1),2);
valuesUsed=0;
for x=[1:size(highToLow,1)]
    for y=[highToLow(x,1):highToLow(x,2)]
        valuesUsed=valuesUsed+1;
        values(valuesUsed,1)=data(y,species);
        values(valuesUsed,2)=pdata(y,2);
    end
end
values=values(1:valuesUsed,:);
[N,X,Y]=hist2d(values(:,1),values(:,2),[0:max(values(:,1))],50);
pdfs{i,j,k,l,m}{1,2}=X;
pdfs{i,j,k,l,m}{2,2}=Y;
pdfs{i,j,k,l,m}{3,2}=N./(sum(sum(N))*(X(2)-X(1))*(Y(2)-Y(1)));
pdfs{i,j,k,l,m}{4,2}=highToLow;

% Save the data.
cellsave(outputFilename, pdfs, d);
disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished.',clock));

