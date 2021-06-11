function calc_switch(outputFilename, inputFilename, species, lowerbound, upperbound, i, j, k, l, m)

if nargin < 8 || nargin > 10
    error('Usage: calc_switch(outputFilename, inputFilename, species, lowerbound, upperbound, i, j, k, [l, m])');
end

d=nargin-5;
if d == 3
    output=cell(i,j,k);
    l=1;
    m=1;
elseif d == 4
    output=cell(i,j,k,l);
    m=1;
elseif d == 5
    output=cell(i,j,k,l,m);
else
    error('Unsupported number of dimensions');
end

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loading %s into position %d,%d,%d,%d,%d',clock,inputFilename,i,j,k,l,m));


transitions=cell(100,1);
R=0;
try
    while 1
    
        % Try to load the next replicate.
        ts=cast(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCountTimes',R+1)),'double');
        data=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/SpeciesCounts',R+1)),[2,1]),'double');
        
        % If we loaded it, see if we need to extend the data set.
        R=R+1;
        disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loaded replicate %d',clock,R));
        if R > size(transitions,1)
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanding data size.',clock));
            newSize=size(transitions,1)*2;
            transitions(end+1:newSize)=cell(newSize-size(transitions,1),1);
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanded data size to %d.',clock,newSize));
        end
        
        % Find all of the transitions.
        crossingsFound=0;
        crossings=zeros(1000000,1);
        crossingIndex=zeros(1000000,1);
        for x=[1:size(data,1)-1]
            if data(x,species) < lowerbound && data(x+1,species) >= lowerbound
                crossingsFound=crossingsFound+1;
                crossings(crossingsFound) = 0;
                crossingIndex(crossingsFound) = x+1;
            elseif data(x,species) >= lowerbound && data(x+1,species) < lowerbound
                crossingsFound=crossingsFound+1;
                crossings(crossingsFound) = 0;
                crossingIndex(crossingsFound) = x;
            elseif data(x,species) <= upperbound && data(x+1,species) > upperbound
                crossingsFound=crossingsFound+1;
                crossings(crossingsFound) = 1;
                crossingIndex(crossingsFound) = x;
            elseif data(x,species) > upperbound && data(x+1,species) <= upperbound
                crossingsFound=crossingsFound+1;
                crossings(crossingsFound) = 1;
                crossingIndex(crossingsFound) = x+1;
            end
            if crossingsFound >= size(crossings,1)
                disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanding crossing data size.',clock));
                newSize=size(crossings,1)*2;
                crossings(end+1:newSize)=zeros(newSize-size(crossings,1),1);
                crossingIndex(end+1:newSize)=zeros(newSize-size(crossingIndex,1),1);
                disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanded crossing data size to %d.',clock,newSize));
            end
        end
        crossings=crossings(1:crossingsFound);
        crossingIndex=crossingIndex(1:crossingsFound);
        disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Found %d crossings.',clock, length(crossings)));
        
        transitions{R}=zeros(0,3);
        for x=[1:size(crossings,1)-1]
            if crossings(x) == 0 && crossings(x+1) == 1
                index=round((crossingIndex(x)+crossingIndex(x+1))/2);
                transitions{R}(end+1,:)=[1 index ts(index)];
            elseif crossings(x) == 1 && crossings(x+1) == 0
                index=round((crossingIndex(x)+crossingIndex(x+1))/2);
                transitions{R}(end+1,:)=[0 index ts(index)];
            end
        end
        disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Found %d transitions.',clock, length(transitions{R})));
        
    end
catch err
    if (strcmp(err.identifier,'MATLAB:hdf5readc:notAttributeOrDataset') == 0)
        rethrow(err);
    end
end
if R > 0
    % Remove any unused replicate columns.
    transitions=transitions(1:R);
    
    % Save the data.
    output{i,j,k,l,m}=transitions;
    cellsave(outputFilename, output, d);
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished %d replicates.',clock, R));
else
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): No replicates in file.',clock));
end

