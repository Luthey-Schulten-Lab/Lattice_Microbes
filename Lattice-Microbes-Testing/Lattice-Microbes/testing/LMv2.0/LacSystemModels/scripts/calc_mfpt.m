function calc_mfpt(outputFilename, inputFilename, species, tf, p1, p2, p3, i, j, k, l, m)

if nargin < 10 || nargin > 12
    error('Usage: (outputDirname, inputFilename, species, tf, p1, p2, p3, i, j, k, l, m)');
end

d=nargin-7;
if d == 3
    mfpt=cell(i,j,k);
    l=1;
    m=1;
    onOffFlag=k;
elseif d == 4
    mfpt=cell(i,j,k,l);
    m=1;
    onOffFlag=l;
elseif d == 5
    mfpt=cell(i,j,k,l,m);
    onOffFlag=m;
else
    error('Unsupported number of dimensions');
end

p1i=p1+1; p2i=p2+1; p3i=p3+1; pMax=2*p3; pMaxi=pMax+1;

disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Loading %s into position %d,%d,%d,%d,%d',clock,inputFilename,i,j,k,l,m));
fpts=zeros(pMaxi,1000);
truncated=ones(pMaxi,1000);
R=0;
try
    while 1
    
        % Try to load the next replicate.
        counts=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/FirstPassageTimes/%02d/Counts',R+1,species-1)),[2,1]),'double');
        fpt=cast(permute(hdf5read(inputFilename,sprintf('/Simulations/%07d/FirstPassageTimes/%02d/Times',R+1,species-1)),[2,1]),'double');
        
        % If we loaded it, see if we need to extend the data set.
        R=R+1;
        if mod(R,1000) == 0, disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Processing replicate %d',clock,R));, end
        if R > size(fpts,2)
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanding data size.',clock));
            newSize=size(fpts,2)*2;
            fpts(:,end+1:newSize)=zeros(pMaxi,newSize-size(fpts,2));
            truncated(:,end+1:newSize)=ones(pMaxi,newSize-size(truncated,2));
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Expanded data size to %d.',clock,newSize));
        end
        
        % Process the replicate.
        minCount=counts(1); maxCount=counts(end);
        if maxCount > pMax, fpt=fpt(1:end-(maxCount-pMax));, maxCount=pMax;, end
        if onOffFlag == 1 && maxCount < p3
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Run %d was truncated: max value %d < %d',clock,R,maxCount,p3));
            fpts(:,R)=tf;
        elseif onOffFlag == 2 && minCount > p1
            disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Run %d was truncated: min value %d > %d',clock,R,minCount,p1));
            fpts(:,R)=tf;
        else
            fpts(:,R)=max(fpt);
        end
        minIndex=minCount+1; maxIndex=maxCount+1;
        fpts(minIndex:maxIndex,R)=fpt;
        truncated(minIndex:maxIndex,R)=0;
    end
catch
end
if R > 0

    % Remove any unused replicate columns.
    fpts=fpts(:,1:R);
    truncated=truncated(:,1:R);
    
    mfpt{i,j,k,l,m}=zeros(7,pMaxi);
    mfpt{i,j,k,l,m}(1,:)=[0:pMax];
    mfpt{i,j,k,l,m}(2,:)=mean(fpts,2);
    mfpt{i,j,k,l,m}(3,:)=max(fpts,[],2);
    mfpt{i,j,k,l,m}(4,:)=min(fpts,[],2);
    for r=[1:size(fpts,1)]
        [tmpm,tmpc]=expfit(fpts(r,:),0.01,truncated(r,:));
        mfpt{i,j,k,l,m}(5,r)=tmpm;
        mfpt{i,j,k,l,m}(6,r)=tmpc(1);
        mfpt{i,j,k,l,m}(7,r)=tmpc(2);
    end
    cellsave(outputFilename, mfpt, d);
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished %d replicates.',clock, R));
else
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): No replicates in file.',clock));
end

