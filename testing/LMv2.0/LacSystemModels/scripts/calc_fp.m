function calc_fp(outputDirname, pdfDirname, i, j, k, l, m)

if nargin < 5 || nargin > 7
    error('Usage: calc_fp(outputDirname, pdfDirname, i, j, k, [l, m])');
end

d=nargin-2;
if d == 3
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finding fixed points in pdf %d at position %d,%d,%d',clock,pdfDirname,i,j,k));
    fps=cell(i,j,k);
    pdf=cellload_single(pdfDirname,i,j,k);
elseif d == 4
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finding fixed points in pdf %d at position %d,%d,%d',clock,pdfDirname,i,j,k,l));
    fps=cell(i,j,k,l);
    pdf=cellload_single(pdfDirname,i,j,k,l);
elseif d == 5
    disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finding fixed points in pdf %d at position %d,%d,%d',clock,pdfDirname,i,j,k,l,m));
    fps=cell(i,j,k,l,m);
    pdf=cellload_single(pdfDirname,i,j,k,l,m);
else
    error('Unsupported number of dimensions');
end

showPlots=1;

fp1=-1;
regionEnd=200;
if regionEnd > size(pdf,2), regionEnd=size(pdf,2);, end
if size(find(pdf(2,1:regionEnd)>0,1,'first'),2) > 0
    regionStart=find(pdf(2,1:regionEnd)>0,1,'first');
    if regionEnd-regionStart>=40 && size(find(pdf(2,regionStart:regionEnd)>0),2) >= 0.8*(regionEnd-regionStart)
        fpi=findlocalmax(pdf(1,:),pdf(2,:),regionStart,regionEnd);
        if fpi > 0
            fp1=pdf(1,fpi);
        end
    end
end

fp3=-1;
regionStart=1000;
if size(pdf,2) > regionStart
    regionEnd=size(pdf,2);
    if regionEnd-regionStart>=40
        fpi=findlocalmax(pdf(1,:),pdf(2,:),regionStart,regionEnd);
        if fpi > 0
            fp3=pdf(1,fpi);
        end
    end
end

fp2=-1;
if fp1>0 && fp3>0
    fpi=findlocalmin(pdf(1,:),pdf(2,:),fp1+1,fp3+1, 20, 1);
    if fpi > 0
        fp2=pdf(1,fpi);
    end
end

if d == 3
    fps{i,j,k}=[fp1 fp2 fp3];
elseif d == 4
    fps{i,j,k,l}=[fp1 fp2 fp3];
elseif d == 5
    fps{i,j,k,l,m}=[fp1 fp2 fp3];
end

cellsave(outputDirname, fps, d);
disp(sprintf('(%04d-%02d-%02d %02d:%02d:%05.02f): Finished.',clock));


if showPlots == 1
    figure
    loglog(pdf(1,:),pdf(2,:), '-');
    if fp1 ~= -1
        vline(fp1, 'g');
    end
    if fp2 ~= -1
        vline(fp2, 'r');
    end
    if fp3 ~= -1
        vline(fp3, 'g');
    end
end


