% get input A as ref 
%B = [1,2,3,4];
A = [3,2,1,4,3,2,5];
%A = [1,2,4,3,1,4,4,2,1];
% get input B as contig
%A = [1,2,2,3,4];
B = [3,1,4,4,1,2,1];
%B = [1,2,3,2,1,4,2,1]
% sequence converted in numbers A->1, C -> 2, G -> 3, T -> 4, U -> 5
nA = length(A);
nB = length(B);
% first row and first colmn
NWMatrix = zeros(nB+1,nA+1);

NWMatrix(:,1) = 0:-1:-nB;
NWMatrix(1,:) = 0:-1:-nA;

% Creating the Needleman Wunsh Algo matrix 
for iB=2:nB+1
    for iA=2:nA+1
        % Indel
        MIM = [NWMatrix(iB-1,iA)-1,NWMatrix(iB,iA-1)-1]; 
        
        if A(iA-1)-B(iB-1) == 0
           % Match
            MIM(3) = NWMatrix(iB-1,iA-1)+1;
        else
            % Missmatch
            MIM(3) = NWMatrix(iB-1,iA-1)-1;
        end
        [a,b] = max(MIM);
        NWMatrix(iB,iA) = a;
        NWindx(iB-1,iA-1) = b;
        clear a b 
    end
    MIM = 0; 
end

% Trace back to find the best match 
iB = nB;
iA = nA; 
alnA = [];
alnB = [];
while iB>0 && iA>0
     if NWindx(iB,iA)==1
         alnA(length(alnA)+1) = 0;
         alnB(length(alnB)+1) = B(iB);
         iB = iB-1;
     elseif NWindx(iB,iA)==2
         alnA(length(alnA)+1) = A(iA);
         alnB(length(alnB)+1) = 0;
         iA = iA-1;
     elseif NWindx(iB,iA)
         alnA(length(alnA)+1) = A(iA);
         alnB(length(alnB)+1) = B(iB);
         iA = iA-1; 
         iB = iB-1;
     end
end

alnA = alnA(end:-1:1);
alnB = alnB(end:-1:1);


%
Start = round(10000*rand(10,2));
End = Start+ round(1000*rand(10,2));
Stop=0;
iPath = 1;
PossiblePath = num2cell(1:length(Start));
PathExtension = ones(1,length(Start));
while Stop == 0
    ExtensionPathIndx = find(PathExtension==1);
    nPossiblePath = length(PossiblePath);
    if ~isempty(ExtensionPathIndx)
        PresentPath = PossiblePath{ExtensionPathIndx(1)};
        PresentDiagEnd = End(PresentPath(end),:);
        PotentialDiag = Start(:,2)>PresentDiagEnd(1,2) & Start(:,1) > PresentDiagEnd(1,1);
        lenPotential = sum(PotentialDiag);
        if lenPotential ==0 
            PathExtension(ExtensionPathIndx(1)) = 0;
        else   # [ [1], [2], [3] ]
            PossiblePath(ExtensionPathIndx(1)+lenPotential:nPossiblePath+lenPotential-1) = PossiblePath(ExtensionPathIndx(1)+1:nPossiblePath) ; # EXTEND! # [ [1], [2], [2], [3]]
            trash = repmat(PossiblePath{ExtensionPathIndx(1)},lenPotential,1);  # repmat() --> [[1], [1]]
            trash = [trash find(PotentialDiag==1)];  # find [[1 3], [1 4]]
            [row,col] = size(trash);  # replaces here, outputs [[1 3], [1 4], [2], [3]]
            PossiblePath(ExtensionPathIndx(1):ExtensionPathIndx(1)+lenPotential-1) = mat2cell(trash,repmat(1,1,lenPotential),col);
        end
    else
     Stop = 1;
    end
    clear trash 
end

nPossiblePath = length(PossiblePath); 
StartPositionContig = 1;
EndPositionContig = 10000;
for i=1:nPossiblePath
    Path = PossiblePath{i};
    PathStarts = Start(Path,:);
    PathEnds = End(Path,:);
    Score(i) = PathStarts(1,1)-StartPositionContig+sum(PathStarts(2:end,1)-PathEnds(1:end-1,1))+sum(PathStarts(2:end,2)-PathEnds(1:end-1,2))+ ...
        EndPositionContig-PathEnds(end,1);
end
