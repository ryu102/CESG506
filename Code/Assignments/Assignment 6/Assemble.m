function [Psys, Fsys, Ksys] = Assemble(nodes, mesh, P, U,idxFree,idxFixed);

    %% define input format

    % %        nodeID  X       Y
    % nodes=[  1       -5.5    0.0;
    %          2        0.0    0.5;
    %          3        4.0    0.0];

    % %      Member       i       j        EA
    % mesh=[  1           1       2       2100;
    %         2           2       3       2100];

    % %          Node    X       Y
    % support=[  1       1       1;
    %            3       1       1];

    %% Do the assembly

    % initialize
    [nDOF, nNds] = size(nodes(:,2:end)');
    nDOF = 3;
    
    % initialize the internal force vector
    FsysC = mat2cell(zeros(nDOF*nNds,1),nDOF*ones(nNds,1));

    % create an empty stiffness matrix as a cell array of nDOF x nDOF matrices
    KsysC = mat2cell(zeros(nDOF*nNds),nDOF*ones(nNds,1),nDOF*ones(nNds,1));

    % list of initial nodal positions (coordinates)
    posList = nodes(:,2:end)';

    % element loop
    for e=1:length(mesh(:,1))
	i=mesh(e,2);
	j=mesh(e,3);
	EA=mesh(e,4);
    EI=mesh(e,5);

	% get element nodal force and nodal stiffness
	[fe,ke] = curvedbeam(posList(:,i), posList(:,j), U(:,i), U(:,j), EA, EI);

	% assemble element e to system
	FsysC{i} = FsysC{i} + fe(1:3);
	FsysC{j} = FsysC{j} + fe(4:6);

	KsysC{i,i} = KsysC{i,i} + ke(1:3,1:3);
	KsysC{i,j} = KsysC{i,j} + ke(1:3,4:6);
	KsysC{j,i} = KsysC{j,i} + ke(4:6,1:3);
	KsysC{j,j} = KsysC{j,j} + ke(4:6,4:6);
    end

    % flatten the stiffness matrix from cell array to matrix
    Psys = reshape(P,nDOF*nNds,1);
    Fsys=cell2mat(FsysC);
    Ksys=cell2mat(KsysC);

    % get index lists for free and supported DOFs
%     [idxFree, idxFixed] = GetIndices(nodes, support);

    % reduce the stiffness matrix and the redual to free DOFs only
    Psys(idxFixed)  =[];
    Fsys(idxFixed)  =[];
    Ksys(:,idxFixed)=[];
    Ksys(idxFixed,:)=[];

end
