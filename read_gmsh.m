fid   = fopen('ex3.msh','r');%give filename
tline = fgetl(fid);%read 1st line

T.Edges      = [];
T.FBndyEdges = zeros(0,1);
T.EdgeEls    = zeros(0,2);

while ischar(tline)
    disp(tline)%display the line
    tline = fgetl(fid);%read another line
    % read points
    if strcmp(tline,'$Nodes')
        tline = fgetl(fid);
        N = str2num(tline);%number of nodes
        pt = zeros(N,2);
        for i=1:N        
            tline = fgetl(fid);
            data = str2num(tline);
            pt(i,1) = data(2);%x-coord of the node            
            pt(i,2) = data(3);%y-coord of the node
        end
    end
    
    % read elements
    if strcmp(tline,'$Elements')
        tline = fgetl(fid);
        N = str2num(tline);%number of lines and triangles
        iele=1;%number of triangles
        ibd_dirchlet = 1;%number of dirichlet edges
        ibd_nuemann  = 1;%number of nuemann edges
        for i=1:N
            tline = fgetl(fid);
            data = str2num(tline);
            if (data(2) == 1)%if it is edge 
                if (data(4)==1)% if it is Dirichlet
                    bd(ibd_dirchlet,:) = data([6,7]);%append to bd
                    ibd_dirchlet = ibd_dirchlet+1;
                elseif(data(4)==2) % if it is Nuemann
                    k = size(T.Edges,1)+1;
                    T.Edges(k,:) = data([6,7]);%append to T.Edges
                    l = size(T.FBndyEdges,1)+1;
                    T.FBndyEdges(l) = k;%append to T.FBndyEdges
                    T.EdgeEls(k,2)  = -l;%add info to T.EdgeEls at (k,2)
                end
            elseif data(2) == 2%if it is triangle
                ele(iele,:) = data([6,7,8]);%nodes of triangle
                iele = iele+1;
            end
        end
        break;
    end
end
fclose(fid);%close file

%%
T.Nodes = pt;
T.Elements = zeros(size(ele,1),3);
for i=1:size(ele,1)
    res = check_edges(T.Edges,ele(i,1),ele(i,2));%test the edge ele(i,1)---ele(i,2)
    T.Elements(i,1) = res(2);
    if(res(1)==0)%if it is a new edge
        T.Edges(res(2),:) = [ele(i,1),ele(i,2)];% add it to T.Edges
        T.EdgeEls(res(2),:) = [i,0];% i is the first triangle having it as an edge
    else% it is a known edge
        if(T.EdgeEls(abs(res(2)),2)==0)% if it is the initial value 0
            T.EdgeEls(abs(res(2)),2) = i;% i is the 2nd triangle having it as an edge
        else % if it is not 0, and it should be <0,
            T.EdgeEls(abs(res(2)),1) = i;% it is a nuemann edge, so i is the only one triangle having it as an edge
        end
    end
    
    res = check_edges(T.Edges,ele(i,2),ele(i,3));%test the edge ele(i,2)---ele(i,3)
    T.Elements(i,2) = res(2);
    if(res(1)==0)
        T.Edges(res(2),:) = [ele(i,2),ele(i,3)];
        T.EdgeEls(res(2),:) = [i,0];
    else
        if(T.EdgeEls(abs(res(2)),2)==0)
            T.EdgeEls(abs(res(2)),2) = i;
        else 
            T.EdgeEls(abs(res(2)),1) = i;
        end
    end
    
    res = check_edges(T.Edges,ele(i,3),ele(i,1));%test the edge ele(i,3)---ele(i,1)
    T.Elements(i,3) = res(2);
    if(res(1)==0)
        T.Edges(res(2),:) = [ele(i,3),ele(i,1)];
        T.EdgeEls(res(2),:) = [i,0];
    else
        if(T.EdgeEls(abs(res(2)),2)==0)
            T.EdgeEls(abs(res(2)),2) = i;
        else 
            T.EdgeEls(abs(res(2)),1) = i;
        end
    end
end
T.Elements = int32(T.Elements);
T.Edges = int32(T.Edges);
T.EdgeEls = int32(T.EdgeEls);
T.EdgeCFlags = int32(zeros(size(T.Edges,1),1));%straight edge(0) or curved(1)

T.NodePtrs = zeros(size(T.Nodes,1),1);
T.CNodePtrs= [];%constraint node
T.FNodePtrs= [];%free node
for i=1:size(bd,1)%check all nodes on dirichlet boundary
    for j=1:2
        if(T.NodePtrs(bd(i,j))==0)%if it is not considered
            if(check_node(T.Edges,T.EdgeEls,bd(i,j))==1)%it at least one edge if dirichlet
                k = length(T.CNodePtrs);
                T.CNodePtrs(k+1) = bd(i,j);% add it to T.CNodePtrs            
                T.NodePtrs(bd(i,j)) = -(k+1);%add info to T.NodePtrs
            end
        end
    end
end
for i=1:size(T.Nodes,1)
    if(T.NodePtrs(i)==0) % all rest nodes with 0 is free
        k = length(T.FNodePtrs);
        T.FNodePtrs(k+1) = i;% add it to T.FNodePtrs
        T.NodePtrs(i) = (k+1);%add info to T.NodePtrs 
    end
end
T.CNodePtrs=int32(T.CNodePtrs');
T.FNodePtrs=int32(T.FNodePtrs');
T.NodePtrs=int32(T.NodePtrs);
T.Degree =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = check_edges(edges,v1,v2)
% result=[a,b]
% a=0 means v1--v2 can be added into edges as a new edge. 
%       In this case, b=size(edges,1)+1>0
% otherwise, a=1:
%              if b>0, v1--v2 = edges[b,:]
%              if b<0, v2--v1 = edges[-b,:]


k = find(edges==v1);
Nedge = size(edges,1);

if(isempty(k))
    result = [0,Nedge+1];
    return;
end

for j=1:length(k)
    ie = k(j);
    if(ie<=Nedge && edges(ie,2)==v2)
        result = [1,ie];
        return;
    elseif(ie>Nedge && edges(ie-Nedge,1)==v2)
        result = [1,Nedge-ie];
        return;
    end    
end

result = [0,Nedge+1];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = check_node(edges,edge_ele,v1)

res = 0;
k = find(edges==v1);
Nedge = size(edges,1);

for j=1:length(k)%k should >= 2
    ie = k(j);
    if(ie<=Nedge)
        edge_index = ie;
    else
        edge_index = ie-Nedge;
    end
    if(edge_ele(edge_index,2)==0)
        res = 1;
        return
    end
end
end