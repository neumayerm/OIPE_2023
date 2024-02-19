function L = fem_make_reg_matrix_1st_order(inzid,index_in_matval,node)
% assembles the first order n x n regularization matrix L, where
% n is the number of finite elements. This type of regularization
% matrix incorporates some kind of smoothness assumption about the
% quantities to be reconstructed.
% inzid is the n x 3 connectivity matrix for triangular finite elements
% with linear shape functions.
%
% An entry into L is generated only if 2 vertices of triangles are identical.
%
% Bernhard Brandstätter (brand@ieee.org) 1998
% last modified: 2002

L=sparse(zeros(size(index_in_matval,1),size(index_in_matval,1)));

for i=1:size(index_in_matval,1)
    for j=i+1:size(index_in_matval,1)
      h=ismember(inzid(index_in_matval(i),:),inzid(index_in_matval(j),:));
      h1=sum(h);
      if h1==2
         L(i,j)=-1;
         L(j,i)=-1;
      end
    end
    
   L(i,i)=abs(sum(L(i,:)));
end

return

L = sparse(1, size(index_in_matval,1));
k = 1;
for i=1:size(index_in_matval,1)
  for j=(i+1):size(index_in_matval,1)
    h=ismember(inzid(index_in_matval(i),:),inzid(index_in_matval(j),:));
    h1=sum(h);
    if h1==2     
      p = node(inzid(index_in_matval(i),find(h==1)),:);
      a = norm(p(1,:)-p(2,:));    
      L(k,i)=a;
      L(k,j)=-a;
      k = k+1;
    end
  end
end
L = L./max(abs(L(:)));
return    
