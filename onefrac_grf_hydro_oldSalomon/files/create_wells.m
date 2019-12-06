function [ node,elem,bdflag,fractures] = create_wells( node,elem,bdflag,fractures,positions,width)
%CREATE_WELLS Summary of this function goes here
%   Detailed explanation goes here

elem_t(:,1)=node(elem(:,1),1)+node(elem(:,2),1)+node(elem(:,3),1);
elem_t(:,2)=node(elem(:,1),2)+node(elem(:,2),2)+node(elem(:,3),2);
elem_t=elem_t/3;



for i=1:size(positions,1)
   center=positions(i,:);
   xtop=center(1)+width;
   xbot=center(1)-width;
   ytop=center(2)+width;
   ybot=center(2)-width;
   
   nodes_inin=node(:,1)>xbot&node(:,1)<xtop&node(:,2)>ybot&node(:,2)<ytop;
   triangles_in=elem_t(:,1)>xbot&elem_t(:,1)<xtop&elem_t(:,2)>ybot&elem_t(:,2)<ytop;
   elem(triangles_in,:)=[];
   elem_t(triangles_in,:)=[];
   bdflag(triangles_in,:)=[];
   
   old_idx=1:(length(node));
   node(nodes_inin,:)=[];
   new_idx=old_idx;
   new_idx(~nodes_inin)=1:(length(node));
   elem(:,1)=new_idx(elem(:,1));
   elem(:,2)=new_idx(elem(:,2));
   elem(:,3)=new_idx(elem(:,3));
   for j=1:length(fractures)
       fractures{j}=new_idx(fractures{j});
   end
   
   
   nodes_boundary=node(:,1)>=xbot&node(:,1)<=xtop&node(:,2)>=ybot&node(:,2)<=ytop;
   
   tmp=find(nodes_boundary);
   equality_elem=double(elem==tmp(1));
   for j=2:length(tmp)
       equality_elem=equality_elem+double(elem==tmp(j));
   end
   tmp2=sum(equality_elem,2);
   id_elem=tmp2==2;
   
   tmp=find(id_elem);
   
   for j=1:length(tmp)
       pattern=equality_elem(tmp(j),:);
       if pattern(1)==1&&pattern(2)==1
           bdflag(tmp(j),3)=i;
       end
       if pattern(1)==1&&pattern(3)==1
           bdflag(tmp(j),2)=i;
       end
       if pattern(3)==1&&pattern(2)==1
           bdflag(tmp(j),1)=i;
       end
       
       
   end
   
end

end

