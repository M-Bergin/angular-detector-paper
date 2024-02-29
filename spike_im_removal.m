function [filtered_image,L]=spike_im_removal(im_in,n_IQR)

%Dirty function to perform a basic image filter and remove spikes
%MB 25/9/20

%Size of the image
N=size(im_in);


%Find the spikes

r = iqr(im_in(:));      %IQR of the data

m = median(im_in(:),"omitnan");  %Median of the data 

%Original Code line (MB) - above updated for Matlab2022 based on help
%recommendation (AF - 25-07-23)
%m=nanmedian(im_in(:));     %Median of the data 

%Test histogram
% figure;histogram(im_in(:))
% hold on
% xline(m,'LineWidth',1)
% xline(m+n_IQR*r,'LineWidth',1)

%n_IQR=6;
%Find points much higher than the rest of the data
inds=find(im_in>m+n_IQR*r);

%Find corresponding locations in the matrix
[row, col]=ind2sub(N,inds);

filtered_image=im_in;


L=length(inds);

%Original code from MB - version below the commented section replaces
%'nanmean' with the recommended 'mean' function with options to omit NaN
%values activated. AF, 25-07-23

% for l=1:L
%     %Take average over adjacent points if box is still in the image
%     if row(l)==1
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l):row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);
%     elseif row(l)==N(1)
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l),col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);       
%     elseif col(l)==1
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l):col(l)+1)']);
%     elseif col(l)==N(2)
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l))']);
%     else
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);
%     end
%     
% end
% 
% %Apply again to ensure they are all removed correctly
% 
% for l=1:L
%     %Take average over adjacent points if box is still in the image
%     if row(l)==1
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l):row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);
%     elseif row(l)==N(1)
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l),col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);       
%     elseif col(l)==1
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l):col(l)+1)']);
%     elseif col(l)==N(2)
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l))']);
%     else
%         filtered_image(row(l),col(l))=NaN;
%         filtered_image(row(l),col(l))=nanmean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)']);
%     end
%     
% end


for l=1:L
    %Take average over adjacent points if box is still in the image
    if row(l)==1
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l):row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");
    elseif row(l)==N(1)
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l),col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");       
    elseif col(l)==1
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l):col(l)+1)'],"omitnan");
    elseif col(l)==N(2)
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l))'],"omitnan");
    else
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");
    end
    
end

%Apply again to ensure they are all removed correctly

for l=1:L
    %Take average over adjacent points if box is still in the image
    if row(l)==1
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l):row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");
    elseif row(l)==N(1)
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l),col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");       
    elseif col(l)==1
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l):col(l)+1)'],"omitnan");
    elseif col(l)==N(2)
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l))'],"omitnan");
    else
        filtered_image(row(l),col(l))=NaN;
        filtered_image(row(l),col(l))=mean([filtered_image(row(l)-1:row(l)+1,col(l));filtered_image(row(l),col(l)-1:col(l)+1)'],"omitnan");
    end
    
end




