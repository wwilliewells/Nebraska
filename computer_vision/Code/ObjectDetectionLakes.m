function row = ObjectDetectionLakes(image,divisor)
close all
clc
% 6900 - 6905,6907 - 10, 6944 - 53, 6957 - 58, 6965 - 89, 6991, 6993,6996,
% 7010 - 7017
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC6900','JPEG'));6900
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC6958','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC6981','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC6987','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC6993','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC7024','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera1/MFDC7083','JPEG'));7102
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC6973','JPEG'));6973
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC6991','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7003','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7021','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7037','JPEG'));7037
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7043','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7053','JPEG'));
%colorImage=im2single(imread('/ImageAnalysis/cameraImages/Camera2/MFDC7099','JPEG'));7173
limg=imread(sprintf('/ImageAnalysis/cameraImages/Camera1/MFDC%s',image),'JPEG');%7083
timg=im2single(limg);%03,21,37
% figure,imshow(limg)

% grayscale
img=rgb2gray(timg);

[m,n]=size(img);
img=img(1:m-100,:);
[m,n]=size(img);
figure, imshow(img)

% g = fspecial('log', 11, 0.5);
% 
% f = imfilter(img,g,'conv'); 
% figure, imshow(f)


simg = img(m/divisor -m/(2*divisor):m/divisor+m/(2*divisor),:);
[ms,ns]= size(simg);
cimg = edge(simg,'canny',0.99,1.0);
figure, imshow(cimg)
mi = 0;
for i=1:1:ms
    if sum(cimg(i,:)) > mi 
        mi = sum(cimg(i,:));
        row = i;
    end
end
row+m/divisor -m/(2*divisor)
simg = img(row+m/divisor -m/(2*divisor):m,:);
figure, imshow(simg)

% %% mean shift
% % W = window % sum_{x\in W} (xH(x)) % H(x) = histogram % N(x) = H(x) /
% % sum(H(x))
% h=imhist(f+img);
% hs=sum(h);
% % normalized histogram
% Nx=h/hs;
% 
% segment=zeros(m,n);
% wn=1;
% wm=1;
% ws=3;
% for i=1:1:m
%     wn=1;
%     for j=1:1:n
%         for k=wm:1:wm+ws
%             for p=wn:1:wn+ws
%                 segment(i,j)=segment(i,j)+img(k,p)*Nx(floor(img(k,p)*255)+1);
%             end
%         end
%         if j > ws
%             wn=wn+1;
%         end
%     end
%     if i > ws
%         wm=wm+1;
%     end
% end
% 
% color=0;
% flag=1;
% num=-1;
% pcolor=0;
% xcount=0;
% 
% while flag==1
%     if xcount == n*m
%         break;
%     end
%     flag2=0;
%     for i=1:1:m
%         for j=1:1:n
%             if segment(i,j)< 1 && flag2 == 0
%                 num=segment(i,j);
%                 color=color+1;
%                 flag2=1;
%             end
%             if segment(i,j) <= num + 0.0004 && segment(i,j) >= num - 0.0004
%                 segment(i,j)=round(color*3.2278);
%             end
%         end
%     end
%     if pcolor==color
%         flag=0;
%     end
%     pcolor=color;
%     xcount=xcount+1;
% end
% segment=segment/255;
% figure, imshow(segment)
% 
% %%
% % Majority vote row as horizon
% row = m;
% row2 = 1;
% wmax = 0;
% wmin = inf;
% pcount = m;
% irow=zeros(3,1);
% 
% mcount = 0;
% hcount=0;
% for i = floor(3*m/8):1:ceil(m-3*m/8)
%     ucount = sum(segment(i,:) - segment(i-1,:)); 
%     dcount=sum(segment(i,:)-segment(i+1,:));
%     if  abs(ucount + dcount) <  wmin && abs(ucount) > wmax 
%         
%             wmin = abs(ucount + dcount);
%             wmax = abs(ucount);
%  
%             irow(3) = irow(2);
%             irow(2) = irow(1);
%             irow(1) = i
%             [ucount, dcount, wmin, wmax]
%     end
% 
%         %         for j = 1:1:n
% %             if j == 1
% %                 hcount=0;
% %                 ccount = 1;
% %             elseif segment(i,j) - segment(i,j-1) == 0
% %                 ccount = ccount + 1;
% %             else
% %                 if ccount > hcount
% %                     hcount = ccount;
% %                 end
% %                 ccount = 0;
% %             end
% %         end
% %             row = i;
% %             wmin = abs(ucount);
% %         if hcount > mcount
% %             mcount = hcount;
% %             row = i;
% %             wmin = abs(ucount);
% %         end
%         
%         %row2 = i; %[hcount ccount i]
% %     elseif ucount > wmax && i < m - m/4 && i > m/4 %&& abs(ucount + dcount) < 0.5 && abs(pcount) < 1
% %         for j = 1:1:n
% %             if j == 1
% %                 hcount=0;
% %                 ccount = 1;
% %             elseif segment(i,j) - segment(i,j-1) == 0
% %                 ccount = ccount + 1;
% %             else
% %                 if ccount > hcount
% %                     hcount = ccount;
% %                 end
% %                 ccount = 0;
% %             end
% %         end
% %         wmax = ucount;
% %         
% %         row2 = i; [2 hcount ccount i]
%     %end %
%     %pcount = ucount;
%     
% end
% 
% mrow = abs(irow - m/divisor);
% if (mrow(3) < mrow(2) && mrow(3) < mrow(1)) 
%     row = irow(3);
%     pc =3;
% elseif mrow(2) < mrow(3) && mrow(2) < mrow(1) && irow(2) < m/divisor
%     row = irow(2);pc = 2;
% else
%     row = irow(1);pc =1;
% end
% 
% if irow(3) < m/divisor && pc ~= 3
%     if pc == 1 && irow(1) > m/divisor
%         row = irow(3);
%     elseif pc == 2 && irow(2) > m/divisor
%         row = irow(3);
%     end
% elseif irow(2) < m/divisor && pc ~= 2
%     if pc == 1 && irow(1) > m/divisor
%         row = irow(2);
%     elseif pc == 3 && irow(3) > m/divisor
%         row = irow(2);
%     end
% elseif irow(1) < m/2 && pc ~= 1
%     if pc == 3 && irow(3) > m/2
%         row = irow(1);
%     elseif pc == 2 && irow(2) > m/2
%         row = irow(1);
%     end
% end
% figure, imshow(img(row:m,:))
% %row =(row+row2)/2
% 
% %% Smooth uneven horizon
% ssi = 0;
% ssimg = img(row-30:row+30,:); % img(row-30:row+30,:);
% [ssm,ssn] = size(ssimg);
% for j = 1:1:ssn
%     maxi = 0;
%     if j == 1
%         for k = 1:1:ssm - 1
%             if ssimg(k,j) - ssimg(k+1,j) < maxi
%                 maxi= ssimg(k,j) - ssimg(k+1,j);
%                 ssi = k + 1;% [maxi,ssi]
%             end
%         end
%         ssip = ssi;
%         ssi = 1;
%     else
%         for k = ssip - 3 :1: ssip +3
%             if k > 0 && k < ssm - 1
%                 if ssimg(k,j) - ssimg(k+1,j) < maxi
%                     maxi = ssimg(k,j) - ssimg(k+1,j);
%                     ssi = k + 1; % [k,j,maxi,ssi]
%                 end
%             end
%         end
%         ssip = ssi;
%     end
%     if ssip > 1
%         img(row-30:row-30+ssip-1,j) = img(row-30+ssip,j);
%     end
% end   
% figure, imshow(ssimg)
% 
% simg = img(row:m,:);
% figure, imshow(simg)
%% edge detection
cimg = edge(simg,'canny',0.57,0.6);
figure, imshow(cimg)

%% search for majority edges
swn=50;
vswn=swn;

[cm,cn]=size(cimg);
bu=cm;
br=cn;
box=[1,1];
wmax=-Inf;
for i=1:1:cm-swn
    if mod(i,10) == 0
        vswn=vswn+6;
    end
    for j=1:1:cn-swn
        if i+vswn-1 < cm-swn && j+vswn-1 < cn-swn
            bimg=cimg(i:i+vswn-1,j:j+vswn-1);
            bvu=vswn-1;bvr=vswn-1;
        elseif i+vswn-1 < cm - swn && j+vswn - 1 >=cn-swn
            bimg=cimg(i:i+vswn-1,j:cn);
            bvu=vswn-1;bvr=cn-j;
        elseif j+vswn-1 < cn - swn && i+vswn - 1 >=cm-swn
            bimg=cimg(i:cm,j:j+vswn-1);
            bvu=cm-i;bvr=vswn-1;
        else
            bimg=cimg(i:cm,j:cn);
            bvu=cm-i;bvr=cn-j;
        end
        ucount=log(sum(sum(bimg == 1))/(vswn^2));
        if ucount > wmax
            box=[i,j];
            bu=bvu;
            br=bvr; % [count box;wmax bu br];
            wmax=ucount;
        end
    end
end
bimg=cimg(box(1):box(1)+bu,box(2):box(2)+br);
%figure, imshow(bimg)
[bu,br]

%% search for connected edges
truncate=[1,swn];
flag=false;
ucount=0;
pcount=0;
[bm,bn]=size(bimg);
for bi=bm:-1:1
    anticount=sum(bimg(bi,:) == 1);
    if anticount > 0
        pcount = pcount + 1;
    else
        pcount = 0;
    end
    if pcount > 5 
        flag=true;% object identified
        truncate(1)=1;
        truncate(2)=bi+pcount;
    end
    ucount=sum(bimg(bi,:) == 0);
    if ucount == br+1 && flag == true
        truncate(1)=bi;
        break;
    end
end
[box(1) + m - cm,br]

%% Results
bd = box(1) + m - cm + truncate(1);
bu = box(1) + m - cm + truncate(2);
bl = box(2);%+n-cn
br = box(2) + br;
%[bu-bd br-bl truncate(1) truncate(2);bu bd br bl];

figure, imshow(timg);
hold on

rectangle('Position',[bl bd br-bl bu-bd], 'EdgeColor','c');
hold off
