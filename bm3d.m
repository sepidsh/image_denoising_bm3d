clc;clear;
img = imread('lena.jpg');
img = imresize(img,[218,218]);
img = padarray(img,[19 19],'symmetric','both');
z = imnoise(img,'gaussian', 0, 0.016);
figure(2);
imshow(z(19:end-19,19:end-19,:));
d=psnr(uint8(z),img)
title('noisy image');
%%
%%%%%%%%%%%%%%%%%%%%%%%% BM3d %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Hard theresholding %%%%%%%%%%%%%%%%%%%%%
sigma=sqrt(0.016); 
Nht1= 8;  %parameters intialization 
N_ht_2= 16;
N_ht_step= 4;
Nhts= 19;
N_ht_es= 1;
n_ht_pr= 0;
betha_ht= 2;
lambda_2d= 0;
lambda_3d= 2.7;
th_ht_math= 2500;
image_size = size(z);
ntw = 8;
numerator_h = zeros(image_size(1),image_size(2));
denomerator_h = zeros(image_size(1),image_size(2));
for ch= 1:3
for i = Nhts+1:image_size(1)-Nhts
    i
    for j = Nhts+1:image_size(2)-Nhts
        window = z(i-Nhts:i+Nhts, j-Nhts:j+Nhts,ch);
        blocks = im2col(window(:,:), [Nht1 Nht1], 'sliding');
        d = zeros(size(blocks,2),1);
        d = d + 50000;
        for k =1:6:size(blocks,2)
            tmp =double(blocks(:,529))-double(blocks(:,k));
            tmp = reshape(tmp, [Nht1 Nht1]);
            tmp = norm(tmp,2);
            tmp = tmp.*tmp;
            d(k) = tmp/(Nht1*Nht1);
        end
        [sorted_d, I] = sort(d);
        I = I(1:ntw);
        blocks = blocks(:, I);
        P = zeros([Nht1 Nht1 ntw]);
        for k = 1:ntw
            P(:,:,k) = reshape(blocks(:,k), [Nht1 Nht1]);
        end
        P = dctn(P);
        P = wthresh(P,'h',sigma*lambda_3d);
        wP = 1/sum(P(:)>0);
        fP = idctn(P);
        for k = 1:ntw
			q = fP(:,:,k);
            tmpx = max(1,i-Nhts) + floor((529-1)/32);
            tmpy = max(1,j-Nhts) + (mod(529-1,32));
            numerator_h(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) = numerator_h(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) + (wP * q);
            denomerator_h(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) = wP + denomerator_h(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1);
        end
    end
end 
u_basic(:,:,ch) = numerator_h./denomerator_h;

end

imwrite(uint8(u_basic(Nhts+1:end-Nhts, Nhts+1:end-Nhts,:)),'res_phase1.jpg');
save('tmp');
load('tmp');
u_basic=u_basic(Nhts+1:end-Nhts, Nhts+1:end-Nhts,:);

a=psnr(uint8(u_basic),img(Nhts+1:end-Nhts, Nhts+1:end-Nhts,:))
u_basic = padarray(u_basic,[19 19],'symmetric','both');
%%

numerator_wi = zeros(image_size(1),image_size(2));
denomerator_wi = zeros(image_size(1),image_size(2));
for ch =1:3
for i = Nhts+1:image_size(1)-Nhts
    step=i
    for j = Nhts+1:image_size(2)-Nhts
        window = z(i-Nhts:i+Nhts, j-Nhts:j+Nhts,ch);
        window2 = u_basic(i-Nhts:i+Nhts, j-Nhts:j+Nhts,ch);
        blocks = im2col(window(:,:), [Nht1 Nht1], 'sliding');
        blocks2 = im2col(window2(:,:), [Nht1 Nht1], 'sliding');
        d = zeros(size(blocks,2),1);
        d=d+50000;
        for k  =1:3:size(blocks,2)
            tmp = wthresh(double(blocks(:,529)),'h',sigma*lambda_2d)-wthresh(double(blocks(:,k)),'h',sigma*lambda_2d);
            tmp = reshape(tmp, [Nht1 Nht1]);
            tmp = norm(tmp,2);
            tmp = tmp.*tmp;
            d(k) = tmp/(Nht1*Nht1);
        end
        [sorted_d, I] = sort(d);
        I = I(1:ntw);
        blocks = blocks(:, I);
        blocks2 = blocks2(:, I);
        P = zeros([Nht1 Nht1 ntw]);
        P_basic = zeros([Nht1 Nht1 ntw]);
        for k = 1:ntw
            P(:,:,k) = reshape(blocks(:,k), [Nht1 Nht1]);
            P_basic(:,:,k) = reshape(blocks2(:,k), [Nht1 Nht1]);
        end
        
        P_basic = dctn(P_basic);
        for k = 1:ntw
            tmp = P_basic(:,:,k);
            tmp = norm(tmp);
            tmp = tmp.*tmp;
            wP(k) = tmp/(tmp+(sigma*sigma));
        end
        P = dctn(P);
        for k = 1:ntw
            P(:,:,k) = P(:,:,k)*wP(k);
        end
        P = idctn(P);
        wP = wP.*wP;
        wP = sum(wP(:));
        wP = 1/wP;
        for k = 1:ntw
            q = P(:,:,k);
                tmpx = max(1,i-Nhts) + floor((529-1)/32);
                tmpy = max(1,j-Nhts) + (mod(529-1,32));
                numerator_wi(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) = numerator_wi(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) + (wP * q);
                denomerator_wi(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1) = wP + denomerator_wi(tmpx:tmpx+Nht1-1 , tmpy:tmpy+Nht1-1);
        end
    end
end 
result(:,:,ch) = numerator_wi./denomerator_wi;
end
result = result(Nhts+1:end-Nhts, Nhts+1:end-Nhts,:);
imwrite(uint8(result),'res_phase2.jpg');
b=psnr(uint8(result),img(Nhts+1:end-Nhts, Nhts+1:end-Nhts));
save('tmp2');