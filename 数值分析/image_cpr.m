clear all;clc;
Image=imread('lena.jpg');
Image=double(Image);
Size = size(Image);
Image_cpr = uint8(zeros(Size));
L =512;
U1 = zeros(Size(1),L,3);
S1 = zeros(L,L,3);
V1 = zeros(Size(2),L,3);
for i = 1:3
    [U1(:,:,i),S1(:,:,i),V1(:,:,i)]=svds(Image(:,:,i),L);
    Image_cpr(:,:,i) = U1(:,:,i)*S1(:,:,i)*V1(:,:,i)';
end
imwrite(Image_cpr,'lena2.jpg');
imshow(Image_cpr);

