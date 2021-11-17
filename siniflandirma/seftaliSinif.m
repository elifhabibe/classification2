
%% SABÝTLER
clear; clc;
DizinResimler='03.Data';
SinifExcelDosya='SýnýfAdlarý2.xlsx';

EgitimSetiYuzde=80;
YSATekrarSayisi=5;%YSA modeli baþtan 5 kez tekrar edilecek demek. Nedeni ise YSA nýn baþlangýç þartlarýna baðlý olmasý.Yani sistem 5 kez deneyip en iyisini seçecek
 
%% DOSYALAR
%ResimDosyalar = getirResimDosyalar(DizinResimler);
[~,~,DataSinif] = xlsread(SinifExcelDosya); 
DataSinif=DataSinif(2:end,:);

[ResimSinif1, ResimSinif2, ResimSinif3] = sinifAta(DataSinif);
clear DataSinif ResimDosyalar SinifExcelDosya 

% kod içinde daha sonra kullanmayacaðým deðiþkenler, hafizade yer tutmasýn
% diye sildim zaten bu deðiþkeler ilerde kullanýlmayacak o yüzden sildim

%% EGITIM VE TEST DOSYALARIN ATANMASI
[EgitimResim,EgitimResimSinif,TestResim,TestResimSinif] = ysaGirisCikis(EgitimSetiYuzde,ResimSinif1,ResimSinif2,ResimSinif3); 
clear Resim* TrainingPercentage
%excel dosyasýndaki resim sýnýflarý bilgisi, ve eðitim yüzdesi parametrelerine göre resimleri, Eðitim ve Test olarak ikiye ayýrýyor.
%% YSA EÐÝTÝM
[EgitimData,EgitimSinif] = resim2data(DizinResimler,EgitimResim,EgitimResimSinif);
if ~isequal(size(EgitimData,2),size(EgitimSinif,2)), error('YSA EÐÝTÝM veri boyutlarý uyuþmuyor'); end
clear EgitimResim EgitimResimSinif EgitimSetiYuzde
%Herbir Eðitim Resim dosyasýnýn verisini ve sýnýf bilgisini ilgili parametreleri kullanarak okuyor
%resim2data fonk bu iþi yapýyor aþþaðýda tamýnladým

HataOpt=1e10;  % Baþlangýç için uydurma deðer 
for k=1:YSATekrarSayisi %Eðitim baþlýyor.Bu bloktaki iþlem 5 kez tekrarlanýyor demek ayný modeli, farklý "ilk aðýrlýk" deðerleri ile yeniden çalýþtýrmak
    net = patternnet(5); %3 satýrda YSA modeli tanýmlanýyor patternnet sýnýflandýrmaya özel að"patternnet(5)" demek, bir giriþ bir çýkýþ katmaný ve bir de gizli katmandan oluþuyor demek.gizli katmandaki möron sayýsý 5 demek Aðýrlýk deðerleri burda rastgele belirleniyor her seferinde.
    net.trainParam.showWindow=false;
    net = train(net,EgitimData,EgitimSinif);%oluþturulan net, eðitim data ve sýnýflarý ile eðitiliyor
    yy = net(EgitimData); %eðitilen net, eðitim datasý için bir sonuç veriyor
    % yy, herbir eðitim giriþi için 3 çýktý üretiyor. Bunlarýn herbiri 0-1 aralýðýnda.
    if ~isequal(unique(yy(:)),[0; 1])
        yy = sonuc2sinif(yy); %bu fonk aþþ tanýmladým. 100 001 gibi deðerler döndürüyor.Bunu yaparken de ilgili sýnýf için en yüksek "üyelik deðerini" kullanýyor.[0.6 0.1 0.3] için [1 0 0] birinci sýnýf
    end
    
    [~,HataYuzde] = TestDegerlendirme(EgitimSinif,yy);%Elimizde YSA çýktý deðerleri ve gerçek deðerler var."TestDegerlendirme" fonksiyonu da ikisini kýyaslayarak HataYuzde deðeri hesaplýyor Yani modelin, gerçeðe ne kadar uzak olduðunu buluyor
    %HataSTD=std(EgitimSinif(:)-yy(:));
    fprintf(1,'YSA Modelleme Tekrar: %02d  > Hata (%%): %.4f\n',k,HataYuzde);
    if HataYuzde<HataOpt %bu hata deðerini en küçük olacak þekilde kontrol ediyor ve eðer en küçükten küçük ise hem YSA modelini (netOpt olarak) hem de hatayý(HataOpt) saklýyor.en uygun modeli, en az hata ile buluyor.
        netOpt=net;
        HataOpt=HataYuzde;        
    end
    clear net yy HataSTD 
end %for un sonu 5 kere döncek
fprintf(1,'\n');
fprintf(1,'Min Hata(%%): %.4f\n\n',HataOpt);
yyEgitimOpt=netOpt(EgitimData);
yyEgitimOpt = sonuc2sinif(yyEgitimOpt);
if ~isequal(unique(yyEgitimOpt(:)),[0; 1])
    yyEgitimOpt = sonuc2sinif(yyEgitimOpt);
end
clear EgitimData EgitimSinif k HataOpt yyEgitimOpt YSATekrarSayisi k
%eðitilmiþ bir YSA modelimiz var.Bunu elimizdeki test verisi ile deneyebilriz.
%% YSA TEST
[TestData,TestSinif] = resim2data(DizinResimler,TestResim,TestResimSinif);% bu fonk ile  test resimlerine ait veri ve sýnýf bilgileri çaðrýlýyor
if ~isequal(size(TestData,2),size(TestSinif,2)), error('YSA TEST veri boyutlarý uyuþmuyor'); end
clear EgitimResim EgitimResimSinif
yyTahmin=netOpt(TestData);%test datasýný girdi olarak alarak test kümesindeki herbir resim için yyTahmin üretiyor.
if ~isequal(unique(yyTahmin(:)),[0; 1])
    yyTahmin = sonuc2sinif(yyTahmin);
end
[TestKiyasSonuc,HataYuzde] = TestDegerlendirme(TestSinif,yyTahmin);%elde edilen sonuçlar ile gerçek deðerler kýyaslanarak Hata miktarý elde ediliyor.burada elde edilen çýktý, herbir test resmine ait.

%% SONUÇLARIN YAZDIRILMASI
TestKiyasSonucStr=cell(size(TestKiyasSonuc));
for k=1:size(TestKiyasSonuc,1)
    for j=1:size(TestKiyasSonuc,2)-1
        if TestKiyasSonuc(k,j)==1% 1=k, 2=O, 3=B demek
            TestKiyasSonucStr{k,j}='Küçük';
        elseif TestKiyasSonuc(k,j)==2
            TestKiyasSonucStr{k,j}='Orta';
        elseif TestKiyasSonuc(k,j)==3
            TestKiyasSonucStr{k,j}='Büyük';
        else
            error('Sýnýf 1, 2, 3 den farklý!');
        end
    end
    TestKiyasSonucStr{k,3}=TestKiyasSonuc(k,3);
end
            
fprintf('\n');
fprintf(1,'YSA TEST SONUÇLARI:\n\n');
fprintf(1,'      Dosya    ORIGINAL SINIF    YSA SINIF   BAÞARI\n');
fprintf(1,'     -------  ---------------    ---------   -------\n');

for k=1:length(TestKiyasSonucStr)
    fprintf('%02d  %8s %10s %16s %8d\n', k, TestResim{k},TestKiyasSonucStr{k,1},TestKiyasSonucStr{k,2},TestKiyasSonucStr{k,3});
end

fprintf('\n');
fprintf('Tüm Test resimleri Baþarý Oraný(%%): %.2f [%d/%d] \n\n', 100.*(sum(TestKiyasSonuc(:,3))./numel(TestKiyasSonuc(:,3))), sum(TestKiyasSonuc(:,3)),numel(TestKiyasSonuc(:,3)));


function yySinif = sonuc2sinif(yy)
yySinif=zeros(size(yy));
for k=1:length(yy)
    [~,i]=max(yy(:,k));
    if isempty(i), error('MAX verisi belirlenemedi!'); end
    if ~isequal(size(i),[1 1]), error('MAX verisi birden fazla sýnýfa ait!'); end
    yySinif(i,k)=1;
    clear i
end
end

function [TestKiyasSonuc,HataYuzde] = TestDegerlendirme(TestSinif,yyTahmin)
if ~isequal(size(TestSinif),size(yyTahmin)), error('Test ve Tahmin veri boyutlarý ayný deðil!'); end

% TestKiyasSonuc=[Gerçek Test Sonuç];
TestKiyasSonuc=nan(length(TestSinif),3);
for k=1:length(TestSinif)
    p=find(TestSinif(:,k)==1);
    if isempty(p), error('TEST verisinde Sýnýf Belirtilmemiþ!'); end
    if ~isequal(size(p),[1 1]), error('TEST verisinde birden fazla sýnýf atanmýþ!'); end
    TestKiyasSonuc(k,1)=p;
    
    p=find(yyTahmin(:,k)==1);
    if isempty(p), error('TAHMÝN verisinde Sýnýf Belirtilmemiþ!'); end
    if ~isequal(size(p),[1 1]), error('TAHMÝN verisinde birden fazla sýnýf atanmýþ!'); end
    TestKiyasSonuc(k,2)=p;
    clear p
end
TestKiyasSonuc(:,3)=1; %varsayýlan olarak atanmasý
mask= TestKiyasSonuc(:,1) ~= TestKiyasSonuc(:,2);  
TestKiyasSonuc(mask,3)=0;
HataYuzde = 100 - (sum(TestKiyasSonuc(:,3))/numel(TestKiyasSonuc(:,3))).*100;
end

function [Data,Sinif] = resim2data(DizinResimler,Resim,ResimSinif)
Data1=imread(fullfile(DizinResimler,Resim{1}));
[N,M,T]=size(Data1); clear Data1;

Data=nan(N*M*T,length(Resim));

for k=1:length(Resim)
    A=imread(fullfile(DizinResimler,Resim{k}));
    if ~isequal(size(A),[N,M,T]), error('Resim boyutu standard deðil!'); end
    Data(:,k) = double(A(:));
    clear A
end
Sinif=ResimSinif';
end

function [EgitimResim,EgitimResimSinif,TestResim,TestResimSinif] = ysaGirisCikis(TrainingPercentage,ResimSinif1,ResimSinif2,ResimSinif3)

d=min([length(ResimSinif1) length(ResimSinif2) length(ResimSinif3)]);
nEgitim=floor(d*(TrainingPercentage/100));

% A.) EGITIM KUMESÝ
EgitimResim=[];
EgitimResimSinif=[];
% Sýnýf1 (K)  atamasý
EgitimResim=[EgitimResim; ResimSinif1(1:nEgitim)];
EgitimResimSinif=[EgitimResimSinif;repmat([ 1 0 0],nEgitim,1)];

% Sýnýf2 (O) atamasý
EgitimResim=[EgitimResim; ResimSinif2(1:nEgitim)];
EgitimResimSinif=[EgitimResimSinif;repmat([ 0 1 0],nEgitim,1)];

% Sýnýf3 (B) atamasý
EgitimResim=[EgitimResim; ResimSinif3(1:nEgitim)];
EgitimResimSinif=[EgitimResimSinif;repmat([ 0 0 1],nEgitim,1)];
if ~isequal(size(EgitimResim,1),size(EgitimResimSinif,1)), error('EðitimResim dosya ve sýnýf sayýsýnda sorun var'); end
if ~isequal(length(EgitimResim),length(unique(EgitimResim))), error('EðitimResim atamasýnda sorun var'); end

%-------------------------------------------------------------
% B.) TEST KUMESÝ
TestResim=[];
TestResimSinif=[];
% Sýnýf1 (K)  atamasý
TestFiles=ResimSinif1(nEgitim+1:end);
TestResim=[TestResim; TestFiles];
TestResimSinif=[TestResimSinif;repmat([ 1 0 0],length(TestFiles),1)];

% Sýnýf2 (O) atamasý
TestFiles=ResimSinif2(nEgitim+1:end);
TestResim=[TestResim; TestFiles];
TestResimSinif=[TestResimSinif;repmat([ 0 1 0],length(TestFiles),1)];

% Sýnýf3 (B) atamasý
TestFiles=ResimSinif3(nEgitim+1:end);
TestResim=[TestResim; TestFiles];
TestResimSinif=[TestResimSinif;repmat([ 0 0 1],length(TestFiles),1)];
if ~isequal(size(TestResim,1),size(TestResimSinif,1)), error('TestResim dosya ve sýnýf sayýsýnda sorun var'); end
if ~isequal(length(TestResim),length(unique(TestResim))), error('TestResim atamasýnda sorun var'); end
end

function [ResimSinif1, ResimSinif2, ResimSinif3] = sinifAta(DataSinif)
m1=0;
m2=0;
m3=0;
for k=1:length(DataSinif)
    if strcmpi(DataSinif{k,2},'K')
        m1=m1+1;
        ResimSinif1{m1,1}=DataSinif{k,1};
    elseif strcmpi(DataSinif{k,2},'O')
         m2=m2+1;
        ResimSinif2{m2,1}=DataSinif{k,1};
    elseif strcmpi(DataSinif{k,2},'B')
         m3=m3+1;
        ResimSinif3{m3,1}=DataSinif{k,1};
    else
        error('Sýnýf K, O, veya B deðerlerinden biri deðil!')
    end
end
end

function ResimDosyalar = getirResimDosyalar(DizinResimler)

% Training images
ResimDosyalar=dir(DizinResimler);
ResimDosyalar=ResimDosyalar(3:end,:);
ResimDosyalar=struct2cell(ResimDosyalar);
ResimDosyalar=ResimDosyalar(1,:);
ResimDosyalar=ResimDosyalar';
end