
% -A script for merging multiple edi files and deleting the noisy outlies. 
% -This function is written initially for the outputs from SigMT.
% -INPUT: files content the frequencies,impedances and tippers and their
% corresponding errors.
% -Dependence:calc_MT  calc_Z 
% Copyright 2023 Feng Jiang, South China Sea Institute of Oceanology, CAS.
close all; clear all; clc;
[file,filepath,indx] = uigetfile('*.txt','select multi-freq-bands files:','MultiSelect','on');
data = [];
if ~iscell(file)
    nfile = 1;
    data = load([filepath,file]);
elseif iscell(file)
    nfile = length(file);
    for ifile = 1:nfile
        idata = load([filepath,file{ifile}]);
        data = [data;idata];
    end
else
    error('You didnot choice any file!');
end

flist = data(:,1);
Z = nan(2,2,length(flist));
Zvar = nan(2,2,length(flist));
T = nan(2,1,length(flist));
Tvar = nan(2,1,length(flist));
Z(1,1,:) = data(:,2) + 1i*data(:,3); %Zxx
Z(1,2,:) = data(:,4) + 1i*data(:,5); %Zxy
Z(2,1,:) = data(:,6) + 1i*data(:,7); %Zyx
Z(2,2,:) = data(:,8) + 1i*data(:,9); %Zyy
Zvar(1,1,:) = data(:,10); %zxxvar
Zvar(1,2,:) = data(:,11); %zxyvar
Zvar(2,1,:) = data(:,12); %zyxvar
Zvar(2,2,:) = data(:,13); %zyyvar
T(1,1,:) = data(:,14) + 1i*data(:,15); %Tx
T(2,1,:) = data(:,16) + 1i*data(:,17); %Ty
Tvar(1,1,:) = data(:,18); %txvar
Tvar(2,1,:) = data(:,19); %tyvar

mu = 4*pi*10^-7;
Z_SI = Z*(mu*1000);
Zvar_SI = Zvar*(mu*1000);
f = flist;
[rhoa,rhoaerr,phs,phserr] = calc_MT(Z_SI,Zvar_SI,f);

rhoa_xy = squeeze(rhoa(1,2,:));
rhoa_yx = squeeze(rhoa(2,1,:));
rhoa_xyerr = squeeze(rhoaerr(1,2,:));
rhoa_yxerr = squeeze(rhoaerr(2,1,:));
phs_xy = squeeze(phs(1,2,:));
phs_yx = squeeze(phs(2,1,:));
phs_xyerr = squeeze(phserr(1,2,:));
phs_yxerr = squeeze(phserr(2,1,:));
tx = squeeze(abs(T(1,1,:)));
ty = squeeze(abs(T(2,1,:)));
txvar = squeeze(Tvar(1,1,:));
tyvar = squeeze(Tvar(2,1,:));

ax = figure(1);
ax(1) = subplot(12,1,1:3);
ax_rhoxy = errorbar(flist,rhoa_xy,rhoa_xyerr,'ro'); hold on
ax_rhoyx = errorbar(flist,rhoa_yx,rhoa_yxerr,'bo');
set(gca,'xdir','reverse','xscale','log','yscale','log');
ylabel('App.Res.')
grid on;
ax(2) = subplot(12,1,4:9);
ax_phsxy = errorbar(flist,phs_xy,phs_xyerr,'ro'); hold on
ax_phsyx = errorbar(flist,phs_yx,phs_yxerr,'bo');
set(gca,'xdir','reverse','xscale','log');
ylim([-180,180]);
grid on;
ylabel('Phase');
ax(3) = subplot(12,1,10:12);
ax_tx = errorbar(flist,tx,txvar,'ro'); hold on
ax_ty = errorbar(flist,ty,tyvar,'bo');
set(gca,'xdir','reverse','xscale','log','yscale','log');
xlabel('Frequency(Hz)');
ylabel('Tipper');
grid on;

%rho
ax_rhoxy.XDataSource = 'flist';
ax_rhoxy.YDataSource = 'rhoa_xy';
ax_rhoxy.YNegativeDeltaSource = 'rhoa_xyerr';
ax_rhoyx.XDataSource = 'flist';
ax_rhoyx.YDataSource = 'rhoa_yx';
ax_rhoyx.YNegativeDeltaSource = 'rhoa_yxerr';
%phs
ax_phsxy.XDataSource = 'flist';
ax_phsxy.YDataSource = 'phs_xy';
ax_phsxy.YNegativeDeltaSource = 'phs_xyerr';
ax_phsyx.XDataSource = 'flist';
ax_phsyx.YDataSource = 'phs_yx';
ax_phsyx.YNegativeDeltaSource = 'phs_yxerr';
%tips
ax_tx.XDataSource = 'flist';
ax_tx.YDataSource = 'tx';
ax_tx.YNegativeDeltaSource = 'txvar';
ax_ty.XDataSource = 'flist';
ax_ty.YDataSource = 'ty';
ax_ty.YNegativeDeltaSource = 'tyvar';

%backup for re-do
bkflist = flist;
bkrho_xy = rhoa_xy;
bkrho_yx = rhoa_yx;
bkrho_xyerr = rhoa_xyerr;
bkrho_yxerr = rhoa_yxerr;
bkphs_xy = phs_xy;
bkphs_yx = phs_yx;
bkphs_xyerr = phs_xyerr;
bkphs_yxerr = phs_yxerr;
bktx = tx;
bkty = ty;
bktxvar = txvar;
bktyvar = tyvar;

pts = [];
comp = [];
count = 1;
click = 1;
while click~=0
    [x,y,button] = ginput(1);
    if button == 1
        disrhoxy = nan;
        disrhoyx = nan;
        disphsxy = nan;
        disphsyx = nan;
        distx = nan;
        disty = nan;
        [~,axnum] = ismember(gca,ax);
        if axnum == 1
            disrhoxy = sqrt(abs(log10(flist)-log10(x)).^2 + abs(log10(rhoa_xy)-log10(y)).^2);
            disrhoyx = sqrt(abs(log10(flist)-log10(x)).^2 + abs(log10(rhoa_yx)-log10(y)).^2);
        elseif axnum == 2
            disphsxy = sqrt(abs(log10(flist)-log10(x)).^2 + abs(deg2rad(phs_xy-y)).^2);
            disphsyx = sqrt(abs(log10(flist)-log10(x)).^2 + abs(deg2rad(phs_yx-y)).^2);
        elseif axnum == 3
            distx = sqrt(abs(log10(flist)-log10(x)).^2 + abs(tx-y).^2);
            disty = sqrt(abs(log10(flist)-log10(x)).^2 + abs(ty-y).^2);
        else
            warning('!!!Wrong Click!!!')
        end
        
        mindis = nan(6,2);
        [mindis(1,1),mindis(1,2)]= min(disrhoxy); %[minrhoxy,indrhoxy]
        [mindis(2,1),mindis(2,2)]= min(disrhoyx); %[minrhoyx,indrhoyx]
        [mindis(3,1),mindis(3,2)]= min(disphsxy);%[minphsxy,indphsxy]
        [mindis(4,1),mindis(4,2)]= min(disphsyx); %[minphsyx,indphsyx]
        [mindis(5,1),mindis(5,2)]= min(distx); %[mintx,indtx]
        [mindis(6,1),mindis(6,2)]= min(disty); %[minty,indty]
        [mintarget,mincomp] = min(mindis,[],1);
        ind = mindis(mincomp(1),2);
        pts(count) = ind;
        comp(count) = mincomp(1);
        disp(['ind = ',num2str(ind),'; comp = ',num2str(mincomp(1))]);
        count = count + 1;
        %flist(pts) = nan;
        for i = 1:length(pts)
            if comp(i) == 1||comp(i) == 3
                rhoa_xy(pts(i)) = nan;
                rhoa_xyerr(pts(i)) = nan;
                phs_xy(pts(i)) = nan;
                phs_xyerr(pts(i)) = nan;
            elseif comp(i) == 2||comp(i) == 4
                rhoa_yx(pts(i)) = nan;
                rhoa_yxerr(pts(i)) = nan;
                phs_yx(pts(i)) = nan;
                phs_yxerr(pts(i)) = nan;
            elseif comp(i) ==5||comp(i)==6
                tx(pts(i)) = nan;
                ty(pts(i)) = nan;
                txvar(pts(i)) = nan;
                tyvar(pts(i)) = nan;
            end
        end
        refreshdata;
        drawnow;
    elseif button == 3 % undo the delete clicks.
        if ~isempty(pts)
            %flist(pts(end)) = bkflist(pts(end));
            if comp(end) == 1||comp(end) == 3
                rhoa_xy(pts(end)) = bkrho_xy(pts(end));
                rhoa_xyerr(pts(end)) = bkrho_xyerr(pts(end));
                phs_xy(pts(end)) = bkphs_xy(pts(end));
                phs_xyerr(pts(end)) = bkphs_xyerr(pts(end));
            elseif comp(end) ==2||comp(end) == 4
                rhoa_yx(pts(end)) = bkrho_yx(pts(end));
                rhoa_yxerr(pts(end)) = bkrho_yxerr(pts(end));
                phs_yx(pts(end)) = bkphs_yx(pts(end));
                phs_yxerr(pts(end)) = bkphs_yxerr(pts(end));
            elseif comp(end) == 5 || comp(end) == 6
                tx(pts(end)) = bktx(pts(end));
                ty(pts(end)) = bkty(pts(end));
                txvar(pts(end)) = bktxvar(pts(end));
                tyvar(pts(end)) = bktyvar(pts(end));
            end
            refreshdata; 
            drawnow;
            pts = pts(1:end-1);
            comp = comp(1:end-1);
            count = count - 1;
        else
            pts = [];
            count = 1;
        end
    else
        click = 0; % exit by click the middle key on the mouse
    end
end


%----------------------------------------------------
rhoa(1,1,isnan(rhoa_xy)) = nan;
rhoaerr(1,1,isnan(rhoa_xy)) = nan;
phs(1,1,isnan(phs_xy)) = nan;
phserr(1,1,isnan(phs_xy)) = nan;
rhoa(2,2,isnan(rhoa_yx)) = nan;
rhoaerr(2,2,isnan(rhoa_yx)) = nan;
phs(2,2,isnan(phs_yx)) = nan;
phserr(2,2,isnan(phs_yx)) = nan;

rhoa(1,2,:) = rhoa_xy;
rhoaerr(1,2,:) = rhoa_xyerr;
rhoa(2,1,:) = rhoa_yx;
rhoaerr(2,1,:) = rhoa_yxerr;

T_new = T; Tvar_new = Tvar;
T_new(1,1,isnan(tx)) = nan + 1i*nan;
T_new(2,1,isnan(ty)) = nan + 1i*nan;
Tvar_new(1,1,isnan(txvar)) = nan;
Tvar_new(2,1,isnan(tyvar)) = nan;

[Z_new,Zerr_new] = calc_Z(rhoa,rhoaerr,phs,phserr,flist);
% mean
f = unique(flist);
for i = 1:length(f)
    Z_final(1,1,i) = mean(squeeze(Z_new(1,1,flist == f(i))),'omitnan');
    Z_final(1,2,i) = mean(squeeze(Z_new(1,2,flist == f(i))),'omitnan');
    Z_final(2,1,i) = mean(squeeze(Z_new(2,1,flist == f(i))),'omitnan');
    Z_final(2,2,i) = mean(squeeze(Z_new(1,1,flist == f(i))),'omitnan');
    Zerr_final(1,1,i) = mean(squeeze(Zerr_new(1,1,flist == f(i))),'omitnan');
    Zerr_final(1,2,i) = mean(squeeze(Zerr_new(1,2,flist == f(i))),'omitnan');
    Zerr_final(2,1,i) = mean(squeeze(Zerr_new(2,1,flist == f(i))),'omitnan');
    Zerr_final(2,2,i) = mean(squeeze(Zerr_new(1,1,flist == f(i))),'omitnan');
    
    T_final(1,1,i) = mean(squeeze(T_new(1,1,flist == f(i))),'omitnan');
    T_final(2,1,i) = mean(squeeze(T_new(2,1,flist == f(i))),'omitnan');
    Tvar_final(1,1,i) = mean(squeeze(Tvar_new(1,1,flist == f(i))),'omitnan');
    Tvar_final(2,1,i) = mean(squeeze(Tvar_new(2,1,flist == f(i))),'omitnan');
end


figure(2) % plot the final data
Z_SI = Z_final*(mu*1000);
Zvar_SI = Zerr_final*(mu*1000);
[rhoa,rhoaerr,phs,phserr] = calc_MT(Z_SI,Zvar_SI,f);
rhoa_xy = squeeze(rhoa(1,2,:));
rhoa_yx = squeeze(rhoa(2,1,:));
rhoa_xyerr = squeeze(rhoaerr(1,2,:));
rhoa_yxerr = squeeze(rhoaerr(2,1,:));
phs_xy = squeeze(phs(1,2,:));
phs_yx = squeeze(phs(2,1,:));
phs_xyerr = squeeze(phserr(1,2,:));
phs_yxerr = squeeze(phserr(2,1,:));
tx = squeeze(abs(T_final(1,1,:)));
ty = squeeze(abs(T_final(2,1,:)));
txvar = squeeze(Tvar_final(1,1,:));
tyvar = squeeze(Tvar_final(2,1,:));

subplot(12,1,1:3);
errorbar(f,rhoa_xy,rhoa_xyerr,'ro'); hold on
errorbar(f,rhoa_yx,rhoa_yxerr,'bo');
set(gca,'xdir','reverse','xscale','log','yscale','log');
ylabel('App.Res.')
grid on;
subplot(12,1,4:9);
errorbar(f,phs_xy,phs_xyerr,'ro'); hold on
errorbar(f,phs_yx,phs_yxerr,'bo');
set(gca,'xdir','reverse','xscale','log');
ylim([-180,180]);
grid on;
ylabel('Phase');
subplot(12,1,10:12);
errorbar(f,tx,txvar,'ro'); hold on
errorbar(f,ty,tyvar,'bo');
set(gca,'xdir','reverse','xscale','log','yscale','log');
xlabel('Frequency(Hz)');
ylabel('Tipper');
grid on;

% write out the final data
if exist([filepath,'merged.txt'],'file')    
    warning('Delete the exist merged.file before writing a new one!')
    delete([filepath,'merged.txt']);
end
fid = fopen([filepath,'merged.txt'],'a+');
for irow = 1:length(f)
    freq = f(irow);
    zxxr = real(squeeze(Z_final(1,1,irow))); zxxi = imag(squeeze(Z_final(1,1,irow)));
    zxyr = real(squeeze(Z_final(1,2,irow))); zxyi = imag(squeeze(Z_final(1,2,irow)));
    zyxr = real(squeeze(Z_final(2,1,irow))); zyxi = imag(squeeze(Z_final(2,1,irow)));
    zyyr = real(squeeze(Z_final(2,2,irow))); zyyi = imag(squeeze(Z_final(2,2,irow)));
    zxxerr = squeeze(Zerr_final(1,1,irow));
    zxyerr = squeeze(Zerr_final(1,2,irow));
    zyxerr = squeeze(Zerr_final(2,1,irow));
    zyyerr = squeeze(Zerr_final(2,2,irow));
    txr = real(squeeze(T_final(1,1,irow))); txi = imag(squeeze(T_final(1,1,irow)));
    tyr = real(squeeze(T_final(2,1,irow))); tyi = imag(squeeze(T_final(2,1,irow)));
    txvar = squeeze(Tvar_final(1,1,irow));  tyvar = squeeze(Tvar_final(2,1,irow));
    datarow = [freq,zxxr,zxxi,zxyr,zxyi,zyxr,zyxi,zyyr,zyyi,...
               zxxerr,zxyerr,zyxerr,zyyerr,...
               txr,txi,tyr,tyi,txvar,tyvar];
    fprintf(fid,'%12.10e\b\b',datarow);
    fprintf(fid,'\n');
end
fclose(fid);

    





























