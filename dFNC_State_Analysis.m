clc
clear all
close all

cl_path_gm='/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_Cluster_FDR/new/cl_map_sFC_HC_lrTLE_gm.nii';
cl_path_wm='/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_Cluster_FDR/new/cl_map_sFC_HC_lrTLE_wm.nii';
% cl_path_gm='/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_dFC_analysis/Kmeans/new_GM_VoxFC_clustering_K10.nii';
% cl_path_wm='/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_dFC_analysis/Kmeans/new_WM_VoxFC_clustering_K10.nii';

dp_hc_tle = '/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC_TLE';
dp_hc='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC';
dp_ltle='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/TLE_Left';
dp_rtle='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/TLE_Right';
output_folder = '/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_dFC_analysis/';

dgp_hc='/Users/sukeshdas/Documents/SKD/Data/ECP/Release_1.0_CTRL_Demographic_Behav_Data.txt';
dgt_hc=readtable(dgp_hc);
var_hc={'SubjID','Age','Sex'};
dgts_hc=dgt_hc(:,var_hc);
dgp_tle='/Users/sukeshdas/Documents/SKD/Data/ECP/Release_1.0_TLE_Demographic_Behav_Data.txt';
dgt_tle=readtable(dgp_tle); % Demographics table
vari={'SubjID','Sex','Age','TLEside','Aura_freq','SP_freq','CP_freq','SG_freq','OtherSz_freq','EEGSzcount'};
dgts_tle=dgt_tle(:,vari); % Demographics table selected var

addpath('/Users/sukeshdas/Documents/SKD/WM/');
fdn1='wsfREG_rc_';
fdn2='_task-rest_dir-AP_run-1_bold.nii';
fun_file='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC/sub-EC1006/wsfREG_rc_sub-EC1006_task-rest_dir-AP_run-1_bold.nii';
head=niftiinfo(fun_file);
mkdir(output_folder)

FDRc=1; % If FDR correction
RegOut=1;
FZT=1;
TR=0.802;
C=10;
Win=45;%30,90,50
Slide=1;
alpha=0.9;
beta=1;
lamda=0.5; %0.2
sfc=1;
State_Analysis=1;

cn_check=1;
cn=21;  %21,16 active voxel connectivity check

gm_cl_vol=double(niftiread(cl_path_gm));
[nx,ny,nz]=size(gm_cl_vol);
gm_cl_mat=reshape(gm_cl_vol,nx*ny*nz,1);
cl_gm=unique(sort(gm_cl_mat(:)));
cl_gm=cl_gm(2:end);
noc_gm=length(cl_gm);
cl_mat=zeros(nx*ny*nz,1);
for icl=1:1:noc_gm
    cl_mat(find(gm_cl_mat==cl_gm(icl)))=icl;
end

wm_cl_vol=double(niftiread(cl_path_wm));
wm_cl_mat=reshape(wm_cl_vol,nx*ny*nz,1);
cl_wm=unique(sort(wm_cl_mat(:)));
cl_wm=cl_wm(2:end);
noc_wm=length(cl_wm);
for icl=1:1:noc_wm
    cl_mat(find(wm_cl_mat==cl_wm(icl)))=noc_gm+icl;
end
cls=sort(unique(cl_mat(:)));
cls=cls(2:end);
noc=length(cls); % No of total clusters
cl_ids=cell(noc,1);
ncl_ids=zeros(noc,1);
for icl=1:1:noc
    ids=find(cl_mat==cls(icl));
    cl_ids{icl}=ids;
    ncl_ids(icl)=length(ids);
end

dir_hc_tle=dir(dp_hc_tle);
dir_hc_tle(1:3)=[];
nos1=length(dir_hc_tle); % No of initial subjects
dir_hc=dir(dp_hc);
dir_hc(1:3)=[];
dir_ltle=dir(dp_ltle);
dir_ltle(1:3)=[];
dir_rtle=dir(dp_rtle);
dir_rtle(1:3)=[];
disp('Checking Culsters in Subjects...')
reject_cl=[];
reject_sub=[];
for isub=1:1:nos1
    sub=dir_hc_tle(isub).name;
    disp(strcat('Reading Subject for Effective Clusters:: ',sub));
    func_data_name = strcat(fdn1,sub,fdn2);
    fdp=fullfile(dp_hc_tle,sub,func_data_name);
    func_data = double(niftiread(fdp));
    not=size(func_data,4);
    tvmat=reshape(func_data,nx*ny*nz,not);
    flag=1;
    for icl=1:1:noc
        tv_mat_id=tvmat(cl_ids{icl},:);
        tv_sum=std(tv_mat_id,0,2);
        zid=find(tv_sum==0);
        nid=find(isnan(tv_sum));
        reject_vid=[zid;nid];
        reject_vid=unique(reject_vid);
        if((ncl_ids(icl)-length(reject_vid))<3) % If Cluster does not have have altest 3 active voxel then rejet the cluster
            reject_cl=[reject_cl;icl];
            if(flag==1)
                reject_sub=[reject_sub;isub];
                flag=0;
            end
        end
    end
end
reject_cl=unique(sort(reject_cl));
cl_ids(reject_cl)=[];
noc=size(cl_ids,1); % New Clusters

fns_hc={dir_hc.name};
fns_ltle={dir_ltle.name};
fns_rtle={dir_rtle.name};
fns_tle=[fns_ltle fns_rtle];
fns_hc_tle={dir_hc_tle.name};
[~,re_hc]=ismember({dir_hc_tle(reject_sub).name},fns_hc);
[~,re_ltle]=ismember({dir_hc_tle(reject_sub).name},fns_ltle);
[~,re_rtle]=ismember({dir_hc_tle(reject_sub).name},fns_rtle);
re_hc(re_hc==0)=[];
re_ltle(re_ltle==0)=[];
re_rtle(re_rtle==0)=[];

dir_hc_tle(reject_sub)=[];
fns_hc(re_hc)=[];
fns_ltle(re_ltle)=[];
fns_rtle(re_rtle)=[];
nos=length(dir_hc_tle); % New no of subjects
[~,lhc]=ismember(fns_hc,fns_hc_tle);
[~,lltle]=ismember(fns_ltle,fns_hc_tle);
[~,lrtle]=ismember(fns_rtle,fns_hc_tle);
nhc=length(lhc);
nltle=length(lltle);
nrtle=length(lrtle);

subns_hc=char(fns_hc); % To get the subject names according to the names in the Behavioral data
subns_hc=(cellstr(subns_hc(:,5:end)))';
[~,hc2]=ismember(subns_hc,table2cell(dgts_hc(:,1))); % Maching the names of hc data considered with the Behavioral list
re_var_hc=table2cell(dgts_hc(hc2,["Sex" "Age"]));
subns_tle=char(fns_tle); % To get the subject names according to the names in the Behavioral data
subns_tle=(cellstr(subns_tle(:,5:end)))';
[~,ltle2]=ismember(subns_tle,table2cell(dgts_tle(:,1))); % Maching the names of tle data considered with the Behavioral list
re_var=table2cell(dgts_tle(ltle2,["Sex" "Age" "TLEside" "CP_freq" "EEGSzcount"])); % Required variables
disp(strcat('Left TLE Nos=',num2str(length(find(strcmp(re_var(:,3),'Left'))))));
disp(strcat('Right TLE Nos=',num2str(length(find(strcmp(re_var(:,3),'Right'))))));
disp(strcat('Bilateral TLE Nos=',num2str(length(find(strcmp(re_var(:,3),'Bilateral'))))));
disp(strcat('Uncertain TLE Nos=',num2str(length(find(strcmp(re_var(:,3),'Uncertain'))))));
disp(strcat('Male Nos=',num2str(length(find(strcmp(re_var(:,1),'Male'))))));
disp(strcat('Female Nos=',num2str(length(find(strcmp(re_var(:,1),'Female'))))));
disp(strcat('Avg Age=',num2str(mean(cell2mat(re_var(:,2))))));

idv=(1:noc).' >= (1:noc);
n=noc*(noc-1)/2;

disp(strcat('Sliding Window= ',num2str(Win)))
Gr=[];
disp(strcat('No of subjects = ',num2str(nos)));
hcs=0;
ltles=0;
rtles=0;

sFCv=zeros(n,nos);
dFCv_stat=[];

c_var_ltle=zeros(nltle,2);
c_var_rtle=zeros(nrtle,2);
c_sex_age=[];

T = table('Size',[nos 3],'VariableTypes',{'string','string','double'},'VariableNames',["Group","Sex","Age"]);
sub_frms=zeros(nos,1); % To keep the start and end index of the connectivity

for isub=1:1:nos
    sub=dir_hc_tle(isub).name;
    disp(strcat('Reading Subject for dFC computation:: ',sub));
    func_data_name = strcat(fdn1,sub,fdn2);
    fdp=fullfile(dp_hc_tle,sub,func_data_name);
    func_data = double(niftiread(fdp));
    not=size(func_data,4);
    tvmat=reshape(func_data,nx*ny*nz,not);
    sub_cl_tc=zeros(not,noc);
    for icl=1:1:noc
        tv_mat_id=tvmat(cl_ids{icl},:);
        tv_sum=std(tv_mat_id,0,2);
        zid=find(tv_sum==0);
        nid=find(isnan(tv_sum));
        reject_vid=[zid;nid];
        reject_vid=unique(reject_vid);
        tv_mat_id(reject_vid,:)=[];
        sub_cl_tc(:,icl)=(mean(tv_mat_id,1))';
    end
    if(sfc==1)
        sFCmat=corr(sub_cl_tc);

        if(FZT==1)
            sFCv(:,isub)=atanh(sFCmat(idv==0));
        else
            sFCv(:,isub)=sFCmat(idv==0);
        end
    end
    nframes=floor((not-Win)/Slide)+1;
    sub_frms(isub)=nframes;
    dFCmat=zeros(noc,noc,nframes);

    ifr1=0;
    for i=1:Slide:(not-Win)
        ifr1=ifr1+1;
        tc_frm = sub_cl_tc(i:i+Win-1,:);
        if(~isempty(find(isnan(corr(tc_frm)))))
            tc_frm(tc_frm==0)=eps;
        end
        corr_mat=corr(tc_frm);%tanh(inv(L1precisionBCD(cov(tc_frm),lamda)));%corr(tc_frm);
        dFCmat(:,:,ifr1)=corr_mat;
    end
    dFCv=reshape(dFCmat,noc*noc,nframes);
    dFCv=dFCv(idv==0,:);
    if(FZT==1)
        dFCvFT=atanh(dFCv);
    else
        dFCvFT=dFCv;
    end
    dFCv_stat=[dFCv_stat;dFCvFT'];

    if(find(isub==lhc))
        grp=0;
        hcs=hcs+1;
        T.Group(isub)='HC';
        T.Sex(isub)=cell2mat(re_var_hc(hcs,1));
        T.Age(isub)=cell2mat(re_var_hc(hcs,2));
        if(strcmp(cell2mat(re_var_hc(hcs,1)),'Male')==1)
            sex=1;
        elseif(strcmp(cell2mat(re_var_hc(hcs,1)),'Female')==1)
            sex=-1;
        end
        c_sex_age=[c_sex_age;repmat([sex,cell2mat(re_var_hc(hcs,2))],nframes,1)];
    elseif(find(isub==lltle))
        grp=-1;
        ltles=ltles+1;
        c_var_ltle(ltles,:)=[cell2mat(re_var(ltles+rtles,4)),str2num(cell2mat(re_var(ltles+rtles,5)))];
        T.Group(isub)='lTLE';
        T.Sex(isub)=cell2mat(re_var(ltles+rtles,1));
        T.Age(isub)=cell2mat(re_var(ltles+rtles,2));
        if(strcmp(cell2mat(re_var(ltles+rtles,1)),'Male')==1)
            sex=1;
        elseif(strcmp(cell2mat(re_var(ltles+rtles,1)),'Female')==1)
            sex=-1;
        end
        c_sex_age=[c_sex_age;repmat([sex,cell2mat(re_var(ltles+rtles,2))],nframes,1)];
    elseif(find(isub==lrtle))
        grp=1;
        rtles=rtles+1;
        c_var_rtle(rtles,:)=[cell2mat(re_var(ltles+rtles,4)),str2num(cell2mat(re_var(ltles+rtles,5)))];
        T.Group(isub)='rTLE';
        T.Sex(isub)=cell2mat(re_var(ltles+rtles,1));
        T.Age(isub)=cell2mat(re_var(ltles+rtles,2));
        if(strcmp(cell2mat(re_var(ltles+rtles,1)),'Male')==1)
            sex=1;
        elseif(strcmp(cell2mat(re_var(ltles+rtles,1)),'Female')==1)
            sex=-1;
        end
        c_sex_age=[c_sex_age;repmat([sex,cell2mat(re_var(ltles+rtles,2))],nframes,1)];
    end

    Gr=[Gr;repmat(grp,nframes,1)];
end
Sub_G=Gr(cumsum(sub_frms));
if(RegOut==1)
    for j=1:1:n
        [~,~,sFCv(j,:)]=regress((sFCv(j,:))',[c_sex_age(cumsum(sub_frms),:),ones(nos,1)]);
        [~,~,dFCv_stat(:,j)]=regress((dFCv_stat(:,j)),[c_sex_age,ones(sum(sub_frms),1)]);
    end
end

sFCv_hc=sFCv(:,find(Sub_G==0));
sFCv_ltle=sFCv(:,find(Sub_G==-1));
sFCv_rtle=sFCv(:,find(Sub_G==1));
sFCv_tle=[sFCv_ltle,sFCv_rtle];
asFCv=mean([sFCv_hc,sFCv_tle],2);

dFCv_stat_hc=dFCv_stat(find(Gr==0),:);
dFCv_stat_tle=dFCv_stat(find(Gr==-1 | Gr==1),:);
% dFCv_stat_rtle=dFCv_stat(find(Gr==1),:);
% dFCv_stat_tle=[dFCv_stat_ltle;dFCv_stat_rtle];

if(State_Analysis==1)
    dm='cosine';%'sqeuclidean';% 'cosine','correlation'
    num_replicates=1;%50;

    K=2:15;
    Elb=zeros(length(K),2);
    for ki=1:1:length(K)
        k=K(ki);
        disp(['K = ' num2str(k)]);
        % load('wm_atlas1_w30_seed_8.mat');
        % seed=zeros(k,size(all_subs_dFCmat,2));
        % for j=1:1:k
        %     id11=find(id1==j);
        %     seed(j,:)=mean(all_subs_dFCmat(id11,:),1);
        % end
        % [id1,c1,sd]=kmeans(all_subs_dFCmat,k,'Distance',dm,'start',seed); %Time point clustering
        [id1,c1,sd]=kmeans(dFCv_stat,k,'Distance',dm,'replicates',num_replicates); %Time point clustering
        inter_cld=0;
        for jkr=1:1:k-1
            for jkc=jkr+1:1:k
                inter_cld=inter_cld+sum((c1(jkr,:)./sqrt(sum(c1(jkr,:).^2))-c1(jkc,:)./sqrt(sum(c1(jkc,:).^2))).^2);
            end
        end
        Elb(ki,:)=[k sum(sd.^2)/inter_cld];
        % save(strcat('wm_atlas1_w30_seed_',num2str(k)),'id1');
    end
    figure(1)
    plot(Elb(:,1),Elb(:,2),'Color',[0 0 1],'Linewidth',2)
    hold on
    plot(Elb(:,1),Elb(:,2),'*','Color',[1 0 1],'Linewidth',2)
    ylabel('Distortion')
    xlabel('ks','fontsize',11,'Interpreter','latex')
    title('Elbow Method')
    % saveas(gcf,fullfile(output_folder,'K_elbow_Vox_dFNC_GM_WM.png'))

    K = str2double(cell2mat(inputdlg('Choose the K-value','K-Value')));
    close gcf
    [id,cl,seed] = kmeans(dFCv_stat,K,'distance',dm,'replicates',2*num_replicates);            % K-means clustering
    % save([output_folder,'ROI_CL_dFNC_W45_st_id_1'],'id');
    % save([output_folder,'ROI_CL_dFNC_W45_st_centers_1'],'cl');
    % save([output_folder,'ROI_CL_dFNC_W90_Lpre_L02_st_id'],'id');
    % save([output_folder,'ROI_CL_dFNC_W90_Lpre_L02_st_centers'],'cl');

    load([output_folder,'ROI_CL_dFNC_W45_st_id_1']);
    load([output_folder,'ROI_CL_dFNC_W45_st_centers_1']);

    id_hc=id(find(Gr==0));
    id_tle=id(find(Gr==-1 | Gr==1));
    id_ltle=id(find(Gr==-1));
    id_rtle=id(find(Gr==1));
    st_tf_gr=zeros(K,5); % probability of time frames in every states and for 4 groups (total probability and prob of hc,tle,ltle,rtle)
    states=zeros(noc,noc,K);
    % states_hc=zeros(noc,noc,K); % cluster means of states of 4 groups(HC,TLE,lTLE,rTLE)
    % states_tle=zeros(noc,noc,K);
    for ik=1:1:K
        state=zeros(noc,noc);
        state(find(idv==0)) = mean(dFCv_stat(find(id==ik),:),1);%cl(ik,:);
        state=state';
        state(find(idv==0))= mean(dFCv_stat(find(id==ik),:),1);%cl(ik,:);
        states(:,:,ik)=state;

        % state_hc=zeros(noc,noc);
        % state_hc(find(idv==0))=mean(dFCv_stat_hc(find(id_hc==ik),:)./repmat(sqrt(sum(dFCv_stat_hc(find(id_hc==ik),:).^2,2)),1,n),1);
        % state_hc=state_hc';
        % state_hc(find(idv==0))=mean(dFCv_stat_hc(find(id_hc==ik),:)./repmat(sqrt(sum(dFCv_stat_hc(find(id_hc==ik),:).^2,2)),1,n),1);
        % state_tle=zeros(noc,noc);
        % state_tle(find(idv==0))=mean(dFCv_stat_tle(find(id_tle==ik),:)./repmat(sqrt(sum(dFCv_stat_tle(find(id_tle==ik),:).^2,2)),1,n),1);
        % state_tle=state_tle';
        % state_tle(find(idv==0))=mean(dFCv_stat_tle(find(id_tle==ik),:)./repmat(sqrt(sum(dFCv_stat_tle(find(id_tle==ik),:).^2,2)),1,n),1);
        %
        % states_hc(:,:,ik)=state_hc;
        % states_tle(:,:,ik)=state_tle;

        nf_st_hc=length(find(id_hc==ik));
        nf_st_ltle=length(find(id_ltle==ik));
        nf_st_rtle=length(find(id_rtle==ik));
        nf_st=nf_st_hc+nf_st_ltle+nf_st_rtle;

        st_tf_gr(ik,1)=nf_st/length(Gr);
        st_tf_gr(ik,2)=nf_st_hc/length(id_hc);
        st_tf_gr(ik,3)=(nf_st_ltle+nf_st_rtle)/(length(id_ltle)+length(id_rtle));
        st_tf_gr(ik,4)=nf_st_ltle/length(id_ltle);
        st_tf_gr(ik,5)=nf_st_rtle/length(id_rtle);
    end

    lbls_gm="GM" +string(1:noc_gm);
    lbls_wm="WM" +string(1:noc_wm);
    lbls=[lbls_gm,lbls_wm];

    figure(2);
    subplot(2,3,1)
    h1=heatmap(states(:,:,1),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(a) State 1')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(2,3,2)
    h2=heatmap(states(:,:,2),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(b) State 2')
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(2,3,3)
    h3=heatmap(states(:,:,3),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(c) State 3')
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;
    subplot(2,3,4)
    h1=heatmap(states(:,:,4),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(d) State 4')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    xlabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(2,3,5)
    h2=heatmap(states(:,:,5),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(e) State 5')
    xlabel('FNs');
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(2,3,6)
    h3=heatmap(states(:,:,6),'ColorLimits',[min(states(:)) max(states(:))]);
    title('(f) State 6')
    xlabel('FNs');
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dFC_state_sFC_corr=zeros(K,1);
    for istate=1:1:K
        x=states(:,:,istate);
        sig=x(idv==0);
        dFC_state_sFC_corr(istate)=corr(asFCv,sig);
        mx=max(sig);
        mn=min(sig);
        levels_sig=linspace(mn,mx,50);
        [n_sig bin_sig]=histc(sig,levels_sig);
        bl_sig=length(n_sig);
        sl_r=0.99;
        sl_rv=round(length(sig)*sl_r);
        cdf1=cumsum(n_sig);
        rli=sum(cdf1<sl_rv)+1;
        rl=levels_sig(rli(1));
        n_sig_l=flip(n_sig);
        levels_sig_l=flip(levels_sig);
        cdf2=cumsum(n_sig_l);
        lli=sum(cdf2<sl_rv)+1;
        ll=levels_sig_l(lli(1));
        C=[];
        for i=1:1:size(x,1)-1
            for j=i+1:1:size(x,2)
                if(x(i,j)>=rl || x(i,j)<=ll)
                    C=[C;[i,j,x(i,j)]];
                end
            end
        end
        Cs=flip(C(:,3));
        C(:,3)=abs(C(:,3));
        tcn=size(C,1);
        C(:,3)=C(:,3)*10;
        CT = array2table(C);
        groups=1:1:size(x,1);
        p=chordPlot(groups,CT);
        h = get(gca,'Children');
        % set(h(1:tcn), 'Color', 'r');
        for icn=1:1:tcn
            if(Cs(icn)<0)
                set(h(icn), 'Color', 'b');
            else
                set(h(icn), 'Color', 'r');
            end
        end
        title(strcat('State ',num2str(istate)))
        p.NodeLabel=lbls;
    end

    figure
    xtl=cellstr("State" +string(1:K));
    bar(find(dFC_state_sFC_corr>=0),dFC_state_sFC_corr(dFC_state_sFC_corr>=0),'FaceColor','r')
    hold on
    bar(find(dFC_state_sFC_corr<0),dFC_state_sFC_corr(dFC_state_sFC_corr<0),'FaceColor','b')
    hold off
    title('Correlation with sFC')
    set(gca,'Xtick',1:K,'XTickLabel',xtl)
    ylabel('Correlation')

    figure
    bp=bar(st_tf_gr,'FaceColor','flat');
    set(bp, {'DisplayName'}, {'Total probability','Probability HC','Probability TLE','Probability lTLE','Probability rTLE'}')
    legend();
    xticklabels(xtl)
    ylabel('Probability of time frames of groups')

    stfrm=0;
    scp_hc=zeros(nhc,1);
    scp_ltle=zeros(nltle,1);
    scp_rtle=zeros(nrtle,1);
    stm_hc=zeros(K,K,nhc);%State transition mat
    stm_ltle=zeros(K,K,nltle);
    stm_rtle=zeros(K,K,nrtle);
    st_sub_gr=zeros(3,K); % subject in each state
    % st_gr=zeros(K,n,3); % Groupwise states
    cen_sub_states=zeros(n,nos,K); % Centeroid of every subject and every state (every column is a centeroid of every sub and 3 dim is for states )
    hcs=0;ltles=0;rtles=0;
    for isub=1:1:nos
        nofs=sub_frms(isub);
        sub_id=id(stfrm+1:stfrm+nofs);
        dFCv_sub=dFCv_stat(stfrm+1:stfrm+nofs,:);
        % Dual regression for estimating subject-wise states
        beta1=zeros(nofs,K);
        beta2=zeros(K,n);
        for ir1=1:1:nofs
            beta1(ir1,:)=regress((dFCv_sub(ir1,:))',cl'); %[cl',ones(n,1)]
        end
        for ir2=1:1:n
            beta2(:,ir2)=regress(dFCv_sub(:,ir2),beta1);%[beta1,ones(nofs,1)]
        end

        st_gr_tem=zeros(K,n);
        for ik=1:1:K
            % kid=find(sub_id==ik);
            % cen_sub_states(:,isub,ik)=mean(dFCv_sub(kid,:)./repmat(sqrt(sum(dFCv_sub(kid,:).^2,2)),1,n),1);
            cen_sub_states(:,isub,ik)=beta2(ik,:);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stfrm=stfrm+nofs;
        pc=length(find(diff(sub_id)~=0))/(nofs-1); % Probability of changes of states
        gr_id=Gr(stfrm);
        stm=zeros(K,K); % temporary state transition matrix
        sub_id_p=sub_id(1:end-1); % previous frame's state ids
        sub_id_n=sub_id(2:end);% next frame's state ids
        for istate1=1:1:length(sub_id_p)
            p_st=sub_id_p(istate1);
            n_st=sub_id_n(istate1);
            stm(p_st,n_st)=stm(p_st,n_st)+1;
        end
        stm=stm./length(sub_id_p); % conditional Probability
        for istate=1:1:K
            stm(istate,:)=stm(istate,:)./(length(find(sub_id_p==istate))./length(sub_id_p)); % Multiplying prior with conditional prob
        end
        stm(isnan(stm))=0;
        sids=unique(sort(sub_id)); % current subject associated with which clusters/states
        if(gr_id==0)
            hcs=hcs+1;
            scp_hc(hcs)=pc;
            stm_hc(:,:,hcs)=stm;
            st_sub_gr(1,sids)=st_sub_gr(1,sids)+1;
            % st_gr(:,:,1)=st_gr(:,:,1)+st_gr_tem;
        elseif(gr_id==-1)
            ltles=ltles+1;
            scp_ltle(ltles)=pc;
            stm_ltle(:,:,ltles)=stm;
            st_sub_gr(2,sids)=st_sub_gr(2,sids)+1;
            % st_gr(:,:,2)=st_gr(:,:,2)+st_gr_tem;
        elseif(gr_id==1)
            rtles=rtles+1;
            scp_rtle(rtles)=pc;
            stm_rtle(:,:,rtles)=stm;
            st_sub_gr(3,sids)=st_sub_gr(3,sids)+1;
            % st_gr(:,:,3)=st_gr(:,:,3)+st_gr_tem;
        end
    end
    stm_tle=cat(3,stm_ltle,stm_rtle);
    astm_hc=mean(stm_hc,3);
    astm_tle=mean(stm_tle,3);
    astm_ltle=mean(stm_ltle,3);
    astm_rtle=mean(stm_rtle,3);

    [~,p_hc_tle_stm,~,~]=ttest2(reshape(stm_hc,nhc,K*K),reshape(stm_tle,(nltle+nrtle),K*K));
    if(FDRc==1)
        cp_hc_tle_stm=mafdr(p_hc_tle_stm,'BHFDR','true');
    else
        cp_hc_tle_stm=p_hc_tle_stm;
    end
    cp_hc_tle_stm=reshape(cp_hc_tle_stm,K,K);

    cen_sub_states_hc=cen_sub_states(:,find(Sub_G==0),:);
    cen_sub_states_tle=cen_sub_states(:,find(Sub_G==-1 | Sub_G==1),:);
    hc_dFNC_states=zeros(noc,noc,K);
    tle_dFNC_states=zeros(noc,noc,K);
    for ik=1:1:K
        hc_dFNC_state=zeros(noc,noc);
        hc_dFNC_state(find(idv==0))= mean(cen_sub_states_hc(:,:,ik),2);
        hc_dFNC_state=hc_dFNC_state';
        hc_dFNC_state(find(idv==0))= mean(cen_sub_states_hc(:,:,ik),2);
        hc_dFNC_states(:,:,ik)=hc_dFNC_state;

        tle_dFNC_state=zeros(noc,noc);
        tle_dFNC_state(find(idv==0))= mean(cen_sub_states_tle(:,:,ik),2);
        tle_dFNC_state=tle_dFNC_state';
        tle_dFNC_state(find(idv==0))= mean(cen_sub_states_tle(:,:,ik),2);
        tle_dFNC_states(:,:,ik)=tle_dFNC_state;
    end
    min_level=min([hc_dFNC_states(:);tle_dFNC_states(:)]);
    max_level=max([hc_dFNC_states(:);tle_dFNC_states(:)]);
    figure
    subplot(4,3,1)
    h1=heatmap(hc_dFNC_states(:,:,1),'ColorLimits',[min_level max_level]);
    title('(a) State 1 (HC)')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(4,3,2)
    h2=heatmap(hc_dFNC_states(:,:,2),'ColorLimits',[min_level max_level]);
    title('(b) State 2 (HC)')
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(4,3,3)
    h3=heatmap(hc_dFNC_states(:,:,3),'ColorLimits',[min_level max_level]);
    title('(c) State 3 (HC)')
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;
    subplot(4,3,4)
    h1=heatmap(tle_dFNC_states(:,:,1),'ColorLimits',[min_level max_level]);
    title('(d) State 1 (TLE)')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    % xlabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(4,3,5)
    h2=heatmap(tle_dFNC_states(:,:,2),'ColorLimits',[min_level max_level]);
    title('(e) State 2 (TLE)')
    % xlabel('FNs');
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(4,3,6)
    h3=heatmap(tle_dFNC_states(:,:,3),'ColorLimits',[min_level max_level]);
    title('(f) State 3 (TLE)')
    % xlabel('FNs');
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;
    subplot(4,3,7)
    h1=heatmap(hc_dFNC_states(:,:,4),'ColorLimits',[min_level max_level]);
    title('(g) State 4 (HC)')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(4,3,8)
    h2=heatmap(hc_dFNC_states(:,:,5),'ColorLimits',[min_level max_level]);
    title('(h) State 5 (HC)')
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(4,3,9)
    h3=heatmap(hc_dFNC_states(:,:,6),'ColorLimits',[min_level max_level]);
    title('(i) State 6 (HC)')
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;
    subplot(4,3,10)
    h1=heatmap(tle_dFNC_states(:,:,4),'ColorLimits',[min_level max_level]);
    title('(j) State 4 (TLE)')
    colorbar;
    colormap(jet)
    ylabel('FNs');
    xlabel('FNs');
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    subplot(4,3,11)
    h2=heatmap(tle_dFNC_states(:,:,5),'ColorLimits',[min_level max_level]);
    title('(k) State 5(TLE)')
    xlabel('FNs');
    colorbar;
    colormap(jet)
    h2.XDisplayLabels = lbls;
    h2.YDisplayLabels = lbls;
    subplot(4,3,12)
    h3=heatmap(tle_dFNC_states(:,:,6),'ColorLimits',[min_level max_level]);
    title('(l) State 6 (TLE)')
    xlabel('FNs');
    colorbar;
    colormap(jet)
    h3.XDisplayLabels = lbls;
    h3.YDisplayLabels = lbls;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(1,3,1)
    bar([1:1:nhc],scp_hc)
    ylim([0 max([scp_hc;scp_ltle;scp_rtle])])
    hold on
    line([0 nhc],[mean(scp_hc) mean(scp_hc)],'color','m','Linewidth',3);
    title("(a) HC")
    xlabel('Subjects')
    ylabel('Probability of state transition')
    subplot(1,3,2)
    bar([1:1:nltle],scp_ltle)
    ylim([0 max([scp_hc;scp_ltle;scp_rtle])])
    hold on
    line([0 nltle],[mean(scp_ltle) mean(scp_ltle)],'color','m','Linewidth',3);
    title("(b) lTLE")
    xlabel('Subjects')
    subplot(1,3,3)
    bar([1:1:nrtle],scp_rtle)
    ylim([0 max([scp_hc;scp_ltle;scp_rtle])])
    hold on
    line([0 nrtle],[mean(scp_rtle) mean(scp_rtle)],'color','m','Linewidth',3);
    title("(c) rTLE")
    xlabel('Subjects')

    figure
    bp=bar(st_sub_gr','FaceColor','flat');
    set(bp, {'DisplayName'}, {'HC','lTLE','rTLE'}')
    legend();
    xticklabels(xtl)
    ylabel('No of subjects')

    xtl1="State" +string(1:K);
    
    figure
    subplot(1,3,1)
    h11=heatmap(round(astm_hc,3),'ColorLimits',[0 .7]);
    colorbar;
    title(strcat('(a) Average state transition (HC)'))
    h11.XDisplayLabels = xtl1;
    h11.YDisplayLabels = xtl1;
    xlabel('To state')
    ylabel('From state')

    subplot(1,3,2)
    h12=heatmap(round(astm_ltle,3),'ColorLimits',[0 .7]);
    colorbar;
    title(strcat('(b) Average state transition (lTLE)'))
    h12.XDisplayLabels = xtl1;
    h12.YDisplayLabels = xtl1;
    xlabel('To state')
    ylabel('From state')

    subplot(1,3,3)
    h13=heatmap(round(astm_rtle,3),'ColorLimits',[0 .7]);
    colorbar;
    title(strcat('(c) Average state transition (rTLE)'))
    h13.XDisplayLabels = xtl1;
    h13.YDisplayLabels = xtl1;
    xlabel('To state')
    ylabel('From state')

    % cen_sub_states=zeros(n,nos,K);
    % for isub=1:1:nos
    %     frms=sub_frms(isub);
    %
    %
    % end

    for ik=1:1:K
        [~,p_hc_tle,~,stats_hc_tle]=ttest2((cen_sub_states_hc(:,:,ik))',(cen_sub_states_tle(:,:,ik))');
        t_hc_tle= stats_hc_tle.tstat;

        % [~,p_hc_ltle,~,stats_hc_ltle]=ttest2(dFCv_stat_hc(find(id_hc==ik),:),dFCv_stat_ltle(find(id_ltle==ik),:));
        % t_hc_ltle= stats_hc_ltle.tstat;
        % [~,p_hc_rtle,~,stats_hc_rtle]=ttest2(dFCv_stat_hc(find(id_hc==ik),:),dFCv_stat_rtle(find(id_rtle==ik),:));
        % t_hc_rtle= stats_hc_rtle.tstat;
        if(FDRc==1)
            cp_hc_tle=mafdr(p_hc_tle,'BHFDR','true');
            % cp_hc_ltle=mafdr(p_hc_ltle,'BHFDR','true');
            % cp_hc_rtle=mafdr(p_hc_rtle,'BHFDR','true');

        else
            cp_hc_tle=p_hc_tle;
            % cp_hc_ltle=p_hc_ltle;
            % cp_hc_rtle=p_hc_rtle;
        end
        hc_dFNC_mat=zeros(noc,noc);
        % hc_dFNC_mat(find(idv==0))= mean(dFCv_stat_hc(find(id_hc==ik),:),1);
        hc_dFNC_mat(find(idv==0))= mean(cen_sub_states_hc(:,:,ik),2);
        hc_dFNC_mat=hc_dFNC_mat';
        % hc_dFNC_mat(find(idv==0))=mean(dFCv_stat_hc(find(id_hc==ik),:),1);
        hc_dFNC_mat(find(idv==0))= mean(cen_sub_states_hc(:,:,ik),2);
        %
        tle_dFNC_mat=zeros(noc,noc);
        % tle_dFNC_mat(find(idv==0))=mean(dFCv_stat_tle(find(id_tle==ik),:),1);
        tle_dFNC_mat(find(idv==0))= mean(cen_sub_states_tle(:,:,ik),2);
        tle_dFNC_mat=tle_dFNC_mat';
        % tle_dFNC_mat(find(idv==0))=mean(dFCv_stat_tle(find(id_tle==ik),:),1);
        tle_dFNC_mat(find(idv==0))= mean(cen_sub_states_tle(:,:,ik),2);


        cp_hc_tle_mat=zeros(noc,noc);
        cp_hc_tle_mat(find(idv==0))= cp_hc_tle;
        cp_hc_tle_mat=cp_hc_tle_mat';
        cp_hc_tle_mat(find(idv==0))= cp_hc_tle;

        p_min=0.05;

        figure
        % subplot(1,3,1)
        % imagesc(hc_dFNC_mat);
        % colorbar;
        % colormap(jet)
        % ylabel('FNs');
        % title(strcat('(a) State ',num2str(ik), ', dFNC (HC)'))
        subplot(1,3,1)
        h1=heatmap(hc_dFNC_mat,'ColorLimits',[min([hc_dFNC_mat(:);tle_dFNC_mat(:)]) max([hc_dFNC_mat(:);tle_dFNC_mat(:)])]);
        colorbar;
        colormap(jet)
        ylabel('FNs');
        h1.XDisplayLabels = lbls;
        h1.YDisplayLabels = lbls;
        title(strcat('(a) State ',num2str(ik), ', dFNC (HC)'))

        subplot(1,3,2)
        h2=heatmap(tle_dFNC_mat,'ColorLimits',[min([hc_dFNC_mat(:);tle_dFNC_mat(:)]) max([hc_dFNC_mat(:);tle_dFNC_mat(:)])]);
        colorbar;
        colormap(jet)
        ylabel('FNs');
        h2.XDisplayLabels = lbls;
        h2.YDisplayLabels = lbls;
        title(strcat('(b) State ',num2str(ik), ', dFNC (TLE)'))

        subplot(1,3,3)
        h3=heatmap(cp_hc_tle_mat,'ColorLimits',[0 p_min]);
        colorbar;
        colormap(jet)
        ylabel('FNs');
        h3.XDisplayLabels = lbls;
        h3.YDisplayLabels = lbls;
        title(strcat('(c) State ',num2str(ik), ', p value'))

        x=tle_dFNC_mat-hc_dFNC_mat;
        C=[];
        for i=1:1:size(x,1)-1
            for j=i+1:1:size(x,2)
                if(cp_hc_tle_mat(i,j)<p_min)
                    C=[C;[i,j,x(i,j)]];
                end
            end
        end
        if(~isempty(C))
            Cs=flip(C(:,3));
            C(:,3)=abs(C(:,3));
            tcn=size(C,1);
            C(:,3)=C(:,3)*10;
            CT = array2table(C);
            groups=1:1:size(x,1);
            p=chordPlot(groups,CT);
            h = get(gca,'Children');
            % set(h(1:tcn), 'Color', 'r');
            for icn=1:1:tcn
                if(Cs(icn)<0)
                    set(h(icn), 'Color', 'b');
                else
                    set(h(icn), 'Color', 'r');
                end
            end
            title(strcat('State ',num2str(ik)))
            p.NodeLabel=lbls;
        end
    end
end