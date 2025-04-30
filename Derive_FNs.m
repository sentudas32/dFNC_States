clc
clear all
close all
GM=1;
WM=0;
FDR=1;
AdMat=1; % using n fold cross validation of adjacency matrices
Elbow=0;

atlasp_gm='/Users/sukeshdas/Documents/SKD/Atlas/atlas400_Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii';
atlasp_wm='/Users/sukeshdas/Documents/SKD/Atlas/Eve_Atlas-master/JHU_MNI_SS_WMPM_Type-III_to_MNI_brain.nii.gz';
atlas_xlsp='/Users/sukeshdas/Documents/SKD/Atlas/atlas400_Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.xlsx';
dp_hc_tle = '/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC_TLE';
dp_hc='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC';
dp_ltle='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/TLE_Left';
dp_rtle='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/TLE_Right';
output_folder = '/Users/sukeshdas/Documents/SKD/ECP_Exp/OutPut_ECP_SmoothSeparate/sFC_Cluster_FDR/new/';
addpath('/Users/sukeshdas/Documents/SKD/WM/');
fdn1='wsfREG_rc_';
fdn2='_task-rest_dir-AP_run-1_bold.nii';
fun_file='/Users/sukeshdas/Documents/SKD/Data/ECP/ECP_SmoothSeparate/func/HC/sub-EC1006/wsfREG_rc_sub-EC1006_task-rest_dir-AP_run-1_bold.nii';
mkdir(output_folder)

if(GM==1)
    atlasp=atlasp_gm;
    reject_roi=[];%[38,39,40,41,42,43,81,100,127,178,179,180,236,238,239,240,241,242,243,272,280,281,284,299,302,303,380,381];
elseif(WM==1)
    atlasp=atlasp_wm;
    reject_roi=[];%[17,29,30,48,49,62,63,89,90,108,109,118];
end

atlas_vol=double(niftiread(atlasp));
[nx1,ny1,nz1]=size(atlas_vol);

atlas_rs=reslice_data(atlasp,fun_file,0);
[nx,ny,nz]=size(atlas_rs);
atlas_rs_mat=reshape(atlas_rs,nx*ny*nz,1);
roi_select=unique(sort(atlas_rs_mat(:)));
roi_select=roi_select(2:end);
roi_select(reject_roi)=[];

nroi=length(roi_select);
roi_ids=cell(nroi,1);
nroi_ids=zeros(nroi,1);
for iroi=1:1:nroi
    ids=find(atlas_rs_mat==roi_select(iroi));
    roi_ids{iroi}=ids;
    nroi_ids(iroi)=length(ids);
end

dir_hc_tle=dir(dp_hc_tle);
dir_hc_tle(1:3)=[];
nos=length(dir_hc_tle); % No of initial subjects
dir_hc=dir(dp_hc);
dir_hc(1:3)=[];
dir_ltle=dir(dp_ltle);
dir_ltle(1:3)=[];
dir_rtle=dir(dp_rtle);
dir_rtle(1:3)=[];

fns_hc={dir_hc.name};
fns_ltle={dir_ltle.name};
fns_rtle={dir_rtle.name};
fns_hc_tle={dir_hc_tle.name};
[~,lhc]=ismember(fns_hc,fns_hc_tle);
[~,lltle]=ismember(fns_ltle,fns_hc_tle);
[~,lrtle]=ismember(fns_rtle,fns_hc_tle);

G=[];
disp(strcat('No of subjects = ',num2str(nos)));
esub=[]; % Subs to be excluded
esub1=[];
reject_roi1=[];
hcs=0;
ltles=0;
rtles=0;

sFC_hc=zeros(nroi,nroi,length(lhc));
sFC_ltle=zeros(nroi,nroi,length(lltle));
sFC_rtle=zeros(nroi,nroi,length(lrtle));

for isub=1:1:nos
    sub=dir_hc_tle(isub).name;
    disp(strcat('Reading Subject....',sub));
    func_data_name = strcat(fdn1,sub,fdn2);
    fdp=fullfile(dp_hc_tle,sub,func_data_name);
    func_data = double(niftiread(fdp));
    not=size(func_data,4);
    tvmat=reshape(func_data,nx*ny*nz,not);
    % st_id=0;
    sub_roi_tc=zeros(not,nroi);
    for iroi=1:1:nroi
        % nv=nroi_ids(iroi);
        ids=cell2mat(roi_ids(iroi,1));
        tv_mat_id=tvmat(ids,:);
        tv_sum=std(tv_mat_id,0,2);
        zid=find(tv_sum==0);
        nid=find(isnan(tv_sum));
        reject_vid=sort(unique([zid;nid]));
        tv_mat_id(reject_vid,:)=[];
        if(isempty(tv_mat_id))
            reject_roi1=[reject_roi1 iroi];
            % break;
        end
        tv_mean=mean(tv_mat_id);
        sub_roi_tc(:,iroi)=tv_mean;
        % st_id=st_id+nv;
    end
    if(iroi==nroi)
        sFC=corr(sub_roi_tc);
        nid1=find(isnan(sFC));
        if(~isempty(nid1))
           sFC(nid1)=eps;
        end
        if(find(isub==lhc))
            grp=0;
            hcs=hcs+1;
            sFC_hc(:,:,hcs)=sFC;
        elseif(find(isub==lltle))
            grp=1;
            ltles=ltles+1;
            sFC_ltle(:,:,ltles)=sFC;
        elseif(find(isub==lrtle))
            grp=1;
            rtles=rtles+1;
            sFC_rtle(:,:,rtles)=sFC;
        end
        G=[G;grp];
        if(find(isnan(sFC(:))))
            esub1=[esub1,isub];
        end

    else
        esub=[esub,isub];
        
    end
    
end
sFC_hc_tle=cat(3,sFC_hc,sFC_ltle,sFC_rtle);
asFC_hc_tle=mean(sFC_hc_tle,3);
clear sFC_hc;clear sFC_ltle;clear sFC_rtle;
% asFC_ltle=mean(sFC_ltle,3);
% asFC_rtle=mean(sFC_rtle,3);
% asFC_tle=mean(sFC_tle,3);

if(AdMat==1)
    disp('8 Fold Cross validation in process to find out the correct value of K')
    % separating the data into subsets of features
    num_CV_folds = 8;
    num_replicates = 30;
    IDX_folds_new = cell(1,num_CV_folds);
    K_range_l=5;
    K_range_h=20;

    for K= K_range_l:K_range_h          % going over all possible numbers of clusters, to measure each one's stability
        disp(['K = ' num2str(K)]);
        IDX_folds = zeros(size(asFC_hc_tle,1),num_CV_folds); IDX_folds_new{K} = zeros(size(IDX_folds));
        size_fold = size(asFC_hc_tle,2)/num_CV_folds;
        for c=1:num_CV_folds        % going over folds (sub-matrices)
            disp(c)
            mat_corr_current = asFC_hc_tle(:,round((c-1)*size_fold+1):round(c*size_fold));     % the sub-correlation-matrix
            IDX_folds(:,c) = kmeans(mat_corr_current, K,'distance','correlation','replicates',num_replicates);  % calculating the clustering result for this K
        end

        % computing the difference between adjacency matrices for each fold
        size_chunk = 30;   % need to compute adjacency matrices in parts - otherwise it takes too much memory (~300 million numbers per matrix)
        num_chunks = floor(size(IDX_folds,1) / size_chunk);
        sum_diff_adjmats_folds = zeros(num_CV_folds); sum_common_connections_adjmats = zeros(num_CV_folds); sum_all_connections_adjmats = zeros(num_CV_folds);
        for ch1=1:num_chunks
            for ch2=1:num_chunks    % iterating over all adjacency matrix parts combinations
                current_clustering_adjmats = zeros(size_chunk,size_chunk,num_CV_folds);
                current_chunk1 = (ch1-1)*size_chunk+1 : ch1*size_chunk;
                current_chunk2 = (ch2-1)*size_chunk+1 : ch2*size_chunk;

                 % creating the current adjacency matrix part, for all folds
                for c=1:num_CV_folds
                    for i=1:size_chunk
                        for j=1:size_chunk
                            % In the adjacency matrix, cell (i,j) equals 1 if voxels (i,j) belong to the same cluster and 0 otherwise
                         % This allows comparison of clustering results even if the labels of the same clusters in each result are different
                         % (e.g. if the occipital cluster in solution 1 is labeled as cluster number 4, and in solution 2 it's labeled as cluster 7)
                            if IDX_folds(current_chunk1(i),c)==IDX_folds(current_chunk2(j),c)
                                current_clustering_adjmats(i,j,c)=1;
                            end
                        end
                    end
                end

                % Computing the difference between adjacency matrix for different folds, and adding this difference to the sum matrix
                for c1=1:num_CV_folds
                    for c2=1:num_CV_folds
                        sum_diff_adjmats_folds(c1,c2) = sum_diff_adjmats_folds(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1)~=current_clustering_adjmats(:,:,c2))));
                        sum_common_connections_adjmats(c1,c2) = sum_common_connections_adjmats(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1) & current_clustering_adjmats(:,:,c2))));
                        sum_all_connections_adjmats(c1,c2) = sum_all_connections_adjmats(c1,c2) + (sum(sum(current_clustering_adjmats(:,:,c1) + current_clustering_adjmats(:,:,c2))));
                    end
                end
            end
        end

         % Calculating the average difference between adjacency matrices (across all folds pairs)
            sum_diff_adjmats_folds_all(K) = mean(sum_diff_adjmats_folds(~eye(num_CV_folds)));   % sum of all differences
            Dice_coefficient = sum_common_connections_adjmats * 2 ./ sum_all_connections_adjmats;
            Dice_coefficient_folds_all(K) = mean(Dice_coefficient(~eye(num_CV_folds)));    % Dice's coef is 1 for perfect match, 0 for no commonalities
     end

% plotting the stability results for all K values, to identify peaks
figure;
plot(Dice_coefficient_folds_all,'Color',[0 0 1],'Linewidth',2);
hold on
plot(Dice_coefficient_folds_all,'*','Color',[1 0 1],'Linewidth',2)
hold off
xlim([K_range_l K_range_h])
xlabel('K-values')
ylabel('Dice Coefficients')
saveas(gcf,fullfile(output_folder,'K_grid_search_dice_corficient_sFC_GM.png'))
K = str2double(cell2mat(inputdlg('Choose the K-value','K-Value')));
close gcf
%% K-means clustering of subjects' mean connectivity matrix
% going over all numbers of clusters from K_range_l to K_range_h, and saving the results
temp=zeros(nx*ny*nz,1);
head=niftiinfo(fdp);
disp(['K = ' num2str(K)]);
id = kmeans(asFC_hc_tle, K,'distance','correlation','replicates',2*num_replicates);            % K-means clustering
save([output_folder,'cl_id_sFC_HC_lrTLE_gm'],'id');
for ic1=1:1:K
    roi_idx=find(id==ic1);
    vids=cell2mat(roi_ids(roi_idx));
    temp(vids,:)=ic1;
end
cl_map=reshape(temp,nx,ny,nz);
niftisave(cl_map,[output_folder,'cl_map_sFC_HC_lrTLE_gm.nii'],head);
elseif(Elbow==1)
    K=11;%[4,5,6,7,8,9,11,13,15,17,19,21]; %11(GM),11(WM)
    dm='cosine';%'sqeuclidean';%'distance','correlation', 'cosine'
    nr=10;

    Elb=zeros(length(K),2);
    % Elb_tle=zeros(length(K),2);
    for ki=1:1:length(K)
        k=K(ki);
        load('wm_atlas_dFC_cl_w30_seed.mat');
        seed=zeros(k,size(dFC,2));
        for j=1:1:k
            id11=find(id==j);
            seed(j,:)=mean(dFC(id11,:),1);
        end
        [id,~,sd]=kmeans(dFC,k,'Distance',dm,'start',seed); %Time point clustering
     % [id,~,sd]=kmeans(dFC,k,'Distance',dm,'Replicates',nr); %Time point clustering
        Elb(ki,:)=[k sum(sd.^2)];
        save(strcat('gm_atlas_dFC_cl_w30_seed_',num2str(k)),'id');
     % [id_tle,~,sd_tle]=kmeans(dFC_tle,k,'Distance',dm,'Replicates',nr); %Time point clustering
        % Elb_tle(ki,:)=[k sum(sd_tle.^2)];
    end
    figure(1)
    plot(Elb(:,1),Elb(:,2),'Color',[0 0 1],'Linewidth',2)
    hold on
    plot(Elb(:,1),Elb(:,2),'*','Color',[1 0 1],'Linewidth',2)
    hold off
    % plot(Elb_tle(:,1),Elb_tle(:,2),'Color',[0 1 0],'Linewidth',2)
    % hold on
    % plot(Elb_tle(:,1),Elb_tle(:,2),'*','Color',[1 0 0],'Linewidth',2)
    % legend('HC','HC','TLE','TLE');
    ylabel('Distortion')
    xlabel('k','fontsize',11)
    title('Elbow Method')

    temp=zeros(nx*ny*nz,1);
    head=niftiinfo(fdp);
    for ic1=1:1:K
        roi_idx=find(id==ic1);
        vids=cell2mat(roi_ids(roi_idx));
        temp(vids,:)=ic1;
    end

    cl_map=reshape(temp,nx,ny,nz);
    niftisave(cl_map,[output_folder,'cl_map_gm.nii'],head);
    % cl_map_tle=reshape(temp_tle,nx,ny,nz);
    % niftisave(cl_map_tle,[output_folder,'cl_map_tle.nii'],head);

end

% id = kmeans(asFC_hc_tle, K,'distance','correlation','replicates',num_replicates);            % K-means clustering
% % save([output_folder,'cl_id_sFC_HC_lrTLE_wm'],'id');
% for ic1=1:1:K
%     roi_idx=find(id==ic1);
%     vids=cell2mat(roi_ids(roi_idx));
%     temp(vids,:)=ic1;
% end
% cl_map=reshape(temp,nx,ny,nz);

if(FDR==1)
    G=[];
    disp(strcat('No of subjects = ',num2str(nos)));
    load([output_folder,'cl_id_sFC_HC_lrTLE_gm.mat']);
    esub=[]; % Subs to be excluded
    esub1=[];
    reject_roi1=[];
    hcs=0;
    ltles=0;
    rtles=0;
    sFC_hc=zeros(K,K,length(lhc));
    sFC_ltle=zeros(K,K,length(lltle));
    sFC_rtle=zeros(K,K,length(lrtle));
    for isub=1:1:nos
        sub=dir_hc_tle(isub).name;
        disp(strcat('Reading Subject....',sub));
        func_data_name = strcat(fdn1,sub,fdn2);
        fdp=fullfile(dp_hc_tle,sub,func_data_name);
        func_data = double(niftiread(fdp));
        not=size(func_data,4);
        tvmat=reshape(func_data,nx*ny*nz,not);
        sub_cl_tc=zeros(not,K);
       
        for ic=1:1:K
            if(find(isub==lhc))
                roi_idx=find(id==ic);
            elseif(find(isub==lltle))
                roi_idx=find(id==ic);
            elseif(find(isub==lrtle))
                roi_idx=find(id==ic);
            end
            vids=cell2mat(roi_ids(roi_idx));
            tv_mat_cl=tvmat(vids,:);
            tv_sum=sum(tv_mat_cl,2);
            zid=find(tv_sum==0);
            tv_mat_cl(zid,:)=[];
            nid=find(isnan(tv_sum));
            tv_mat_cl(nid,:)=[];
            if(isempty(tv_mat_cl))
                reject_roi1=[reject_roi1 iroi];
                % break;
            end
            tv_mean=mean(tv_mat_cl);
            sub_cl_tc(:,ic)=tv_mean;
        end
        if(ic==K)

            sFC=corr(sub_cl_tc);
            if(find(isub==lhc))
                grp=0;
                hcs=hcs+1;
                sFC_hc(:,:,hcs)=sFC;
            elseif(find(isub==lltle))
                grp=1;
                ltles=ltles+1;
                sFC_ltle(:,:,ltles)=sFC;
            elseif(find(isub==lrtle))
                grp=1;
                rtles=rtles+1;
                sFC_rtle(:,:,rtles)=sFC;
            end
            G=[G;grp];
            if(find(isnan(sFC(:))))
                esub1=[esub1,isub];
            end

        else
            esub=[esub,isub];
        end
    end
   
    sFC_tle=cat(3,sFC_ltle,sFC_rtle);
    J=(mean(sFC_hc,3)-mean(sFC_tle,3)).^2./((std(sFC_hc,[],3)).^2+(std(sFC_tle,[],3)).^2);
    J_lr=(mean(sFC_ltle,3)-mean(sFC_rtle,3)).^2./((std(sFC_ltle,[],3)).^2+(std(sFC_rtle,[],3)).^2);
    
    lbls="GM" +string(1:K);
    figure(1)
    h=heatmap(round(J,3),'MissingDataColor','w','ColorLimits',[0 .40]);
    h.XDisplayLabels = lbls;
    h.YDisplayLabels = lbls;
    h.FontSize = 14;
    % clim([0 40])
    colormap(jet);
    title('Fisher Discriminant Ratio(J):HC and TLE')
    ylabel("Clusters")
    xlabel("Clusters")

    figure(2)
    h1=heatmap(round(J_lr,3),'MissingDataColor','w','ColorLimits',[0 .40]);
    h1.XDisplayLabels = lbls;
    h1.YDisplayLabels = lbls;
    h1.FontSize = 14;
    % clim([0 40])
    colormap('jet'); 
    title('Fisher Discriminant Ratio(J): lTLE and rTLE')
    ylabel("Clusters")
    xlabel("Clusters")
end



