% f_regroup_all_new allows to group zooplankton function to calculate abundance and spectra
 
%% table_grouping: 

% an xls table with 3 columns
% 1st column = original groups
% 2nd column = new group names living vs non_living, 
% 3rd column = new PFT group names
% 2 abundances : abundance table (ind/m3)
% # rows = # of zoo groups identified at base 
% # columns = # of samples
% 3 biovolumes: biovolume table (mm3/m3)
% # rows = # zoo groups identified at baseline 
% # columns = # of samples
% 4 size_spectra: SStot in output of the creation routine of the base_spectra 

%% function 

function [base_regroup] = f_regroup_all_new(table_groupage,base)

tab=table_groupage;
tab_orig=base.Zoo_groups;
id=tab_orig;

g_planktongroup=tab(:,3);

to_keep=[];

for i=1:size(id,1)
    a=strmatch(id(i),tab(:,1),'exact');
    if isempty(a)==0
        if length(a)>1
        disp(['warning' id(i) 'is in double, please check reference excell file before going further']);
        end
        to_keep=[to_keep;a];
    end
end

tab=tab(to_keep,:);

g_orig=tab(:,1);
g_living=tab(:,2);
g_planktongroup=tab(:,3);
g_trophic=tab(:,4);

Ab=base.Ab;
Ybv_Plain_Area_BV_spectra=base.Ybv_Plain_Area_BV_spectra;
Ybv_Riddled_Area_BV_spectra=base.Ybv_Riddled_Area_BV_spectra;
Bv=base.Bv;
ss=base.Ybv_Ellipsoid_BV_spectra;
Yab=base.Yab;

clear SA
%% total grouping = in 'all' into base 

Ab_t=sum(Ab);
Ybv_Plain_Area_BV_spectra_t=sum(Ybv_Plain_Area_BV_spectra,2);
Ybv_Riddled_Area_BV_spectra_t=sum(Ybv_Riddled_Area_BV_spectra,2);
Bv_t=sum(Bv);
sst=sum(ss,2);
Yab_t=sum(Yab,2);

%% grouping living and no living together 

new_groups=unique(g_living);

for i=1:size(new_groups,1)
    
    ng=new_groups(i);
    f_ng=strmatch(ng,g_living);

    ss_r=ss(:,f_ng);
    Ab_r=Ab(f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    Yab_r=Yab(:,f_ng);
    
    Ab_n1(i)=sum(Ab_r);
    Ybv_Plain_Area_BV_spectra_n1(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n1(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n1(i)=sum(Bv_r);
    ssn1(:,i)=sum( ss_r,2);
    Yab_n1(:,i)=sum(Yab_r,2);
    
    clear SSJ ssj ss_r f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r  Ab_r Yab_r

end

%% grouping plankton groups 

new_groups=unique(g_planktongroup);

for i=1:size(new_groups,1)
    ng=new_groups(i);
    f_ng=strmatch(ng,g_planktongroup);
   
    ss_r=ss(:,f_ng);
    Ab_r=Ab(f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    Yab_r=Yab(:,f_ng);
    
    Ab_n2(i)=sum(Ab_r);
    Ybv_Plain_Area_BV_spectra_n2(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n2(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n2(i)=sum(Bv_r);
    ssn2(:,i)=sum( ss_r,2) ;
    Yab_n2(:,i)=sum(Yab_r,2);
 
    clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r Yab_r

end
    
 %% grouping trophic groups 

new_groups=unique(cell2mat(g_trophic));

for i=1:size(new_groups,1)

    ng=new_groups(i);
    f_ng=ng==cell2mat(g_trophic);
   
    ss_r=ss(:,f_ng);
    Ab_r=Ab(f_ng);
    Ybv_Plain_Area_BV_spectra_r=Ybv_Plain_Area_BV_spectra(:,f_ng);
    Ybv_Riddled_Area_BV_spectra_r=Ybv_Riddled_Area_BV_spectra(:,f_ng);
    Bv_r=Bv(f_ng);
    Yab_r=Yab(:,f_ng);
    
    Ab_n3(i)=sum(Ab_r);
    Ybv_Plain_Area_BV_spectra_n3(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
    Ybv_Riddled_Area_BV_spectra_n3(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
    Bv_n3(i)=sum(Bv_r);
    ssn3(:,i)=sum( ss_r,2);
    Yab_n3(:,i)=sum(Yab_r,2);

    clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r

end

%% regroup 

% create base with abundance, biovolume and spectra: 
Ab_regroup=[Ab' Ab_t Ab_n1 Ab_n2 Ab_n3];
Ybv_Plain_Area_BV_spectra_regroup=[Ybv_Plain_Area_BV_spectra Ybv_Plain_Area_BV_spectra_t Ybv_Plain_Area_BV_spectra_n1 Ybv_Plain_Area_BV_spectra_n2 Ybv_Plain_Area_BV_spectra_n3];
Ybv_Riddled_Area_BV_spectra_regroup=[Ybv_Riddled_Area_BV_spectra Ybv_Riddled_Area_BV_spectra_t Ybv_Riddled_Area_BV_spectra_n1 Ybv_Riddled_Area_BV_spectra_n2 Ybv_Riddled_Area_BV_spectra_n3];
Bv_regroup=[Bv' Bv_t Bv_n1 Bv_n2 Bv_n3];
Ybv_Ellipsoid_BV_spectra_regroup=[ss sst ssn1 ssn2 ssn3];
Yab_regroup=[Yab Yab_t Yab_n1 Yab_n2 Yab_n3];

% add into the base_regroup:
base_regroup=[];
base_regroup.Ab=Ab_regroup;
base_regroup.Ybv_Plain_Area_BV_spectra=Ybv_Plain_Area_BV_spectra_regroup;
base_regroup.Ybv_Riddled_Area_BV_spectra=Ybv_Riddled_Area_BV_spectra_regroup;
base_regroup.Bv=Bv_regroup;
base_regroup.Ybv_Ellipsoid_BV_spectra=Ybv_Ellipsoid_BV_spectra_regroup;
base_regroup.Yab=Yab_regroup;

base_regroup.Zoo_groups=[g_orig;{'all'};unique(g_living);unique(g_planktongroup);cellstr(num2str(unique(cell2mat(g_trophic))))];
base_regroup.originalplace=[1 length(g_orig)];
base_regroup.allplace=length(g_orig)+1;
base_regroup.livingplace=length(g_orig)+2;
base_regroup.notlivingplace=length(g_orig)+3;
base_regroup.planktongroupplace=[length(g_orig)+4 length(g_orig)+3+length(unique(g_planktongroup))];
base_regroup.trophicplace=[length(g_orig)+3+length(unique(g_planktongroup))+1 length(g_orig)+3+length(unique(g_planktongroup))+length(unique(cell2mat(g_trophic)))];

