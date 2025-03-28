
%% process_abundances_biovolumes_NBSSspectra allows to compute abundance and differents biovolumes for all groups
% this extract the following measure on the database:
% all particles abundance and biovolume (in spherical equivalent)
%   abundance in ind/m3 
%   biovolume are in mm3/m3
%   spectra y unit is m3
% warning: to plot the spectra log transform the YY values

function [base_spectres SStot] = process_abundances_biovolumes_NBSSspectra_IFCB(subset,base_spectres,smin,smax,k,uu,zoo_groups,i,fracnb);

test=strcmp(fracnb,'tot');

SStot = [];

m=length(zoo_groups); %m = number of annotated plankton groups 

%cleaning if there is differents annotations 
for j=1:m
    f=find(zoo_groups{j,:}=='-');
    zoo_groups{j,1}(f)='_';
end

%% ---------- dimensions des vecteurs classes de taille ---------
x(1)=smin;
t=1;
while x(t,1)<smax
    t=t+1;
    x(t,1)=x(1)*k^(t-1);
end
x1=diff(x);
x=log(x);
X=(x(1:end-1,1)+x(2:end,1))./2; % choix de la classe milieu
X1=x1; % taille des classes
nb_int = length(X);
%% choice of the compute method

method_list = {'Plain_Area_BV_spectra' 'Ellipsoid_BV_spectra' 'Summed_plain_Area_BV_spectra'};

for meth = 1:3

    clear mj mn sss vol pred ar aarea esd R3 area_int  aareaexc fferet % clear variables that will be reused later

    %extracting data 'area' 'major' 'minor' 'perimferet' 'area_exc':

    mj = subset.major; % (mm)
    mn = subset.minor; % (mm)               
    aarea = subset.area; % (mm2)     
    summedbiovolume=subset.summedbiovolume; %(mm3)

    %calculating volumes:

    if meth == 3;
    %calculation of summed biovolume from IFCB
    sss=summedbiovolume; %calculation of biovolume (mm3/m3)
    
    elseif meth == 2;
    %calculation of area of best fitting ellipse
    ar=pi.*(mj./2).*(mn./2);
    sss=(4/3)*pi.*(mj./2).*(mn./2).*(mn./2); %calculation of biovolume (mm3/m3)

    else
    area_int=aarea./pi;
    esd=2*(sqrt(area_int));
    R3=(esd./2).^3;
    sss=(4/3)*pi.*R3;
    end

    vol=double(sss); %vol=biovolume
    sss=double(log(sss)); %log of biovolume
    
    pred=subset.object_annotation_hierarchy;
    
    for j=1:m  
        id(:,j)=strcmp(zoo_groups{j,1},pred); %id matrix for initially identified groups
        id=double(id);
    end

    ID = id; 
    
    %conversion factor to have all in cubic meter: 

    conver = subset.conver; %per cubic meter

    if test==1 %tot
        eval(['base_spectres(i).tot.conver=conver;'])
    else %d1 and d2
        eval(['base_spectres(i).d' num2str(fracnb) '.conver=conver;'])
    end

%% compute for each annotated plankton groups (m) 

    for j=1:m

    %Abundance

        if test==1 %tot
            
            Ab(j)=sum(ID(:,j)).*conver; %abundance per fraction related to volume
            eval(['base_spectres(i).tot.Ab(j,1)=Ab(j);']) %abundance per fraction related to volume (#/m3)
            Abtot(j)=sum(ID(:,j)); %total abundance over the entire sampled depth (number # of organisms in total not related to the volume filtered)
            eval(['base_spectres(i).tot.Abtot(j,1)=Abtot(j);']) %number # of organisms in total for the scanned fraction not related to the filtered volume            
            
        else %d1 and d2 

            Ab(j)=sum(ID(:,j)).*conver; %abundance per fraction related to volume
            eval(['base_spectres(i).d' num2str(fracnb) '.Ab(j,1)=Ab(j);']) %abundance per fraction related to volume (#/m3)
            Abtot(j)=sum(ID(:,j)); %total abundance over the entire sampled depth (number # of organisms in total not related to the volume filtered)
            eval(['base_spectres(i).d' num2str(fracnb) '.Abtot(j,1)=Abtot(j);']) %number # of organisms in total for the scanned fraction not related to the filtered volume
        end

        if meth == 1
        sizes(:,j)=ID(:,j).*esd; 
        end
        
    %Biovolume

        if meth == 3 %Choose the biovolume methods: Ellipsoid_BV_area method, plain area, riddle area

            if test==1 %tot
                
                vol(isnan(vol)) = 0;
                Bv(j)=sum(ID(:,j).*vol).*conver ; %biovolume per fraction related to volume
                eval(['base_spectres(i).tot.Bv(j,1)=Bv(j);']) %biovolume per fraction related to volume 
                Bvtot(j)=sum(ID(:,j).*vol); %biovolume total over the entire sampled depth (= # of organisms in total not related to the volume filtered)
                eval(['base_spectres(i).tot.Bvtot(j,1)=Bvtot(j);']) %biovolume of organisms in total for the scanned fraction not related to the filtered volume
           
            else %d1 and d2 

                vol(isnan(vol)) = 0;
                Bv(j)=sum(ID(:,j).*vol).*conver ; %biovolume per fraction related to volume 
                eval(['base_spectres(i).d' num2str(fracnb) '.Bv(j,1)=Bv(j);']) %biovolume per fraction related to volume 
                Bvtot(j)=sum(ID(:,j).*vol); %biovolume total over the entire sampled depth (= # of organisms in total not related to the volume filtered)
                eval(['base_spectres(i).d' num2str(fracnb) '.Bvtot(j,1)=Bvtot(j);']) %biovolume of organisms in total for the scanned fraction not related to the filtered volume

            end

        end
        
        aaa=find(ID(:,j)==1);
        f=find(zoo_groups{j,1}=='-');
        zoo_groups{j,1}(f)='_';

        %spectra per fraction:

        %Presence of organisms of this category for this fraction
        if isempty(aaa) == 0;
            
            %---------- Abundance ---------------------
            eval(['[X,Yab_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,1);']);% clear aaa
            eval(['Yab_zoo_groups_',num2str(j),'=Yab_zoo_groups_',num2str(j),'.*conver;']) %'Ygroup'= size distribution for a group for a fraction
            
            % ---------- Biovolume ---------------------
            eval(['[X,Ybv_zoo_groups_',num2str(j),',Z,X1]=f_ZOO_spectrum(sss(aaa,1),vol(aaa,1),smin,smax,k,2);']);% clear aaa
            eval(['Ybv_zoo_groups_',num2str(j),'=Ybv_zoo_groups_',num2str(j),'.*conver;']) %'Ygroup'= size distribution for a group for a fraction
            
        %No organisms in this category for this fraction
        elseif isempty(aaa) == 1;

            %---------- Abundance ---------------------
            eval(['Yab_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);

             % ---------- Biovolume ---------------------
            eval(['Ybv_zoo_groups_',num2str(j),'=zeros(size(X,1),1);']);
        end


        if test==1 %tot 
            %load the data for the fractions into the database

            %---------- Abundance ---------------------
            eval(['base_spectres(i).tot.Yab(:,j) = Yab_zoo_groups_',num2str(j),';']);

            % ---------- Biovolume ---------------------
            eval(['base_spectres(i).tot.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),';']);
            
        else %d1 and d2 
            %load the data for the fractions into the database

            %---------- Abundance ---------------------
            eval(['base_spectres(i).d' num2str(fracnb) '.Yab(:,j) = Yab_zoo_groups_',num2str(j),';']);

            % ---------- Biovolume ---------------------
            eval(['base_spectres(i).d' num2str(fracnb) '.Ybv_',char(method_list(meth)),'(:,j) = Ybv_zoo_groups_',num2str(j),';']);
        end
        
    
    end

end

if test==1 %tot
    eval(['base_spectres(i).tot.X = X;']);
    eval(['base_spectres(i).tot.X1 = X1;']);
    eval(['base_spectres(i).tot.Zoo_groups = zoo_groups;']);
else %d1 and d2 
    eval(['base_spectres(i).d' num2str(fracnb) '.X = X;']);
    eval(['base_spectres(i).d' num2str(fracnb) '.X1 = X1;']);
    eval(['base_spectres(i).d' num2str(fracnb) '.Zoo_groups = zoo_groups;']);
end