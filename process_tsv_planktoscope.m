%process_tsv_zooscan_new.m allows to assemble the data from the tsv files into a database called "base".  

function [base,m,zoo_groups,Idlistshort]=process_tsv_planktoscope

%% create the bases

A=dir('*.tsv');
filenames={A.name};
[no_use,m]=size(filenames);

base=[]; %create the base which gathers the tsv data
Idlist=[]; %create Idlist; Idlist corresponds to the plankton groups list with taxonomic annotation hierarchy included 
Idlistshort=[]; %create Idlistshort; Idlistshort corresponds to the plankton groups list without taxonomic annotation hierarchy 

h=waitbar(0,'Load tsv files and create base...');

for i=1:m

S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
sample=unique(S.sample_id); 
sample_program=unique(S.sample_project);

waitbar(i / m, h, sprintf('Load tsv files and create base: %d of %d', i, m));

%% homogenize the writing of the taxonomic annotation hierarchy
    
[n,no_use]=size(S);
    
    for j=1:n;
        f=find(S.object_annotation_hierarchy{j,:}=='-');
        S.object_annotation_hierarchy{j,1}(f)='_';
        
        f=find(S.object_annotation_hierarchy{j,:}=='>');
        S.object_annotation_hierarchy{j,1}(f)='_';
    end

%% construction of the base

    base(i).sampleid=sample; %name of the sample

    base(i).ship=unique(S.sample_ship); %name of the ship 
    base(i).scientificprog=unique(S.sample_project); %name of the scientific program or project

    base(i).date=unique(S.object_date); %sampling date
    base(i).time=unique(S.object_time); %sampling time 
    base(i).latitude=double(unique(S.object_lat)); %sampling latitude 
    base(i).longitude=double(unique(S.object_lon)); %sampling longitude
    if ismember('object_lat_end', S.Properties.VariableNames) %sampling latitude or latitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).latitude_end = unique(S.object_lat_end);
    end 
    if ismember('object_lon_end', S.Properties.VariableNames) %sampling longitude or longitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).longitude_end=unique(S.object_lon_end); 
    end

    base(i).nettype=unique(S.sample_sampling_gear); %type/name of the net

    base(i).zmax=(unique(S.object_depth_max)); %maximum depth of the net when deployed
    base(i).zmin=(unique(S.object_depth_min)); %minimum depth of the net when deployed 

%% getting pixel size in micrometer and converting in mm

if isnumeric(unique(S.process_pixel)) %dimension of the side of a pixel in the scanned image [mm]
pixelsize=(unique(S.process_pixel))/1000;  %in mm/pixel
else
pixelsize=str2num((cell2mat(unique(S.process_pixel))))/1000;  %in mm/pixel
end

base(i).pixelsize=pixelsize;

%% construction of the base - volume section 
      
    % volconc in m3 
    base(i).volconc=unique(S.sample_concentrated_sample_volume)/1000000; %concentrated or diluted water volume [mL], convert in m3
    
    % dilution_factor
    if isnan(S.sample_dilution_factor)
    base(i).dilution_factor = 1;
    else
    base(i).dilution_factor = 1/unique(S.sample_dilution_factor); % concentrating factor
    end
    
    % nb_frame
    base(i).nb_frame=unique(S.acq_nb_frame);
    
    % celltype
    base(i).celltype=unique(S.acq_celltype)/1000000; %[micrometers], convert in m 
    
    % resolution
    base(i).resolution=unique(S.acq_camera_resolution); %[pixel]
    
    % sample_initial_col_vol_m3
    base(i).sample_initial_col_vol_m3=unique(S.sample_total_volume)./1000; %[l], convert in m3  
    
    %In the PlanktoScope, the volume imaged is calculated directly by the instrument 
    base(i).volume_imaged_processed=unique(S.acq_imaged_volume)/1000000 ; %[mL], convert in m3

%% FracId (d1, d2, ... or tot) 

base(i).fracids=unique(S.acq_id); %see if there is size fraction
[nfrac, no_use]=size(base(i).fracids);

    for fracnb=1:nfrac 
        %nfrac=1, no size fraction (tot)
        %nfrac=2, size fraction was used 

        I=strcmp(S.acq_id,base(i).fracids(fracnb));

        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));']); %taxonomical annotation hierarchy
        
        %operation to create Idlistshort: 
        [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]); 
        Idlistshort=[Idlistshort; S.object_annotation_category(I)];
        Idlistshort=Idlistshort(ia);
        
        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);']);
        I2=strfind(S.object_annotation_status(I),'validated');
        index = (cellfun(@isempty,  I2) ==0);
        eval(['base(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);']); %percentage of validated vignettes in the sample

        if iscell(S.object_major)==0

            eval(['base(i).d' num2str(fracnb) '.major=S.object_major(I)*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
            eval(['base(i).d' num2str(fracnb) '.minor=S.object_minor(I)*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism)
            eval(['base(i).d' num2str(fracnb) '.area_exc=S.object_area_exc(I)*(pixelsize^2);']); %ratio between the distance between the foci of the ellipse and the length of its major axis (elongation of the organism) 
            eval(['base(i).d' num2str(fracnb) '.area=S.object_area(I)*(pixelsize^2);']); %number of pixels in the region (size of the organism)
            eval(['base(i).d' num2str(fracnb) '.ESD=2*(((S.object_area(I)*(pixelsize^2))/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)

        else

            eval(['base_Zooscan(i).d' num2str(fracnb) '.major=cellfun(@str2num,S.object_major(I))*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
            eval(['base_Zooscan(i).d' num2str(fracnb) '.minor=cellfun(@str2num,S.object_minor(I))*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism)
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area_exc=cellfun(@str2num,S.object_area_exc(I))*(pixelsize^2);']); %ratio between the distance between the foci of the ellipse and the length of its major axis (elongation of the organism)
            eval(['base_Zooscan(i).d' num2str(fracnb) '.area=cellfun(@str2num,S.object_area(I))*(pixelsize^2);']); %number of pixels in the region (size of the organism)
            eval(['base_Zooscan(i).d' num2str(fracnb) '.ESD=2*((cellfun(@str2num,S.object_area(I))*(pixelsize^2)/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)          

        end

        eval(['base(i).d' num2str(fracnb) '.conver=  base(i).dilution_factor.*base(i).volconc./(base(i).volume_imaged_processed.*base(i).sample_initial_col_vol_m3);']); %[m3] 

    end

end
      
close(h);

[n,m]=size(base);
zoo_groups=Idlist;  

end
