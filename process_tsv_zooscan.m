%process_tsv_zooscan.m allows to assemble the data from the tsv files into a database called "base".  

function [base,m,zoo_groups,Idlistshort]=process_tsv_zooscan

%% create the base

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
sample_program=unique(S.sample_program);

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
    base(i).scanop=unique(S.sample_scan_operator); %name of the scan operator into the ZooScan 

    base(i).ship=unique(S.sample_ship); %name of the ship 
    base(i).scientificprog=unique(S.sample_program); %name of the scientific program 

    base(i).stationid=(unique(S.sample_stationid)); %name of the station 
    base(i).date=unique(S.object_date); %sampling date
    base(i).time=unique(S.object_time); %sampling time 
    base(i).latitude=double(unique(S.object_lat)); %sampling latitude 
    base(i).longitude=double(unique(S.object_lon)); %sampling longitude
    base(i).depth=(unique(S.sample_bottomdepth)); %depth of the sea at the station 
   
    base(i).CTDref=unique(S.sample_ctdrosettefilename); %if a CTD has been deployed in the same sample station, reference of the associated CTD deployment
    base(i).otherref=unique(S.sample_other_ref); %others references associated to this station or sample

    base(i).barcode=unique(S.sample_barcode); %barcode sample
    
    base(i).nettype=unique(S.sample_net_type); %type/name of the net
    base(i).townb=(unique(S.sample_tow_nb)); %number of tow (number of net deployments)
    base(i).towtype=(unique(S.sample_tow_type)); %type of tow (i.e horizontal, vertical) 
    base(i).netmesh=(unique(S.sample_net_mesh)); %size of the net mesh
    base(i).netsurf=(unique(S.sample_net_surf)); %surface of the net
    base(i).zmax=(unique(S.sample_zmax)); %maximum depth of the net when deployed
    base(i).zmin=(unique(S.sample_zmin)); %minimum depth of the net when deployed 

    base(i).vol=(unique(S.sample_tot_vol)); %volume of seawater filtered by the net during its deployment (m3) 

    base(i).samplecomments=unique(S.sample_comment); %comments 

%% getting pixel size in micrometer and converting in mm
    
pixelsize=(unique(S.process_particle_pixel_size_mm)); %in mm/pixel

%% FracId (d1 and d2 or tot) 
%For a mesh size of 200 µm or less, to avoid underestimation of large and rare organisms, the sample is sieved through a 1000 µm and a 100 µm mesh to obtain 2 size fractions: "d1" (>1000 µm) and "d2" (100-1000 µm). 
%For a mesh size > 200 µm, the d1 and d2 fractions are not necessary. The whole sample is scanned directly: "tot". 

base(i).fracids=unique(S.acq_id); %see if there is size fraction
[nfrac, no_use]=size(base(i).fracids);

    for fracnb=1:nfrac 
        %nfrac=1, no size fraction (tot)
        %nfrac=2, size fraction was used (d1 >1000 µm and d2 100-1000 µm) 

        I=strcmp(S.acq_id,base(i).fracids(fracnb));
        eval(['base(i).d' num2str(fracnb) '.Fracmin=unique(S.acq_min_mesh(I));']); %minimum fraction size (µm)
        eval(['base(i).d' num2str(fracnb) '.Fracsup=unique(S.acq_max_mesh(I));']); %maximum fraction size (µm) 
        eval(['base(i).d' num2str(fracnb) '.Fracnb= unique(S.acq_sub_part(I));']); %sample fraction made by the motoda and analysed into the ZooScan
        eval(['base(i).d' num2str(fracnb) '.Resolution=unique(S.process_img_resolution(I));']); %resolution of the ZooScan 
        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));']); %taxonomical annotation hierarchy
        
        %operation to create Idlistshort: 
        [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]); 
        Idlistshort=[Idlistshort; S.object_annotation_category(I)];
        Idlistshort=Idlistshort(ia);
        
        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);']);
        I2=strfind(S.object_annotation_status(I),'validated');
        index = (cellfun(@isempty,  I2) ==0);
        eval(['base(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);']); %percentage of validated vignettes in the sample

        eval(['base(i).d' num2str(fracnb) '.major=S.object_major(I)*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
        eval(['base(i).d' num2str(fracnb) '.minor=S.object_minor(I)*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism)
        eval(['base(i).d' num2str(fracnb) '.area_exc=S.object_area_exc(I)*(pixelsize^2);']); %ratio between the distance between the foci of the ellipse and the length of its major axis (elongation of the organism) 
        eval(['base(i).d' num2str(fracnb) '.area=S.object_area(I)*(pixelsize^2);']); %number of pixels in the region (size of the organism)
        eval(['base(i).d' num2str(fracnb) '.perimferet=S.object_feret(I)*pixelsize;']); %distance around the border of the region (size of the organism)

        eval(['base(i).d' num2str(fracnb) '.ESD=2*(((S.object_area(I)*(pixelsize^2))/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)
        
        eval(['base(i).d' num2str(fracnb) '.conver=(base(i).d' num2str(fracnb) '.Fracnb)./base(i).vol;']);

    end

end
      
close(h);

[n,m]=size(base);
zoo_groups=Idlist;  

end