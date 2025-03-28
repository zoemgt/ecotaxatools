% process_spectra_and_save.m allows to get abundances/biovolumes and NBSS spectra calculations

%% setting for calculus of the size spectra 

smin=0.000000000001; %set the lower limit of the biovolume spectra that will be calculated
smax=10000; %set the upper limit of the biovolume spectra that will be calculated  (if ESD 26mm)
k=2^(1/4); %set logarithmic base used to calculate bins of the size spectra according to platt & denman's theory (1978): scaling exponent set to 0.25
uu=1; %change the first size class for the regression (and put option = 1 for performing on the spectra from uu values to the end)

%% starting to process the base

h=waitbar(0,'Computing abundance, biovolume and NBSS spectra...')

for i=1:m

    [nfrac, no_use]=size(base(i).fracids);
    
    waitbar(i / m, h, sprintf('Computing abundance, biovolume and NBSS spectra: %d of %d', i, m));

    for fracnb=1:nfrac
        eval(['subset=base(i).d' num2str(fracnb) ';'])
        [base SStot] = process_abundances_biovolumes_NBSSspectra_IFCB(subset,base,smin,smax,k,uu,zoo_groups,i,fracnb);
    end
end
close(h)

%% add the abundance/biovolume and the 3 NBSS biovolumes into the Zooscan_base 

for i=1:m 

    [nfrac, no_use]=size(base(i).fracids);
        
        Ab=[]; %abundance per groups (ind.m-3)
        Bv=[]; %biovolume per groups (mm3.m-3)
        Ybv_Plain_Area_BV_spectra=[]; %NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
        Ybv_Ellipsoid_BV_spectra=[]; %NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
        Yab=[];

    for fracnb=1:nfrac

        if fracnb==1 %tot 

            eval(['Ab=base(i).d' num2str(fracnb) '.Ab;'])
            eval(['Bv=base(i).d' num2str(fracnb) '.Bv;'])
            eval(['Ybv_Plain_Area_BV_spectra=base(i).d' num2str(fracnb) '.Ybv_Plain_Area_BV_spectra;'])
            eval(['Ybv_Ellipsoid_BV_spectra=base(i).d' num2str(fracnb) '.Ybv_Ellipsoid_BV_spectra;'])
            eval(['Yab=base(i).d' num2str(fracnb) '.Yab;'])
            
        else %d1 and d2 

            eval(['Ab=nansum([Ab base(i).d' num2str(fracnb) '.Ab],2);'])
            eval(['Bv=nansum([Bv base(i).d' num2str(fracnb) '.Bv],2);'])
            eval(['Ybv_Plain_Area_BV_spectra=nansum(cat(3,Ybv_Plain_Area_BV_spectra, base(i).d' num2str(fracnb) '.Ybv_Plain_Area_BV_spectra),3);'])
            eval(['Ybv_Ellipsoid_BV_spectra=nansum(cat(3,Ybv_Ellipsoid_BV_spectra, base(i).d' num2str(fracnb) '.Ybv_Ellipsoid_BV_spectra),3);'])
            eval(['Yab=nansum(cat(3,Yab, base(i).d' num2str(fracnb) '.Yab),3);'])
            
        end
    end

%% if there is two size fraction (d1 and d2), assemblate them to have the total abundance, biovolume and NBSS biovolumes of the sample 

    base(i).tot.Ab=Ab; %abundance per fraction related to volume
    
    base(i).tot.Bv=Bv; %biovolume per fraction related to volume
    base(i).tot.Ybv_Plain_Area_BV_spectra = Ybv_Plain_Area_BV_spectra;
    base(i).tot.Ybv_Ellipsoid_BV_spectra = Ybv_Ellipsoid_BV_spectra;
    base(i).tot.Yab = Yab;

%% add the parameters allowing the calculation of the size spectra

    eval(['base(i).tot.X = base(i).d' num2str(fracnb) '.X;'])
    eval(['base(i).tot.X1 = base(i).d' num2str(fracnb) '.X1;'])
    eval(['base(i).tot.Zoo_groups = base(i).d' num2str(fracnb) '.Zoo_groups;'])
    eval(['X = base(i).d' num2str(fracnb) '.X;'])
    sizelist=exp(X);%in mm3
    ESD=2*(sizelist*3/(4*pi)).^(1/3); %in mm
    ESD=ESD*1000; %in µm
    base(i).tot.ESDvector = ESD; 
    
end
%% create regroupped column 
%considering the volumes of water filtered by the nets, planktonic organisms abundance (ind m-3) and biovolume (mm3 m-3) can be also calculated for several regroupement levels of grouping:  
%living or non-living, and a functional/trophic group annotation 

%the affiliation living or non-living, and the functionnal/trophic annotation of each plankton group can be found into the excel: trophic_affiliation_of_organisms.xlsx
table_groupage=readtable('trophic_affiliation_of_organisms.xlsx','ReadVariableNames',false);  

table_groupage=table2cell(table_groupage); %if needing to add new functional/trophic annotation, modify and add your desired within the excel file

f = msgbox({'-1 do not feed';'1 phototrophs';'1.5 mixotrophs';'2 grazers';'2.5 omnivorous';'3 predators';'3.5 unknown trophic group'},'trophic group legend')
%note a 0.5 place is available for bacteria living from dissolved

%regrouping taxa per functional/trophic groups
Zoo_groups=table_groupage(:,1);

%% finding if temporary groups are used
%a number of organisms could not be reliably identified due to a lack of identification criteria and were therefore grouped into temporary groups (t00x) following similar morphological criteria

istemporary=0;
k = strfind(zoo_groups,'temporary_');
test=sum(cell2mat(k)); %number of temporary groups
test2=cellfun(@isempty,k);
test2=test2==0;

if  test>0 %presence of temporary groups 
    
    answer = questdlg('Your files includes one or several temporary "t00X" categories. Do you have any "functional/trophic" mapping existing for those?', ...
        'Temporary categories mapping', ...
        'Yes please load them','No please create them','No please ignore them (not recommended)','Yes please load them');

    switch answer

        case 'No please ignore them (not recommended)'

        zoo_groups(test2)=[];
        [n,p]=size(zoo_groups);

            %updating to remove temporary groups
            for i=1:m
                base(i).tot.Zoo_groups(test2)=[];
                base(i).tot.Ab(test2)=[];            
                base(i).tot.Bv(test2)=[];               
                base(i).tot.Ybv_Plain_Area_BV_spectra (:,test2)=[];
                base(i).tot.Ybv_Ellipsoid_BV_spectra(:,test2)=[];
                base(i).tot.Yab(:,test2)=[];
            end

        case 'Yes please load them'
            
            %uptading to add temporary groups 
            [file,path] = uigetfile('*.xlsx')
            addontemp=readtable([path file],'ReadVariableNames',false);  
            addontemp=table2cell(addontemp);
            Zoo_groups=[Zoo_groups; addontemp(:,1)];
            istemporary=1;

        case 'No please create them'
            
            %associate temporary groups with affiliations:living or non-living, and a functional/trophic group annotation 
            groups=table_groupage(:,2:end);
            [n,p]=size(groups);
            group1=unique(groups(:,1));
            group2=unique(groups(:,2));
            group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));

            I=cellfun(@isempty,k);
            I=I==0;
            new_taxa=zoo_groups(I);
            [p]=length(new_taxa);
            newfunctional=[];

            %Default groups
            defaultGroup1 = 'living';
            defaultGroup2 = 'unidentified';
            defaultGroup3 = '3.5';
            group1 = [defaultGroup1; setdiff(group1, defaultGroup1, 'stable')];
            group2 = [defaultGroup2; setdiff(group2, defaultGroup2, 'stable')];
            group3 = [defaultGroup3; setdiff(group3, defaultGroup3, 'stable')];
            
            for i=1:p
                
                settings = settingsdlg('Description', ['A new temporary taxonomic group have been found ' char(new_taxa(i))],...
                    'title' , 'New taxa functional mapping',...
                    'Alive' , group1 ,...
                    'plankton group' , group2 , ...
                    'trophic group' , group3, ...
                    'WindowWidth' , 800)
                    
                newfunctional=[newfunctional settings];
            end


            newfunctional=struct2table(newfunctional)
            %newfunctional(:,4)=[];

            %updating the new associations into the xls reference list
            addontemp=[new_taxa table2cell(newfunctional)];
            tosave=array2table(addontemp);
            [file,path] = uiputfile('temporarymapping_instrument_net_location.xlsx');
            filename = fullfile(path,file);
            writetable(tosave,filename,'WriteVariableNames',0);
            Zoo_groups=[Zoo_groups; addontemp(:,1)];           
            istemporary=1;

    end

end

%% finding if "new" taxonomic groups are present withtout association of living the affiliation 

[n,p]=size(zoo_groups);
new_taxa={};p=0;

for i=1:n
    J=strcmp(char(zoo_groups(i)),Zoo_groups);
    if sum(J)==0
    p=p+1;
    new_taxa(p,1)=zoo_groups(i);
    end
end

clear Zoo_group

%associate new taxonomic groups with affiliations:living or non-living, and a functional/trophic group annotation
groups=table_groupage(:,2:end);
[n,p]=size(groups);

            group1=unique(groups(:,1));
            group2=unique(groups(:,2));
            group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));

[p]=length(new_taxa);
newfunctional=[];

for i=1:p
    
    settings = settingsdlg('Description', ['A new taxonomic group have been found' char(new_taxa(i))],...
        'title' , 'New taxa functional mapping',...
                    'Alive' , group1 ,...
                    'plankton group' , group2 , ...
                    'trophic group' , group3,  ...
                    'WindowWidth' , 800 )
    newfunctional=[newfunctional settings];
end

if p>0

newfunctional=struct2table(newfunctional);

%updating the new associations into the xls reference list
addon=[new_taxa table2cell(newfunctional)];
addon(:,5)=[]; 
table_groupage=[table_groupage; addon];
tosave=array2table(table_groupage);
cd(directoryoftoolbox);
writetable(tosave,'trophic_affiliation_of_organisms.xlsx','WriteVariableNames',0);
cd(folder);

end

if istemporary==1
   if size(addontemp, 2) >= 5
   addontemp(:,5)=[]; 
   end
   table_groupage=[table_groupage; addontemp];
end

%% function f_regroup_all.m allows to calculate abudance, biovolume and NBSS biovolume for each affilation 

for i=1:m
    [base_regroup] = f_regroup_all_IFCB(table_groupage,base(i).tot);
    base(i).regroupped=base_regroup;
end


%% preparing resume files on abundance / biovolume per taxa per sample

Ab_resume=[];
Bv_resume=[];
samplelist=[];

for i=1:m
    Ab_resume=[Ab_resume;base(i).regroupped.Ab];
    Bv_resume=[Bv_resume;base(i).regroupped.Bv];
    samplelist=[samplelist;base(i).sampleid];
    Zoo_groups=base(i).regroupped.Zoo_groups;
    base(i).regroupped.Idlistshort=Idlistshort;
end

%% computing NBSS associated values: NBSS trunked for computing slopes and intercept 

% Extract NBSS with ellipsoid biovolume

    for i = 1:numel(base)
        base(i).regroupped.NBSSliving = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, base(i).regroupped.livingplace);
        base(i).regroupped.NBSSpft = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, (base(i).regroupped.planktongroupplace(1):base(i).regroupped.planktongroupplace(2)));
        base(i).regroupped.NBSStrophic = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, (base(i).regroupped.trophicplace(1):base(i).regroupped.trophicplace(2)));
    end

% Extract NBSS with plain area biovolume

%     for i = 1:numel(base)
%         base(i).regroupped.NBSSliving = base(i).regroupped.Ybv_Plain_Area_BV_spectra(:, base(i).regroupped.livingplace);
%         base(i).regroupped.NBSSpft = base(i).regroupped.Ybv_Plain_Area_BV_spectra(:, (base(i).regroupped.planktongroupplace(1):base(i).regroupped.planktongroupplace(2)));
%         base(i).regroupped.NBSStrophic = base(i).regroupped.Ybv_Plain_Area_BV_spectra(:, (base(i).regroupped.trophicplace(1):base(i).regroupped.trophicplace(2)));
%     end

% NBSS to NBSS called 'trunked' 

%NBSS of total living plankton
for i = 1:numel(base)
        data2 = base(i).regroupped.NBSSliving;

        % Remove before max
        I = find(data2 == nanmax(data2));
        data2(1:I) = NaN;
        I = data2 == 0;
        data2(I) = NaN;

        % Remove solitary observations
        I2 = find(diff([0; isnan(data2)]) == 1);
        nanlength = [];
        for j = 2:(length(I2)-1)
            test2 = sum(isnan(data2(I2(j):I2(j+1)-1)));
            nanlength = [nanlength; test2];
        end
        nanlength = [0; nanlength; 0];
        I3 = nanlength > 5; % Threshold  at 5, can be change manually here 
        data2(I2(I3):end) = NaN;

        % Save NBSS trunked
        base(i).regroupped.NBSSliving_trunked = data2;
    end

%NBSSpft, NBSS for each plankton groups 
for i = 1:numel(base)
        data2 = base(i).regroupped.NBSSpft;
        
        % Remove before max
        for m = 1:size(data2,2);
        data3 = data2(:,m);
        I = find(data3 == nanmax(data3));
        data3(1:I) = NaN;
        I = data3 == 0;
        data3(I) = NaN;

        % Remove solitary observations
        I2 = find(diff([0; isnan(data3)]) == 1);
        nanlength = [];
        for j = 2:(length(I2)-1)
            test2 = sum(isnan(data3(I2(j):I2(j+1)-1)));
            nanlength = [nanlength; test2];
        end
        nanlength = [0; nanlength; 0];
        I3 = nanlength > 5; % Threshold  at 5, can be change manually here
        data3(I2(I3):end) = NaN;
        data2(:,m) = data3;
        end

        % Save NBSS trunked
        base(i).regroupped.NBSSpft_trunked = data2;
    end

%NBSStrophic, NBSS for each trophic groups
for i = 1:numel(base)
        data2 = base(i).regroupped.NBSStrophic;
        
        % Remove before max
        for m = 1:size(data2,2);
        data3 = data2(:,m);
        I = find(data3 == nanmax(data3));
        data3(1:I) = NaN;
        I = data3 == 0;
        data3(I) = NaN;

        % Remove solitary observations
        I2 = find(diff([0; isnan(data3)]) == 1);
        nanlength = [];
        for j = 2:(length(I2)-1)
            test2 = sum(isnan(data3(I2(j):I2(j+1)-1)));
            nanlength = [nanlength; test2];
        end
        nanlength = [0; nanlength; 0];
        I3 = nanlength > 5; % Threshold  at 5, can be change manually here
        data3(I2(I3):end) = NaN;
        data2(:,m) = data3;
        end

        % Save NBSS trunked        
        base(i).regroupped.NBSStrophic_trunked = data2;
    end

% Computing slopes and intercept for each NBSS 

for i = 1:numel(base)
    
    X = base(i).tot.X;
    NBSS = log10(base(i).regroupped.NBSSliving_trunked);   

    % Exclude zero or NaN values
    valid_indices = ~isnan(X) & ~isnan(NBSS);
    valid_X = X(valid_indices);
    valid_NBSS = NBSS(valid_indices);   

    % Check whether there are enough points to fit a linear regression
    if sum(valid_indices) >= 3
        coef = polyfit(valid_X, valid_NBSS, 1);
        base(i).regroupped.NBSSliving_slope = coef(:,1); 
        base(i).regroupped.NBSSliving_intercept = coef(:,2);
    elseif sum(valid_indices) < 3
        base(i).regroupped.NBSSliving_slope = NaN;
        base(i).regroupped.NBSSliving_intercept = NaN;
    end

end

%% saving the final bases and resume files

instrument='to complete'; 

prompt = {'instrument:','project:'};
title = 'Save base under the name:';
dims = [1 35];
definput = {instrument,'to complete'};
answer = inputdlg(prompt,title,dims,definput)

save(['base_' char(answer(1)) '_' char(answer(2))],'base')