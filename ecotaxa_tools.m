
close all
clear all

[directoryoftoolbox,name,ext] = fileparts(mfilename('fullpath')); 

%% load the tsv files

%select the folder containing all the .tsv files exported from ecotaxa 
%"split in multiples files by sample": one sample equal to one .tsv files

f = msgbox('Please select the folder containing the *.tsv files');
uiwait(f)
folder=uigetdir;
cd(folder)

A=dir('*.tsv');

    %if *tsv.files is uncorrect, warning message:
    a=length(A);
    if a==0
    msgbox('The folder does not contain any *.tsv files from ecotaxa, please extract your data before from the Ecotaxa software with the option: "separate by sample"')
    return
    end

%% select the instrument used

list = {'Zooscan','Flowcam','IFCB','PlanktoScope','UVP'};
    
[indx2,tf] = listdlg('PromptString','What instrument was used?','ListString',list,'SelectionMode','single','ListSize',[300,150]);

%% load the tsv files, create a database of abundance, biovolumes and NBSS spectra by taxonomic groups, plankton groups, and trophic groups
warning('off', 'all');

    if indx2==1 %for Zooscan
        [base,m,zoo_groups,Idlistshort]=process_tsv_zooscan; 
        process_spectra_and_save;  

    elseif indx2==2 %for Flowcam 
        [base,m,zoo_groups,Idlistshort]=process_tsv_flowcam;
        process_spectra_and_save;   

    elseif indx2==3 %for IFCB 
        [base,m,zoo_groups,Idlistshort]=process_tsv_ifcb;
        process_spectra_and_save_IFCB;

     elseif indx2==4 %for PlanktoScope 
        [base,m,zoo_groups,Idlistshort]=process_tsv_planktoscope;
        process_spectra_and_save;  

     elseif indx2==5 %for UVP    
        process_tsv_uvp; %step: add Ecotaxa export (size and other features)
        %process_tsv_uvp_ecopart; %step1: Ecopart detailed export (no size)

    end

