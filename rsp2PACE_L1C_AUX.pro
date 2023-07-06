
; NAME:
;   rsp2PACE_L1C_AUX
;
; PURPOSE:
;   convert RSP L2 data to PACE L1C aux format
;
; :Categories:
;   remote sensing
;
; :program outline
;   -read in data structure .csv files
;   -Loop over dates and files per date
;       -create file
;       -create dimensions and file structure
;       -fill in define variables and write attributes
;       -write global attributes
;       -enter data in file
; :Usage:
;   To run: rsp2PACE_L1C_AUX,settings_file
;
;   Settings_file.csv contains paths and dates to process 
;   (use "all" to run through all dates in campaign folder)
;   Settings_file.csv also contains output folder, etc.
;   Settings_file.csv also contains path to PACE file structure folder
;
;   The PACE file structure is defined by .cvs files:
;       -Dimensions.csv
;       -Folders.csv
;       -Variables.csv
;       -(A file "Global_attributes.csv" can also be read in if items need to be copied, but is not used currently.)
;
;   These structure files are based on example L1C files made available at:
;   https://oceancolor.gsfc.nasa.gov/data/pace/
;
;   These structure files can be altered, but header item names need to stay the same!
;   folder and dimension pointers in the variables.csv file point to those in the folders.csv and dimension.csv files (in order)
;   
;   
; :Author:
;    Bastiaan van Diedenhoven   
;    Senior Scientist
;    SRON Netherlands Institute for Space Research
;    Niels Bohrweg 4
;    Leiden
;    Tel: --
;    Email: b.van.diedenhoven@sron.nl 
;
; :History:
;     original created: 6 July 2023
;     based on rsp2PACE_L1C_HARP2_layout_bandsel.pro
;
; :Notes
;       Uses H5_PARSE to open RSP hdf5 files
;       not all variables in original aux file are available
;       Uses original PACE AUX file variables names, but filled in with similar RSP data 
;
pro rsp2PACE_L1C_AUX,settings_file,$
    switch_stop_after_one=switch_stop_after_one,$; stop after one file to debug/test
    switch_no_talk=switch_no_talk,$;to make it shut up...
    switch_print_pace_variable_list=switch_print_pace_variable_list; to see the list of variables in the pace data

;---------------------------------
;; switches:
;---------------------------------
IF KEYWORD_SET(switch_no_RSP_variables)THEN switch_add_RSP_variables=0 $
    ELSE switch_add_RSP_variables=1; switch to copy additional RSP parameters to the data 

;---------------------------------
;; settings
;---------------------------------
;reads .csv file with all dates 
;Assuming RSP data is stored in folders named "prefix_yyymmdd_end"
;Use dates=all to loop through all dates in the path_rsp directory

settings=read_csv(settings_file,N_TABLE_HEADER=1,TABLE_HEADER=settings_header)
dates=STRCOMPRESS(settings.field1,/REMOVE_ALL)
path_rsp=STRCOMPRESS(settings.field2,/REMOVE_ALL)
path1=STRCOMPRESS(settings.field3,/REMOVE_ALL)
path2=STRCOMPRESS(settings.field4,/REMOVE_ALL)
path_out=STRCOMPRESS(settings.field5,/REMOVE_ALL)
file_out1=STRCOMPRESS(settings.field6,/REMOVE_ALL)
file_out2=STRCOMPRESS(settings.field7,/REMOVE_ALL)
path_PACE_structure_file=STRCOMPRESS(settings.field8[0],/REMOVE_ALL)
bins_across_track_copy=FIX(settings.field9[0])

;---------------------------------
;start of code
;------------------------------------------------------------------------------------------------
;get PACE HARP-2 data structure
IF ~KEYWORD_SET(switch_no_talk) THEN print,'Reading PACE data structure'

PACE_HARP2_L1C_dims=read_csv(path_PACE_structure_file+'Dimensions.csv',N_TABLE_HEADER=1,TABLE_HEADER=dim_header)
;PACE_HARP2_L1C_folders=read_csv(path_PACE_structure_file+'Folders.csv',N_TABLE_HEADER=1,TABLE_HEADER=folders_header) assumes all data in main folder
PACE_HARP2_L1C_vars=read_csv(path_PACE_structure_file+'Variables.csv',N_TABLE_HEADER=1,TABLE_HEADER=var_header)
var_fields=STRSPLIT(var_header,',',/EXTRACT)
PACE_HARP2_L1C_attr=read_csv(path_PACE_structure_file+'Global_attributes.csv',N_TABLE_HEADER=1,TABLE_HEADER=attr_header)

pointer_dim1=where(var_fields eq 'VAR_DIM_POINTER1')
pointer_ndim=where(var_fields eq 'VAR_NDIMS')
;pointer_folder=where(var_fields eq 'FOLDER')
pointer_nattr=where(var_fields eq 'VAR_NATTRIBUTES')
pointer_attr_name1=where(var_fields eq 'VAR_ATTRIBUTE_NAME1')
pointer_attr_value1=where(var_fields eq 'VAR_ATTRIBUTE_VALUE1')


;get all dates if needed 
;---------------------------------
IF(dates[0] eq 'all')THEN BEGIN
    finddates=FILE_BASENAME(FILE_SEARCH(path_rsp[0]+path1[0]+'????????'+path2[0],COUNT=ndates))
    parts=STRSPLIT(finddates,'_',/EXTRACT)
    dates=parts[1]
    path_rsp=REPLICATE(path_rsp[0],ndates)
    path1=REPLICATE(path1[0],ndates)
    path2=REPLICATE(path2[0],ndates)
    file_out1=REPLICATE(file_out1[0],ndates)
    file_out2=REPLICATE(file_out2[0],ndates)

    IF ~KEYWORD_SET(switch_no_talk)THEN print,'processing all dates in directory path_rsp: '
    IF ~KEYWORD_SET(switch_no_talk)THEN print,dates
ENDIF

;Loop over dates
;---------------------------------
ndates=n_elements(dates)
FOR idate=0,ndates-1 DO BEGIN
    IF ~KEYWORD_SET(switch_no_talk)THEN print,'     date',idate+1,'/',ndates,' :',dates[idate]
    path_date=path_rsp[idate]+path1[idate]+dates[idate]+path2[idate]+'/'
    files=FILE_BASENAME(FILE_SEARCH(path_date+'RSP*.h5',COUNT=nfiles))
    
    IF(nfiles eq 0)THEN stop,'no files found at '+path_date

 ;Loop over files within date folder
   FOR ifile=0,nfiles-1 DO BEGIN
            IF ~KEYWORD_SET(switch_no_talk)THEN print,'         file',ifile+1,'/',nfiles
            file_RSP=files[ifile]
        
        parts=STRSPLIT(file_RSP,'_',/EXTRACT)
        time_rsp=parts[2]
        
        ;get RSP data
        IF ~KEYWORD_SET(switch_no_talk)THEN print,'Reading RSP data: ',path_date+file_RSP
        data_RSP=h5_PARSE(path_date+file_RSP,/READ_DATA)
        bands_FWHM=[25.52821700	, 19.16967800	, 18.88557700	, 19.91962200	, 20.79478900,	21.03750800	, 59.67634500	, 80.29510500 ,	126.24256000]
    
        ;file_name
        file_out=file_out1[idate]+time_rsp+file_out2[idate]
        
        IF ~KEYWORD_SET(switch_no_talk)THEN print,'Writing output file:',path_out+file_out
        ;open file output file
        id = NCDF_CREATE(path_out+file_out, /CLOBBER, /NETCDF4_FORMAT)
        
        ;define dimensions and write
        bins_along_track            = data_RSP.dim_scans._nelements
        bins_across_track           = bins_across_track_copy
        nlevels                     = 42 ; could be used later

        dimensions_rsp = {bins_along_track:bins_along_track,bins_across_track:bins_across_track,nlevels:nlevels}
        ndims=n_elements(PACE_HARP2_L1C_dims.(0))
        dim_id=LONARR(ndims)
        FOR idim=0,ndims-1 DO BEGIN
            dim_id[idim]=NCDF_DIMDEF(id, PACE_HARP2_L1C_dims.(0)[idim], dimensions_rsp.(idim))
        ENDFOR  
        
        
        ;write groups
        ; ngroups=n_elements(PACE_HARP2_L1C_folders.(0))
        ; group_id=LONARR(ngroups)
        ; FOR igroup=0,ngroups-1 DO BEGIN
        ;     group_id[igroup]=NCDF_GROUPDEF(id, PACE_HARP2_L1C_folders.(0)[igroup])
        ; ENDFOR
        group_id=id
        
        ;define variables and write attributes
        NVARIABLES=n_elements(PACE_HARP2_L1C_vars.(0))
        var_id=LONARR(NVARIABLES)
        scale_factor=MAKE_ARRAY(NVARIABLES,value=1,/FLOAT)
        Add_offset=MAKE_ARRAY(NVARIABLES,value=0,/FLOAT)
        IF KEYWORD_SET(switch_print_pace_variable_list)THEN print,'PACE variables:'
        FOR ivar=0,NVARIABLES-1 DO BEGIN
            IF KEYWORD_SET(switch_print_pace_variable_list)THEN print,ivar,' ',PACE_HARP2_L1C_vars.(0)[ivar],' ',PACE_HARP2_L1C_vars.(1)[ivar]
            ;for attributes such as fill_values use the same type as the data
            typeswitch=[1,0,0,0]; use float as default
            IF(PACE_HARP2_L1C_vars.(1)[ivar] eq 'DOUBLE')THEN typeswitch=[0,1,0,0]
            IF(PACE_HARP2_L1C_vars.(1)[ivar] eq 'INT')THEN typeswitch=[0,0,1,0]
            
            ndims=PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]
         
            dim_var=INTARR(ndims)  
            FOR idim=0,ndims-1 DO dim_var[idim]=dim_id[PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar]]
            
            group_id_var=group_id;[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]]
            var_id[ivar]=NCDF_VARDEF(group_id_var,PACE_HARP2_L1C_vars.(0)[ivar],dim_var,FLOAT=typeswitch[0],DOUBLE=typeswitch[1],SHORT=typeswitch[2])
            for iatt=0,PACE_HARP2_L1C_vars.(pointer_nattr)[ivar]-1 DO BEGIN
                value1=PACE_HARP2_L1C_vars.(pointer_attr_value1+iatt)[ivar];all strings
                name1=PACE_HARP2_L1C_vars.(pointer_attr_name1+iatt)[ivar]
                typeswitch_write=typeswitch
                IF(name1 eq 'long_name' or name1 eq 'units')THEN typeswitch_write=[0,0,0,1];string
                IF(name1 eq 'scale_factor')THEN BEGIN    
                    typeswitch_write=[1,0,0,0]
                    scale_factor[ivar]=value1
                ENDIF
                IF(name1 eq 'add_offset')THEN BEGIN    
                    typeswitch_write=[1,0,0,0]
                    add_offset[ivar]=value1
                ENDIF            
                value=MAKE_ARRAY(1,value=value1,FLOAT=typeswitch_write[0],DOUBLE=typeswitch_write[1],INTEGER=typeswitch_write[2],STRING=typeswitch_write[3]) 
                NCDF_ATTPUT,group_id_var,var_id[ivar],name1,value
             ENDFOR
        ENDFOR
        
        ;write global attributes
        ;These correspond to attributes in PACE L1C format
        ;Additional RSP-specific attributes are in the RSP_specific_variables folder
        
        CALDAT,SYSTIME(/JULIAN,/UTC),DM, DD, DY, TH, TM, TS
        datetime=STRING(DY,DM,DD,TH,TM,TS,FORMAT='(I4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,"Z")')
        NCDF_ATTPUT,id,/global,'date_created',datetime
        
        NCDF_ATTPUT,id,/global,'instrument',data_RSP.INSTRUMENT._data
        NCDF_ATTPUT,id,/global,'time_coverage_start',data_RSP.START_UTC._data
        NCDF_ATTPUT,id,/global,'time_coverage_end',data_RSP.START_UTC._data
        NCDF_ATTPUT,id,/global,'product_name',file_out

        ;write global attributes from the .csv file
        FOR iattr=0,N_elements(PACE_HARP2_L1C_attr.(0))-1 DO $
            NCDF_ATTPUT,id,/global,PACE_HARP2_L1C_attr.(0)[iattr],PACE_HARP2_L1C_attr.(1)[iattr]
 

        ; go through variables and write RSP data
                
        ;lat, lon
        pace_vars=['latitude','longitude']
        RSP_vars=['COLLOCATED_LATITUDE','COLLOCATED_LONGITUDE']
        nmap=n_elements(RSP_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0]
            rsp_var_map=where(TAG_NAMES(data_RSP.GEOMETRY) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ibin_accross_track=0,bins_across_track-1 DO dataput[ibin_accross_track,*]=(data_RSP.GEOMETRY.(rsp_var_map)._data[*,0]-Add_offset[ivar])/scale_factor[ivar]
            NCDF_VARPUT,group_id,var_id[ivar],dataput
        ENDFOR
        
        ;cth
        pace_vars=['cth_water_cloud','cth_ice_cloud']
        RSP_vars=['CLOUD_TOP_ALTITUDE','CLOUD_TOP_ALTITUDE']
        nmap=n_elements(RSP_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0]
            rsp_var_map=where(TAG_NAMES(data_RSP.DATA) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ibin_accross_track=0,bins_across_track-1 DO BEGIN
                cth_m_to_km=data_RSP.DATA.(rsp_var_map)._data[*,0]/1e3
                dataput[ibin_accross_track,*]=(cth_m_to_km-Add_offset[ivar])/scale_factor[ivar]
            ENDFOR
            NCDF_VARPUT,group_id,var_id[ivar],dataput
        ENDFOR

        pace_vars=['water_cloud_fraction']
        RSP_vars=['CLOUD_TOP_ALTITUDE']; as proxy for cloud mask
        nmap=n_elements(RSP_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0]
            rsp_var_map=where(TAG_NAMES(data_RSP.DATA) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ibin_accross_track=0,bins_across_track-1 DO BEGIN
                cloudmask=MAKE_ARRAY(bins_along_track,/FLOAT)
                icld=where(data_RSP.DATA.(rsp_var_map)._data[*,0] gt 0.,ncld)
                IF(ncld gt 0)THEN cloudmask[icld]=1.
                dataput[ibin_accross_track,*]=cloudmask
            ENDFOR
            NCDF_VARPUT,group_id,var_id[ivar],dataput
        ENDFOR
        ;ice cloud fraction for now is filled with 0.

        ;cer and cot
        pace_vars=['cer_21_water_cloud','cot_21_water_cloud']
        RSP_vars=['CLOUD_BOW_DROPLET_EFFECTIVE_RADIUS','CLOUD_BOW_OPTICAL_THICKNESS']
        nmap=n_elements(RSP_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0]
            rsp_var_map=where(TAG_NAMES(data_RSP.DATA) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ibin_accross_track=0,bins_across_track-1 DO dataput[ibin_accross_track,*]=(data_RSP.DATA.(rsp_var_map)._data[*,0]-Add_offset[ivar])/scale_factor[ivar]
            
            NCDF_VARPUT,group_id,var_id[ivar],dataput
        ENDFOR


        IF ~KEYWORD_SET(switch_no_talk)THEN print,'Closing output file'
        ;ncdf_control, id, /endef
        NCDF_CLOSE, id ;

        IF KEYWORD_SET(switch_stop_after_one)THEN stop
    ENDFOR;files
ENDFOR;dates

end
