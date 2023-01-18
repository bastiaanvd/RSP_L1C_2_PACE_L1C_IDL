
; NAME:
;   rsp2PACE_L1C
;
; PURPOSE:
;   convert RSP L1C data to PACE L1C format
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
;   To run: rsp2PACE_L1C,settings_file
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
;     original created: 1 September 2020
;     12 October 2020: changed to work on ER-2 data (time stamp 20190629T000839Z) instead of P3 data
;     10 December 2020: corrected solar constant definition and units
;     19 February 2021: Recoded for more flexibility; should now work on any current RSP L1C data 
;     10 March 2021: added setting.csv file 
;     17 September 2021: updated contact info
;     17 September 2021: added scattering angle
;     21 December 2022: corrected application of offset and scale factors
;     21 December 2022: added bins_across_track_copy to copy RSP data to emulate swath 
;     11 January 2022: Corrected fill values
;     11 January 2022: Made this version that writes out data using _bands_per_view=1, as HARP2 would have
;
; :Notes
;       Uses H5_PARSE to open RSP hdf5 files
;
;       The variables below are not included because they are missing from RSP data:
;           Viewing angles at instrument
;           altitude_variability
;           I_noise
;           Q_noise
;           U_noise
;           DOLP_noise

;
pro rsp2PACE_L1C_HARP2_layout_bandsel,settings_file,$
    switch_stop_after_one=switch_stop_after_one,$; stop after one file to debug/test
    switch_no_talk=switch_no_talk,$;to make it shut up...
    switch_print_pace_variable_list=switch_print_pace_variable_list, $; to see the list of variables in the pace data
    switch_no_RSP_variables=switch_no_RSP_variables;by default RSP variables are added to the data

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
PACE_HARP2_L1C_folders=read_csv(path_PACE_structure_file+'Folders.csv',N_TABLE_HEADER=1,TABLE_HEADER=folders_header)
PACE_HARP2_L1C_vars=read_csv(path_PACE_structure_file+'Variables.csv',N_TABLE_HEADER=1,TABLE_HEADER=var_header)
var_fields=STRSPLIT(var_header,',',/EXTRACT)
PACE_HARP2_L1C_attr=read_csv(path_PACE_structure_file+'Global_attributes.csv',N_TABLE_HEADER=1,TABLE_HEADER=attr_header)

;read in band selection determine ibandselection
PACE_HARP2_L1C_bands=read_csv(path_PACE_structure_file+'Band_selection.csv',N_TABLE_HEADER=1,TABLE_HEADER=bands_header)
ibandselection=where(PACE_HARP2_L1C_bands.FIELD2 eq 1,number_of_rsp_bands)
ang_reduce_fact=PACE_HARP2_L1C_bands.FIELD3

pointer_dim1=where(var_fields eq 'VAR_DIM_POINTER1')
pointer_ndim=where(var_fields eq 'VAR_NDIMS')
pointer_folder=where(var_fields eq 'FOLDER')
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
        ;number_of_rsp_bands = determined above using band_selection
        istart=data_RSP.Data.Unvignetted_Sector_begin._data
        iend=data_RSP.Data.Unvignetted_Sector_end._data
        number_of_rsp_views = iend-istart+1
        ;make a pointer to select which angles to select per band
        views_select=intarr(number_of_rsp_bands,number_of_rsp_views)
        FOR iwl=0,number_of_rsp_bands-1 DO BEGIN
            fact_reduce=ang_reduce_fact[ibandselection[iwl]]
            nreduced=number_of_rsp_views/fact_reduce
            nleftover=number_of_rsp_views-nreduced*fact_reduce
            FOR ifill=0,number_of_rsp_views-1,fact_reduce DO views_select[iwl,ifill+nleftover/2]=1
        ENDFOR
        find_number_of_views=where(views_select eq 1,number_of_views) 

        intensity_bands_per_view    = 1
        polarization_bands_per_view = 1
        bins_along_track            = data_RSP.dim_scans._nelements
        bins_across_track           = bins_across_track_copy


        dimensions_rsp = {number_of_views:number_of_views,intensity_bands_per_view:intensity_bands_per_view,polarization_bands_per_view:polarization_bands_per_view,bins_along_track:bins_along_track,bins_across_track:bins_across_track}
        ndims=n_elements(PACE_HARP2_L1C_dims.(0))
        dim_id=LONARR(ndims)
        FOR idim=0,ndims-1 DO BEGIN
            dim_id[idim]=NCDF_DIMDEF(id, PACE_HARP2_L1C_dims.(0)[idim], dimensions_rsp.(idim))
        ENDFOR  
        
        
        ;write groups
        ngroups=n_elements(PACE_HARP2_L1C_folders.(0))
        group_id=LONARR(ngroups)
        FOR igroup=0,ngroups-1 DO BEGIN
            group_id[igroup]=NCDF_GROUPDEF(id, PACE_HARP2_L1C_folders.(0)[igroup])
        ENDFOR
        
        
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
            
            group_id_var=group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]]
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
        NCDF_ATTPUT,id,/global,'sun_earth_distance',data_RSP.Platform.Solar_distance._data[0]
        NCDF_ATTPUT,id,/global,'time_coverage_start',data_RSP.START_UTC._data
        NCDF_ATTPUT,id,/global,'time_coverage_end',data_RSP.START_UTC._data
        NCDF_ATTPUT,id,/global,'product_name',file_out
        NCDF_ATTPUT,id,/global,'nadir_bin',MEDIAN(data_RSP.GEOMETRY.NADIR_INDEX._data[*,0]),/SHORT

        ;write global attributes from the .csv file
        FOR iattr=0,N_elements(PACE_HARP2_L1C_attr.(0))-1 DO $
            NCDF_ATTPUT,id,/global,PACE_HARP2_L1C_attr.(0)[iattr],PACE_HARP2_L1C_attr.(1)[iattr]
 

        ; go through variables and write RSP data
                
        ;view_time nadir
        pace_vars='nadir_view_time'
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars)
        ivar=ivar[0]
        
        dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
        dataput=MAKE_ARRAY(dim_ivar,/DOUBLE)
        FOR ipix=0,bins_along_track-1 DO BEGIN
            Nadir_index=data_RSP.Geometry.Nadir_index._data[ipix,0]
            MJD=data_RSP.Geometry.Measurement_Time._data[Nadir_index, ipix, 0]
            sec=(MJD-FLOOR(MJD))*86400.
            dataput[ipix]=sec
        ENDFOR
        NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        
        ;view_time_offsets
        pace_vars='view_time_offsets'
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars)
        ivar=ivar[0]
        
        dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
        FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
        dataput=MAKE_ARRAY(dim_ivar,/DOUBLE)
        FOR ipix=0,bins_along_track-1 DO BEGIN
            Nadir_index=data_RSP.Geometry.Nadir_index._data[ipix,0]
            diff=(data_RSP.Geometry.Measurement_Time._data[istart:iend, ipix, 0]-data_RSP.Geometry.Measurement_Time._data[Nadir_index, ipix, 0])*86400.
            icheckfill=where(ABS(diff) gt 1e4,nfill)
            IF(nfill gt 0)THEN diff[icheckfill]=-999.

            dataput_harp=dblarr(number_of_views)
            i1=0
            for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                i2=i1+number_of_views_wl-1
                dataput_harp[i1:i2]=diff[views_select_wl]
                i1=i2+1
            endfor
            FOR ibin_accross_track=0,bins_across_track-1 DO dataput[*,ibin_accross_track,ipix]=dataput_harp
            
        ENDFOR
        NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        
        ;alt, lat, lon
        pace_vars=['altitude','latitude','longitude']
        RSP_vars=['COLLOCATED_ALTITUDE','COLLOCATED_LATITUDE','COLLOCATED_LONGITUDE']
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
            NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR
        
        
        
        ;viewing geometry and scattering angle and solar geometry 
        pace_vars=['sensor_azimuth','sensor_zenith',$ 
            'scattering_angle',$    
            'solar_azimuth','solar_zenith']
        RSP_vars=['VIEWING_AZIMUTH','VIEWING_ZENITH',$ 
            'SCATTERING_ANGLE',$
            'SOLAR_AZIMUTH','SOLAR_ZENITH']
        convert2=[1.,-1.,1.,1.,1.]
        flip180=[0.,1.,0.,0.,0.]

        nmap=n_elements(pace_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
        
            rsp_var_map=where(TAG_NAMES(data_RSP.GEOMETRY) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
             
             dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ipix=0,bins_along_track-1 DO BEGIN

                dataput1=(((flip180[imap]*180.)+convert2[imap]*data_RSP.GEOMETRY.(rsp_var_map)._data[istart:iend,ipix,0])-Add_offset[ivar])/scale_factor[ivar]
                dataput_harp=dblarr(number_of_views)
                i1=0
                for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                    views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                    i2=i1+number_of_views_wl-1
                    dataput_harp[i1:i2]=dataput1[views_select_wl]
                    i1=i2+1
                endfor

                FOR ibin_accross_track=0,bins_across_track-1 DO dataput[*,ibin_accross_track,ipix]=dataput_harp
            ENDFOR
             NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR

        ;rotation factors
        pace_vars=['cos_rot_scatt_plane','sin_rot_scatt_plane']
        RSP_vars=['COS_ROT_SCATT_PLANE','SIN_ROT_SCATT_PLANE']
        nmap=n_elements(pace_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
        
            rsp_var_map=where(TAG_NAMES(data_RSP.GEOMETRY) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
             
             dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ipix=0,bins_along_track-1 DO BEGIN

                dataput1=data_RSP.GEOMETRY.(rsp_var_map)._data[istart:iend,ipix,0]
                dataput_harp=dblarr(number_of_views)
                i1=0
                for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                    views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                    i2=i1+number_of_views_wl-1
                    dataput_harp[i1:i2]=dataput1[views_select_wl]
                    i1=i2+1
                endfor

                FOR ibin_accross_track=0,bins_across_track-1 DO dataput[*,ibin_accross_track,ipix]=dataput_harp
            ENDFOR
             NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR

        
        ;intensity_wavelengths and polarization_wavelengths        
        pace_vars=['intensity_wavelengths','polarization_wavelengths']
        RSP_vars=['WAVELENGTH','WAVELENGTH']
        nmap=n_elements(pace_vars)
        
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
            rsp_var_map=where(TAG_NAMES(data_RSP.DATA) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]
            
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            dataput_harp=dblarr(number_of_views)

            i1=0
            for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                i2=i1+number_of_views_wl-1    
                dataput_harp[i1:i2]=data_RSP.DATA.(rsp_var_map)._data[ibandselection[iwl]]
                i1=i2+1
            endfor
            dataput[0,*]=dataput_harp

            NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR


        ;intensity__bandpasses and polarization__bandpasses (FWHM?)
        pace_vars=['intensity_bandpasses','polarization_bandpasses']
        nmap=n_elements(pace_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)

            dataput_harp=dblarr(number_of_views)
            i1=0
            for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                i2=i1+number_of_views_wl-1    
                dataput_harp[i1:i2]=bands_FWHM[ibandselection[iwl]]
                i1=i2+1
            endfor
            dataput[0,*]=dataput_harp
        
            NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR
        
        
        ;solar F0
        pace_vars=['intensity_F0','polarization_F0']
        RSP_vars=['SOLAR_CONSTANT','SOLAR_CONSTANT']
        nmap=n_elements(pace_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
            rsp_var_map=where(TAG_NAMES(data_RSP.CALIBRATION) eq RSP_vars[imap])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            IF(PACE_HARP2_L1C_vars.(pointer_ndim)[ivar] gt 1)THEN $
                FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            
            dataput_harp=dblarr(number_of_views)
            data1=(data_RSP.CALIBRATION.(rsp_var_map)._data-Add_offset[ivar])/scale_factor[ivar]
                        i1=0
            for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                i2=i1+number_of_views_wl-1    
                dataput_harp[i1:i2]=data1[ibandselection[iwl]]
                i1=i2+1
            endfor
            dataput[0,*]=dataput_harp

            NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput

        ENDFOR
    
        ;obs_per_view
        ; just fill with 1
        pace_vars='obs_per_view'
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars)
        ivar=ivar[0] 
        dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
        FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
        dataput=MAKE_ARRAY(dim_ivar,/FLOAT,value=1)
        NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        
        ;intensity
        pace_vars='I'
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars)
        ivar=ivar[0]         
        dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
        FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
        dataput=MAKE_ARRAY(dim_ivar,/FLOAT,value=1)
        FOR ipix=0,bins_along_track-1 DO BEGIN
                istart=data_RSP.Data.Unvignetted_Sector_begin._data
                iend=data_RSP.Data.Unvignetted_Sector_end._data
                I_avg=((data_RSP.data.intensity_1._data[*,istart:iend,ipix]+data_RSP.data.intensity_2._data[*,istart:iend,ipix])/2.)
                ifill=where(I_avg lt -900. or ~FINITE(I_avg),nfill)
                For iwl=0,data_RSP.dim_bands._nelements-1 DO I_avg[iwl,*]=(I_avg[iwl,*]*data_RSP.calibration.solar_constant._data[iwl]/!PI-Add_offset[ivar])/scale_factor[ivar]
                IF(nfill gt 0)THEN I_avg[ifill]=-999

                dataput_harp=dblarr(number_of_views)
                
                i1=0
                for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                    views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                    i2=i1+number_of_views_wl-1    
                    dataput_harp[i1:i2]=I_avg[ibandselection[iwl],views_select_wl]
                    i1=i2+1
                endfor
                
                dataput[0,*]=dataput_harp
                FOR ibin_accross_track=0,bins_across_track-1 DO dataput[0,*,ibin_accross_track,ipix]=dataput_harp
        ENDFOR
        NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        
        ;Q and U
        pace_vars=['Q','U']
        RSP_vars=['STOKES_Q','STOKES_U']
        nmap=n_elements(RSP_vars)
        FOR imap=0,nmap-1 DO BEGIN
            ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars[imap])
            ivar=ivar[0] 
            rsp_var_map=where(TAG_NAMES(data_RSP.DATA) eq RSP_vars[imap])
        
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT,value=1)
            FOR ipix=0,bins_along_track-1 DO BEGIN
                    istart=data_RSP.Data.Unvignetted_Sector_begin._data
                    iend=data_RSP.Data.Unvignetted_Sector_end._data
                    QU=data_RSP.data.(rsp_var_map)._data[*,istart:iend,ipix]
                    ifill=where(QU lt -900. or ~FINITE(QU),nfill)
                    For iwl=0,data_RSP.dim_bands._nelements-1 DO QU[iwl,*]=(QU[iwl,*]*data_RSP.calibration.solar_constant._data[iwl]/!PI-Add_offset[ivar])/scale_factor[ivar]
                    IF(nfill gt 0)THEN QU[ifill]=-999

                    dataput_harp=dblarr(number_of_views)

                    i1=0
                    for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                        views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                        i2=i1+number_of_views_wl-1    
                        dataput_harp[i1:i2]=QU[ibandselection[iwl],views_select_wl]
                        i1=i2+1
                    endfor

                    dataput[0,*]=dataput_harp
                    FOR ibin_accross_track=0,bins_across_track-1 DO dataput[0,*,ibin_accross_track,ipix]=dataput_harp

            ENDFOR
            NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput
        ENDFOR
        
        ;DoLP
        pace_vars='DOLP'
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_vars)
        ivar=ivar[0] 
        dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
        FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
        dataput=MAKE_ARRAY(dim_ivar,/FLOAT,value=1)
        FOR ipix=0,bins_along_track-1 DO BEGIN
        ;        Nadir_index=data_RSP.Geometry.Nadir_index._data[ipix,0]
                istart=data_RSP.Data.Unvignetted_Sector_begin._data
                iend=data_RSP.Data.Unvignetted_Sector_end._data
                Dolp1=data_RSP.Data.DOLP._data[*,istart:iend,ipix]
                ifill=where(Dolp1 lt -900. or ~FINITE(Dolp1),nfill)
                Dolp2=(Dolp1/100.-Add_offset[ivar])/scale_factor[ivar]
                IF(nfill gt 0)THEN Dolp2[ifill]=-999

                dataput_harp=dblarr(number_of_views)
                i1=0
                for iwl=0,number_of_rsp_bands-1 DO BEGIN 
                    views_select_wl=where(views_select[iwl,*] eq 1,number_of_views_wl) 
                    i2=i1+number_of_views_wl-1    
                    dataput_harp[i1:i2]=Dolp2[ibandselection[iwl],views_select_wl]
                    i1=i2+1
                endfor
                dataput[0,*]=dataput_harp
                FOR ibin_accross_track=0,bins_across_track-1 DO dataput[0,*,ibin_accross_track,ipix]=dataput_harp

        ENDFOR
        NCDF_VARPUT,group_id[PACE_HARP2_L1C_vars.(pointer_folder)[ivar]],var_id[ivar],dataput

        ; Skipped variables:
        ;0 Viewing angles at instrument 
        ;altitude_variability
        ;19 I_noise
        ;21 Q_noise
        ;23 U_noise
        ;25 DOLP_noise
        
        IF(switch_add_RSP_variables)THEN BEGIN
        ;-----RSP variables
        group_id_rsp=NCDF_GROUPDEF(id, 'RSP_specific_variables')
        
        NCDF_ATTPUT,group_id_rsp,/global,'RSP_original_file',file_RSP
        NCDF_ATTPUT,group_id_rsp,/global,'RSP_experiment',data_RSP.EXPERIMENT._data      
        NCDF_ATTPUT,group_id_rsp,/global,'RSP_aircraft',data_RSP.PLATFORM_DESCRIPTION._data
        NCDF_ATTPUT,group_id_rsp,/global,'RSP_original_read_me',data_RSP.READ_ME._data
         
        ; these are variables only as function of scan number
        RSP_folder=['DATA','DATA','DATA','DATA','DATA',$
            'PLATFORM','PLATFORM','PLATFORM']
        RSP_vars=['CLOUD_TEST_PASSED','CLOUD_TEST_PERFORMED','DATA_QUALITY_FLAGS','WATER_VAPOR_POL','WATER_VAPOR_TOTAL',$
            'PLATFORM_ALTITUDE','LAND_WATER_MASK_AT_NADIR','GROUND_ELEVATION_AT_NADIR']
        NVARIABLES=n_elements(rsp_vars)
        var_id_rsp=LONARR(NVARIABLES)
        
        pace_dimensions='altitude'; used to get dimensions...
        ivar=where(PACE_HARP2_L1C_vars.(0) eq pace_dimensions)
        ivar=ivar[0]        
        ndims=PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]         
        dim_var=INTARR(ndims)  
        FOR idim=0,ndims-1 DO dim_var[idim]=dim_id[PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar]]
        
        FOR imap=0,NVARIABLES-1 DO BEGIN
            rsp_fol_map=where(TAG_NAMES(data_RSP) eq RSP_folder[imap])
            rsp_var_map=where(TAG_NAMES(data_RSP.(rsp_fol_map)) eq RSP_vars[imap])
            
            type=TYPENAME(data_RSP.(rsp_fol_map).(rsp_var_map)._data)
            typeswitch=[1,0]
            IF(Type eq 'BYTE')THEN typeswitch=[0,1]
           
            var_id_rsp[imap]=NCDF_VARDEF(group_id_rsp,data_RSP.(rsp_fol_map).(rsp_var_map)._name,dim_var,FLOAT=typeswitch[0],BYTE=typeswitch[1])
            dim_ivar=dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1)[ivar])
            FOR idim=1,PACE_HARP2_L1C_vars.(pointer_ndim)[ivar]-1 DO dim_ivar=[dim_ivar,dimensions_rsp.(PACE_HARP2_L1C_vars.(pointer_dim1+idim)[ivar])]        
            dataput=MAKE_ARRAY(dim_ivar,FLOAT=typeswitch[0],BYTE=typeswitch[1])            
            
            nattributes_rsp=N_TAGS(data_RSP.(rsp_fol_map).(rsp_var_map))
            FOR iattr=0,nattributes_rsp-1 DO BEGIN
                IF(TYPENAME(data_RSP.(rsp_fol_map).(rsp_var_map).(iattr)) ne 'ANONYMOUS')THEN CONTINUE
                IF(data_RSP.(rsp_fol_map).(rsp_var_map).(iattr)._name eq 'DIMENSION_LIST')THEN CONTINUE
                NCDF_ATTPUT,group_id_rsp,var_id_rsp[imap],data_RSP.(rsp_fol_map).(rsp_var_map).(iattr)._name,data_RSP.(rsp_fol_map).(rsp_var_map).(iattr)._data
            ENDFOR
            dataput=MAKE_ARRAY(dim_ivar,/FLOAT)
            FOR ibin_accross_track=0,bins_across_track-1 DO dataput[ibin_accross_track,*]=REFORM(data_RSP.(rsp_fol_map).(rsp_var_map)._data)
            NCDF_VARPUT,group_id_rsp,var_id_rsp[imap],dataput


        ENDFOR
        ENDIF
        IF ~KEYWORD_SET(switch_no_talk)THEN print,'Closing output file'
        ;ncdf_control, id, /endef
        NCDF_CLOSE, id ;

        IF KEYWORD_SET(switch_stop_after_one)THEN stop
    ENDFOR;files
ENDFOR;dates

end
