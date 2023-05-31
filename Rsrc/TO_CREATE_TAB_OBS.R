### TO FILL TAB_OBS ####

TAB_OBS_ALL_SITES = read_excel("myData/TABELLA_OBS.xlsx")
site_name = 'Eksharad'

###################3

##############
site_file = "C:/Users/39348/OneDrive/Desktop/LAVORO_CHECCO/CALIBRAZIONE/Allometry/Allometry_tab_Eksharad.xlsx"
all_site_data = read_excel(site_file, sheet = 'x')

remaining_site_data = all_site_data [,c(1,2,7,9,27,28,30,31,51:53,29) ]
setnames(remaining_site_data,c('Site_name','Site_ID','Plot_ID','n_month','dbh','stem_n','height','basal_area','volume','biom_stem','biom_root','biom_foliage'))

total_site_data = all_site_data [,c(1,2,7,9,23:26,54:56)]
setnames(total_site_data,c('Site_name','Site_ID','Plot_ID','n_month','dbh','stem_n','basal_area','volume','biom_stem','biom_root','biom_foliage'))

removed_site_data = all_site_data [,c(1,2,7,9,32:35,57:59)]
setnames(removed_site_data,c('Site_name','Site_ID','Plot_ID','n_month','dbh','stem_n','basal_area','volume','biom_stem','biom_root','biom_foliage'))

dead_site_data = all_site_data [,c(1,2,7,9,42:45)]
setnames(dead_site_data,c('Site_name','Site_ID','Plot_ID','n_month','dbh','stem_n','basal_area','volume'))
########
total_site_data$height = NA
removed_site_data$height = NA
dead_site_data$height = NA
dead_site_data$biom_stem = NA 
dead_site_data$biom_root = NA
dead_site_data$biom_foliage = NA 



#########
remaining_site_data$data_type = 'remaining'
total_site_data$data_type = 'total'
removed_site_data$data_type = 'removed'
dead_site_data$data_type = 'dead'

valueX = c('dbh','stem_n','height','basal_area','volume','biom_stem','biom_root','biom_foliage')
sitex_data = rbind(remaining_site_data,total_site_data,removed_site_data,dead_site_data)


TAB_OBS_SITEX = data.table()

for (v in valueX) {
  g = which(i_output[,3]==i) #variable row
  j = i_output[g,c(1,2)] #####group and variable ID####
  TAB_OBS_SITEX = rbind(TAB_OBS_SITEX,data.table(site_id = sitex_data$Site_name,                                       ####unique(all_site_data[which(all_site_data[,1]==site_name),2]),
             Plot_ID = sitex_data$Plot_ID,
             n_month = sitex_data$n_month,
             layers_id = 1,
             group = j[1] ,
             variable = j[2] ,
             obs = sitex_data[,which(variable.names(sitex_data)== v)],
             var_name = v,
             data_type = sitex_data[,13]), fill = T)
}

TAB_OBS_SITEX  
                         
######## ########################################################################################
###########################################
#############################
site_name = 'Eksharad'
site_file = "C:/Users/39348/OneDrive/Desktop/LAVORO_CHECCO/CALIBRAZIONE/Allometry/Allometry_tab_Eksharad.xlsx"
all_site_data = read_excel(site_file, sheet = 'x')


remaining_site_data = all_site_data [,c(which(variable.names(all_site_data)== 'Site_name'),which(variable.names(all_site_data)== 'Site_ID'),which(variable.names(all_site_data)== 'Plot_ID'),which(variable.names(all_site_data)== 'n_mese'),27:31,51:53,29) ]
setnames(remaining_site_data,c('Site_name','Site_ID','Plot_ID','n_month','','','','','','','',''))
remaining_site_data$data_type = 'remaining'
all_site_data$
total_site_data = all_site_data [,c(1,2,7,9,23:26,54:56)]
total_site_data$height = NA
setnames(total_site_data,c('Site_name','Site_ID','Plot_ID','n_month','','','','','','','',''))
total_site_data$data_type = 'total'

removed_site_data = all_site_data [,c(1,2,7,9,32:35,57:59)]
removed_site_data$height = NA
setnames(removed_site_data,c('Site_name','Site_ID','Plot_ID','n_month','','','','','','','',''))
removed_site_data$data_type = 'removed'

dead_site_data = all_site_data [,c(1,2,7,9,42:45)]
dead_site_data$height = NA
dead_site_data$biom_stem = NA 
dead_site_data$biom_root = NA
dead_site_data$biom_foliage = NA 
setnames(dead_site_data,c('Site_name','Site_ID','Plot_ID','n_month','','','','','','','',''))
dead_site_data$data_type = 'dead'

#########
valueX = c('dbh','stem_n','height','basal_area','volume','biom_stem','biom_root','biom_foliage')
sitex_data = rbind(remaining_site_data,total_site_data,removed_site_data,dead_site_data)
setnames(sitex_data,c('Site_name','Site_ID','Plot_ID','n_month','','','','','','','','','data_type'))
o= c('Site_name','Site_ID','Plot_ID','n_month','dbh','stem_n','basal_area','volume','height','biom_stem','biom_root','biom_foliage','data_type')


TAB_OBS_SITEX = data.table()
for (v in valueX) {
g = which(i_output[,3]==v) #variable row
j = i_output[g,c(1,2)] #####group and variable ID####
TAB_OBS_SITEX = rbind(TAB_OBS_SITEX,data.table(site_id = sitex_data$Site_ID,                                       ####unique(all_site_data[which(all_site_data[,1]==site_name),2]),
                                               Plot_ID = sitex_data$Plot_ID,
                                               n_month = sitex_data$n_month,
                                               layers_id = 1,
                                               group = j[1] ,
                                               variable = j[2] ,
                                               obs = sitex_data[,which(o== v)],
                                               var_name = v,
                                               data_type = sitex_data[,13]), fill = T)
}

setnames(TAB_OBS_SITEX,c('site_id','Plot_ID','n_month','layers_id','group','variable','obs','var_name','data_type'))
TAB_OBS_ALL_SITES = rbind(TAB_OBS_ALL_SITES,TAB_OBS_SITEX)

##
writexl::write_xlsx(TAB_OBS_ALL_SITES,"myData/TABELLA_OBS.xlsx")




