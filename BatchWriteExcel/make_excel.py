import os

if __name__ == "__main__":
    
    print('Please refer to info.txt for instructions and tips')
    
    # ----------- MAIN DATA DIRECTORY ---------
    data_dir = 'C:\\Users\\Saima\\Documents\\REL Projects\\Homeira Heterogeneity\\AHP\\'
    #need DOUBLE backslash

    #----------- SUBFOLDERS WITH DATA ---------
    folder_list = ['Total_5_Homeira', 'TotaL_5_Lihua', 'Total2Homeira', 'Total2lihua']

    #------------ EXCEL FILE NAME ------------
    excel_name = 'AHP_04182018.csv' #this should be a csv file extension
    
    #------------ HEADING SETTINGS -----------
    blank = ' ' # can change some fields to blank if needed

    offset = str(0)
    anchan = str(2)
    tags = 'None'
    groupinc = 'yes'
    export_dir = 'ExportDir'     #where to save the figures?
    comment = blank

    '''HEADINGS
    # Dir	                	Offset	Offset	AnChan	Comment	ExportDir	Tags	                                            GroupInc
    #Total5	20131211_600_1_0001f	0	0	2		                Cell 1, RMP - 68.5mv, L5, Tau 3, Gain 20, DC150	    yes

    #this code is problematic if there are nested folders
    #should check for nesting but I'm too lazy right now
    '''

    #------------ FROM HERE BELOW, DO NOT CHANGE ANYTHING!!! ------
    headings = ['Dir', 'Offset', 'Offset', 'AnChan', 'Comment', 'ExportDir', 'Tags', 'GroupInc']
    static_settings = [offset, offset, anchan, comment, export_dir, tags, groupinc]

    comma = ','

    print('---- Writing the excel sheet, patience is a virtue ----')

    f = open(excel_name, 'w')
    f.write(comma.join(headings) + '\n')

    for folder in folder_list:
        files_array = os.listdir(data_dir + folder) #have a list of all the files
        for file in files_array:
            if file[-4:] != '.abf': #not an abf file, skip
                print(file + ' is not the right file type, silly goose!! It has been skipped')
                continue
            
            this_line = [folder, file[0:-4] + 'f'] + static_settings
            f.write(comma.join(this_line) + '\n') #make it comma delimited!

    f.close() #close the file
    
    print('---- All done, friend~ ;-) ----')
