import os

if __name__ == "__main__":
    # where is the data?
    data_dir = 'C:\\Users\\Saima\\Documents\\REL Projects\\Homeira Heterogeneity\\AHP\\' #need double backslash otherwise python thinks its an escape sequence

    #what folders do the data live in?
    folder_list = ['Total_5_Homeira', 'TotaL_5_Lihua', 'Total2Homeira', 'Total2lihua']

    #where to save the figures?
    export_dir = 'ExportDir'

    #HEADINGS
    # Dir	                	Offset	Offset	AnChan	Comment	ExportDir	Tags	                                            GroupInc
    #Total5	20131211_600_1_0001f	0	0	2		                Cell 1, RMP - 68.5mv, L5, Tau 3, Gain 20, DC150	    yes

    #this is problematic if there are nested folders or have other files...
    #should check for nesting but I'm too lazy right now
    
    #some variables for writing the excel file
    offset = str(0)
    anchan = str(2)
    blank = ' '
    tags = 'None'
    groupinc = 'yes'

    comma = ','
    comment = blank
    
    excel_name = 'AHP.csv'
    f = open(excel_name, 'w')

    print('---- Writing the excel sheet, patience is a virtue ----')

    for folder in folder_list:
        files_array = os.listdir(data_dir + folder) #have a list of all the files
        for file in files_array:
            if file[-4:] != '.abf': #not an abf file, skip
                print(file + ' is not the right file type, silly goose!! It has been skipped')
                continue
            
            this_line = [folder, file[0:-4] + 'f', offset, offset, anchan, comment, export_dir, tags, groupinc]

            #just checking my string indexing
            #print(file[0:-4])
            #print(file)
            #me = input()

            f.write(comma.join(this_line) + '\n') #make it comma delimited!

    f.close() #close the file
    
    print('---- All done, friend~ ;-) ----')
