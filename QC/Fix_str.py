import pandas as pd

str_file = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/erisor.str"

# Import the .str file into pandas.
IN_STR = pd.read_csv(str_file, sep="\t", header=None)

# This file has extra columns with no info. I remove them: axis=1 to remove column, how='all' to remove columns where ALL rows are blank.
INN_STR = IN_STR.dropna(axis= 1, how= 'all')
INN_STR = IN_STR.dropna(how='all')

INN_STR.to_csv("/Users/jimblotter/Desktop/Grad_School/Data_Analysis/QC/drop.str", sep='\t',index= False,header= False)



#psuedocode:
#if, after the id, there's nothing other than -9, drop the whole row.



# na_values='-9' in the read_csv will get rid of nas but I need the whole line to be gone and stored elsewhere.

#to get distance matrix:
#read in erisor.str and import to pandas
#remove extra columns with no info
#remove lines with only -9, removing the ENTIRE line. store those lines in another file.

#once I have my matrix:
#meet with Aaron to write a crazy indexing code to get reps

#replace blank with NaN then drop them again?
#instead of dropping create another data frame that selects for them in the origional file.
#name = name.fillna('-9')