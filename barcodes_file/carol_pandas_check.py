import pandas as pd

# Read in the two tables
before = pd.read_csv('/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/barcodes_file/before.csv')
after = pd.read_csv('//Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/barcodes_file/after.csv')

# Making sure before table read in correctly
print(before.head())
# This will give table dimensions in row, column
print(before.shape)
# Making sure after table read in correctly
print(after.head())
# This will give table dimensions in row, column
print(after.shape)

# Merging the tables on the "Well" column. I could do this on both columns and then compare diff in length of dataframe, but think this may be better way to solve
mergedtable = before.merge(after, on="Well")
print(mergedtable.head())

# Now create a new column for True/False if the sample columns match or not
mergedtable['Match'] = mergedtable.Sample_x == mergedtable.Sample_y
# This will give table dimensions in row, column
print(mergedtable.shape)

# NOTE: Make sure all tables have the same number of rows!!!! Otherwise, since we are merging on "Well", that could mean you have a difference in "Well" column

# Then can just save file and open up and see if there are any False in the "Match" column
#mergedtable.to_csv('/Path/to/file.csv')

# Or just return the table where "Match" column contains False
FalseTable = mergedtable[mergedtable.Match == False]
print("This is your False table")
print(FalseTable)