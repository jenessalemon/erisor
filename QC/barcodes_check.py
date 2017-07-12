before = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/before"
after = "/Users/jimblotter/Desktop/Grad_School/Data_Analysis/erisor/QC/after"
import xlrd

workbook = xlrd.open_workbook(r"before.csv")
sheet = workbook.sheet_by_index(0)

col_a = sheet.col_values(0, 1)
col_b = sheet.col_values(1, 1)

my_dict = {a : b for a, b in zip(col_a, col_b)}

print (my_dict)