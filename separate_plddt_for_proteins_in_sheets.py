### Omkar Golatkar ###
### description: takes excel file as input and makes groups of the dataframe based on defined column 
### input file: replace with your filename in place of out_data.tsv; make sure that this py file is in the same directory as your input file or mention the whole path of input file
### 		the data in the excel file will be grouped based on the content in the column 'name'; replace it with the column you want to use for grouping
### output file: output file will be saved in xlsx format with name 'pdata.xlsx'; replace it with your desired filename; it will be saved in the same directory as the '.py' file

import pandas as pd
import openpyxl as op

df = pd.read_csv('out_data.tsv',sep = '\t')

data = df.groupby('name')
writer = pd.ExcelWriter('pdata.xlsx', engine='xlsxwriter')
for row,group in data:
       group.to_excel(writer, sheet_name=row)
writer._save()

wb = op.load_workbook('pdata.xlsx')
res = len(wb.sheetnames)
wb._sheets.sort(key=lambda ws:ws.title)
print('number of protein predictions from alphafold:',res)
writer._save()
