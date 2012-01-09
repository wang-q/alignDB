Set objExcel = CreateObject("Excel.Application")
Set objWorkbook = objExcel.Workbooks.Open(Wscript.Arguments(0))
objExcel.Application.Visible = False
objExcel.Application.DisplayAlerts = False 
objExcel.ActiveWorkbook.SaveAs Wscript.Arguments(1), 43
objExcel.ActiveWorkbook.Close
objExcel.Application.DisplayAlerts = True
objExcel.Application.Quit
WScript.Quit
' XLSX->XLS conversion script by Michael Kuzmin
' http://kuzmin.ca/blog/?p=578
' This is how you can you it:
' c:\xlsx2xls.vbs "C:\inputdocument.xlsx" "C:\outputdocument.xls"