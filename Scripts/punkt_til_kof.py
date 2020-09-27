# Python script: punkt_til_kof.py
# Author: Jan Helge Aalbu, Asplan Viak
# Scriptet tar inn ein feature class med punkter, og skrive ut ei .kof fil med
# navnet <featureclass navn>_export.kof

import arcpy
print("Arcpy importert")
inputfc = arcpy.GetParameterAsText(0)
output_folder = arcpy.GetParameterAsText(1)
#Hardkoder inn typenr og sosikoden
typenr = '05'
sosikode = '2200'


with arcpy.da.SearchCursor(inputfc, ["SHAPE@XY", "Z", "OBJECTID"]) as cursor:
    with open(output_folder + '\\' + inputfc + '_export.kof', 'w') as f: #Setter navn på fil til navn paa featureclass + _export.kof
        for row in cursor:      
            x,y = row[0]
            z = row[1]
            punktnr = row[2]
            print(f"{typenr:^4}{punktnr:<11}{sosikode:<12}{y:<13.3f} {x:<13.3f} {z:0.2f}", file=f) #printer linje for linje til fil
