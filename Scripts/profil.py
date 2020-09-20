import arcpy
import sys
import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


def profil(inputfc, terreng, outputfc='punkter'):
    #Langer punkt kvar 1 meter langs input profil
    arcpy.GeneratePointsAlongLines_management(inputfc, outputfc, 'DISTANCE',
                                          Distance='1 Meters')
    #Legger z koordinater på punkter ut frå valt raster surface
    arcpy.ddd.AddSurfaceInformation(outputfc, terreng, 'Z', 'BILINEAR')
    #Tar ut koordinatlister frå puntkter feature class 
    with arcpy.da.SearchCursor(outputfc, ["SHAPE", 'Z']) as cursor:
        x_list = []
        y_list = []
        z_list = []
        for row in cursor:
            x, y = row[0]
            x_list.append(x)
            y_list.append(y)
            z_list.append(row[1])
    #Etablerer Pandas dataframe for forenkling av vidare databehandling
    df = pd.DataFrame(list(zip(x_list, y_list, z_list)), columns =['X', 'Y', 'Z'])
    
    #Regner ut distanse mellom punkter (muligens unødvening, kan kanskje bruke index
    #til punkter istadenfor? Sidan kvar punkt er etalbert per meter?.
    df['DIST'] = np.sqrt(((df.X - df.X.shift(1))**2)+((df.Y - df.Y.shift(1))**2))
    df['M'] = df.DIST + df.DIST.shift(1)
    df.loc[0, 'M'] = 0
    df.loc[0, 'DIST'] = 0
    df.loc[1, 'M'] = df.loc[1, 'DIST']
    df.loc[0, 'H'] = 0
    
    #Regner ut lengden basert på avstand mellom punkter
    for i in range(2, len(df)):
        df.loc[i, 'M'] = df.loc[i-1, 'DIST'] + df.loc[i-1, 'M']

    #Runder av meterverdien 
    df['M'] = df['M'].round(0)
    
    
    return df

#fgdb = arcpy.GetParameterAsText(2)
#arcpy.env.workspace = fgdb

#Inputdata til beregning
inputfc_profil = arcpy.GetParameterAsText(0)
insurface = arcpy.GetParameterAsText(1)

df = profil(inputfc_profil, insurface)

plt.plot(df['M'], df['Z'])
plt.title(f'Høydeprofil {inputfc_profil}')
plt.xlabel("Lengde (m)")
plt.ylabel("Høyde (m)")
plt.grid(True)
plt.show()

