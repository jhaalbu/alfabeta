#TODO: Betapunkt/linje definert av bruker. Bruke arcgis til å finne skjæringspunkt?? Problem med å treffe meterverdi? Nærmeste punkt??
#TODO: Rekne på profil, og ikkje polylinje? Eller berre bruke høggrads polynom?

import arcpy
import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

class Profil:
    def __init__(self, inputfc, terreng, tempfc_punkter):
        self.inputfc = inputfc
        print(inputfc)
        self.terreng = terreng
        print(terreng)
        self.tempfc_punkter = tempfc_punkter
        print(tempfc_punkter)

        #Langer punkt kvar 1 meter langs input profil
        arcpy.GeneratePointsAlongLines_management(inputfc, tempfc_punkter, 'DISTANCE',
                                            Distance='1 Meters')
        #Legger z koordinater på punkter ut frå valt raster surface
        arcpy.ddd.AddSurfaceInformation(tempfc_punkter, terreng, 'Z', 'BILINEAR')
        #Tar ut koordinatlister frå puntkter feature class 
        with arcpy.da.SearchCursor(tempfc_punkter, ["SHAPE", 'Z']) as cursor:
            x_list = []
            y_list = []
            z_list = []
            for row in cursor:
                x, y = row[0]
                x_list.append(x)
                y_list.append(y)
                z_list.append(row[1])
        #Etablerer Pandas dataframe for forenkling av vidare databehandling
        self.df = pd.DataFrame(list(zip(x_list, y_list, z_list)), columns =['X', 'Y', 'Z'])

            #Regner ut distanse mellom punkter (muligens unødvening, kan kanskje bruke index
        #til punkter istadenfor? Sidan kvar punkt er etalbert per meter?.
        self.df['DIST'] = np.sqrt(((self.df.X - self.df.X.shift(1))**2)+((self.df.Y - self.df.Y.shift(1))**2))
        self.df['M'] = self.df.DIST + self.df.DIST.shift(1)
        self.df.loc[0, 'M'] = 0
        self.df.loc[0, 'DIST'] = 0
        self.df.loc[1, 'M'] = self.df.loc[1, 'DIST']
        self.df.loc[0, 'H'] = 0
        
        #Regner ut lengden basert på avstand mellom punkter
        for i in range(2, len(self.df)):
            self.df.loc[i, 'M'] = self.df.loc[i-1, 'DIST'] + self.df.loc[i-1, 'M']

        #Runder av meterverdien 
        self.df['M'] = self.df['M'].round(0)

    def poly(self, polynom=2):
        '''
        input
        df: dataframe med profil, fra profil funksjonen
        polynom: størrelseorden på polynomet, standard 2. grads
        
        output
        returnerer ein dataframe med kolonnene POLY og H_DEG, som er tilpassa profil og helling i grader
        '''
        self.polynom = polynom
        #Tar meterverdi og høgdeverdi for meterverdi (altså høgdeprofilet) og regner
        #om til ein tilpasse polynom (forenkler geometrien)
        self.p = np.poly1d(np.polyfit(self.df.index, self.df['Z'], self.polynom))
        
        #Etablerer pandas kollonen poly for representerer det forenkla polynomet
        self.df['POLY'] = self.p(self.df.index)

        #Regner ut hellingen langs det forenkla profilet
        for i in range(1, len(self.df)):
            self.df.loc[i, 'H'] = ((self.df.loc[i, 'POLY'] - self.df.loc[i -1, 'POLY'])/(self.df.loc[i, 'M'] - self.df.loc[i - 1, 'M']))

        #Regner om frå hellingstall til vinkel
        self.df['H_DEG'] = np.degrees(np.arctan(self.df['H']))
        
    def get_profil(self):
        return self.df        

    def plot_profil(self):
        plt.plot(self.df['M'], self.df['Z'])
        plt.title(f'Høydeprofil {self.inputfc}')
        plt.xlabel("Lengde (m)")
        plt.ylabel("Høyde (m)")
        plt.grid(True)
        plt.show()

class Skred:
    '''
    Eit skredobjekt kan berre ha ein skredtype
    Profil input må vere ei dataframe frå Profil objekt
    Fra get_profil metoden
    '''

    def __init__(self, profil, skredtype):
        self.profil = profil
        self.df = profil.df
        self.x_start = self.df['X'][0]
        self.y_start = self.df['Y'][0]
        self.m_start = 0
        self.z_Start = self.df['Z'][0]
        self.skredtype = skredtype

    #def betapunkt():
        #Korleis håndtere betapunkt definert av bruker?
        #Etablert metode for å kunne gjere dette etter kvart
        betavinkler = {'sno':-10,'stein':-23,'jord':-20}

        avik = 0.2
        #print(self.df)
        print('beta vinkel', betavinkler[skredtype])
        df10 = self.df.loc[(self.df['H_DEG'] <= (betavinkler[skredtype]+avik)) & (self.df['H_DEG'] >= (betavinkler[skredtype]-avik))]
        index_label = self.df[(self.df['H_DEG'] <= (betavinkler[skredtype]+avik)) & (self.df['H_DEG'] >= (betavinkler[skredtype]-avik))].index.tolist()
        print(index_label)
        skjeringspunkt = index_label[0]

        #Finner "koordinater" for 10 graderspunktet og topp punkt av grafen for å regne beta vinkel

        m_topp = self.df.loc[0, 'M']
        poly_topp = self.df.loc[0, 'POLY']
        m_10 = self.df.loc[skjeringspunkt, 'M']
        poly_10 = self.df.loc[skjeringspunkt, 'POLY']

        #Koordinatpunkt for betapunkt
        self.x_beta = self.df.loc[m_10, 'X']
        self.y_beta = self.df.loc[m_10, 'Y']

        #Regner ut beta vinkelen
        self.beta_helning = abs((poly_10 - poly_topp)/(m_10 - m_topp))
        self.beta_vinkel_radianer = np.arctan(self.beta_helning)
        self.beta_vinkel_grader = np.degrees(self.beta_vinkel_radianer)  

        self.m_beta = m_10
        self.z_beta = poly_10

    def runout(self, sigma=5):
        skredtype = self.skredtype
        alfaparameter = {'sno':[2.3, 0.96, 1.4], 'stein':[2.16, 0.77, -3.9], 'jord':[1.5, 0.96, 4.0]}
        sd = alfaparameter[self.skredtype][0]
        tilpassing = alfaparameter[self.skredtype][1]
        justering = alfaparameter[self.skredtype][2]
        print(f'tilpassing er {tilpassing}, justering er {justering}, standardavik er {sd}, betavinkel er {self.beta_vinkel_grader}')
        self.alfa_vinkelliste = []
        self.alfa_hellingsliste = []
        for i in range(5):
            alfa_vinkel = (tilpassing * self.beta_vinkel_grader) - justering  - (sd * i)
            alfa_helning = np.tan(np.radians(alfa_vinkel))
            print(f'alfavinkelen er {alfa_vinkel}')
            self.alfa_vinkelliste.append(alfa_vinkel)
            self.alfa_hellingsliste.append(alfa_helning)
        print(self.alfa_vinkelliste)
        liste_meterverdi = [0, self.df.loc[len(self.df)-1, 'M']]
        sigmateller = 0
        self.alfa_koordinater = []
        self.alfa_plotverdier = []
        
        while True:
            try:
                print('try' + str(sigmateller))
                liste_alfa = [self.df.loc[0, 'POLY'], self.df.loc[0, 'POLY'] - self.df.loc[len(self.df)-1, 'M']*self.alfa_hellingsliste[sigmateller]]
                #Finer rotpunktet mellom skredbanen (p) (andregradspolynom) og alfa vinkelplanet (q)
                q_alfa = np.polyfit(liste_meterverdi, liste_alfa, 1)
                x_0_alfa = np.roots(self.profil.p - q_alfa)
                alfa_verdi = int(x_0_alfa.max().round())
                alfa_utlop_x = self.df.loc[alfa_verdi, 'X']
                alfa_utlop_y = self.df.loc[alfa_verdi, 'Y']
                alfa_m = self.df.loc[alfa_verdi, 'M']
                alfa_poly = self.df.loc[alfa_verdi, 'POLY']
                self.alfa_koordinater.append((alfa_utlop_x, alfa_utlop_y))
                self.alfa_plotverdier.append((alfa_m, alfa_poly))
            except:
                break
            sigmateller += 1

        return self.alfa_koordinater, self.alfa_plotverdier

def lag_featurepunkt(feature, fgdb, liste_med_punkter):
    fc = profilnavn + '_alfabeta_'+skredtype
    arcpy.CreateFeatureclass_management(fgdb, fc, "Point", "", "", "", 25833)
    arcpy.AddField_management(fc, 'Value', 'NAVN', "TEXT")

    with arcpy.da.InsertCursor(fc, ["NAVN", "SHAPE@"]) as cursor:
        cursor.insertRow(('Alfa', (alfa[0][0], alfa[0][1])))
                             
    #features = [fc_alfa, fc_beta, fc_sigma1, fc_sigma2]
    
    data = fgdb + "\\" + fc
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    aprxMap = aprx.listMaps(aprx.activeMap.name)[0] 
    aprxMap.addDataFromPath(data)

def plot_alfa(profil, skred):
    df = profil.df
    
    # arcpy.AddMessage(str(alfa))
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(df['M'], df['Z'], label='Høgdeprofil') #Høgdeprofilet
    ax.plot(df['M'], df['POLY'], label=f'Tilpasset profil {profil.polynom}. grads') #Forenkla høgdeprofil
    #ax.scatter(beta[5], beta[6], color='r', linewidth='1', label='Punkt med 10 grader helling') # 10 graders punkter
    ax.plot([df['M'][0], skred.m_beta], [df['POLY'][0], skred.z_beta], label=f'Beta {round(skred.beta_vinkel_grader, 1)} \xb0', linestyle='--') #Beta 
    ax.plot([0, skred.alfa_plotverdier[0][0]], [df['POLY'][0], skred.alfa_plotverdier[0][1]], label=f'Alfa {round(skred.alfa_vinkelliste[0], 1)}\xb0') #Må plotte til skjæeringspunkt
    loc = ticker.MultipleLocator(base=100.0)
    ax.xaxis.set_major_locator(loc)
    ax.axis('equal')
    ax.legend()
    return ax

#Setter workspace
#fgdb = arcpy.GetParameterAsText(3)
fgdb = "C:/Users/jan.aalbu/Documents/ArcGIS/Projects/Alfabeta/Alfabeta.gdb"
arcpy.env.workspace = fgdb
arcpy.management.Delete("C:/Users/jan.aalbu/Documents/ArcGIS/Projects/Alfabeta/Alfabeta.gdb/test_punkter")
#Inputdata til beregning
#inputfc_profil = arcpy.GetParameterAsText(0)
inputfc_profil = 'Terrengprofil_Profil_Lines'
desc = arcpy.Describe(inputfc_profil)
profilnavn = desc.name
input_skredtype = 'sno'
outputfc1 = 'test_punkter'
insurface = 'C:/Users/jan.aalbu/Documents/ArcGIS/Projects/Heilevang/Data/DTM Polygon/dtm/NDH Askvoll 5pkt 2018-dtm.tif'

profil1 = Profil(inputfc_profil, insurface, outputfc1)

profil1.poly(4)
#print(profil1.df)     
#profil1.plot_profil()

snoskred = Skred(profil1, 'sno')
#snoskred.betapunkt()
print(snoskred.x_beta, snoskred.y_beta)
print(snoskred.runout())
print(snoskred.alfa_plotverdier)
plot_alfa(profil1, snoskred)
plt.show()