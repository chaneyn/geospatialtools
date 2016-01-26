import numpy as np

def Lambda_Maidment92(phi,clay,sand):

 #Maidment,1992
 Lambda = np.exp(-0.7842831 + 0.0177544*sand - 1.062498*phi - 0.00005304*sand*sand - 0.00273493*clay*clay + 1.11134946*phi*phi - 0.03088295*sand*phi + 0.00026587*sand*sand*phi*phi - 0.00610522*clay*clay*phi*phi - 0.00000235*sand*sand*clay + 0.00798746*clay*clay*phi - 0.00674491*phi*phi*clay) 
 return Lambda

def Residual_Water_Content_Maidment92(phi,clay,sand):

 #Maidment,1992
 SMR = -0.0182482 + 0.00087269*sand + 0.00513488*clay + 0.02939286*phi - 0.00015395*clay*clay - 0.0010827*sand*phi - 0.00018233*clay*clay*phi*phi + 0.00030703*clay*clay*phi - 0.0023584*phi*phi*clay
 return SMR

def Bubbling_Pressure_Maidment92(phi,clay,sand):

 #Maidment,1992
 Psi = np.exp(5.3396738 + 0.1845038*clay - 2.48394546*phi - 0.00213853*clay*clay 
        - 0.04356349*sand*phi - 0.61745089*clay*phi + 0.00143598*sand*sand*phi*phi
        - 0.00855375*clay*clay*phi*phi - 0.00001282*sand*sand*clay + 0.00895359*clay*clay*phi
        - 0.00072472*sand*sand*phi + 0.0000054*clay*clay*sand + 0.50028060*phi*phi*clay)
 return Psi

def Texture_Class(clay,silt,sand):
 
 texture_class = np.copy(clay)
 texture_class[:] = -9999.0
 #SAND (1) ((silt + 1.5*clay) < 15)
 mask = (silt + 1.5*clay) < 15.0
 texture_class[mask] = 1
 #LOAMY SAND (2) ((silt + 1.5*clay >= 15) && (silt + 2*clay < 30))	LOAMY SAND
 mask = ((silt + 1.5*clay >= 15) & (silt + 2*clay < 30))
 texture_class[mask] = 2
 #SANDY LOAM (3) ((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30))
 mask = (((clay >= 7) & (clay < 20) & (sand > 52) & ((silt + 2*clay) >= 30)) | ((clay < 7) & (silt < 50) & ((silt+2*clay) >=30)))
 texture_class[mask] = 3
 #SILT LOAM (4) ((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12))
 mask = ((silt >= 50) & (clay >= 12) & (clay < 27)) | ((silt >= 50) & (silt < 80) & (clay < 12))
 texture_class[mask] = 4
 #SILT (5) (silt >= 80 && clay < 12)	SILT
 mask = (silt >= 80) & (clay < 12)
 texture_class[mask] = 5
 #LOAM (6) ((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52))	LOAM
 mask = ((clay >= 7) & (clay < 27) & (silt >= 28) & (silt < 50) & (sand <= 52))
 texture_class[mask] = 6
 #SANDY CLAY LOAM (7) ((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) 	SANDY CLAY LOAM
 mask = ((clay >= 20) & (clay < 35) & (silt < 28) & (sand > 45))
 texture_class[mask] = 7
 #SILTY CLAY LOAM (8) ((clay >= 27 && clay < 40) && (sand  <= 20))	SILTY CLAY LOAM
 mask = ((clay >= 27) & (clay < 40) & (sand  <= 20))
 texture_class[mask] = 8 
 #CLAY LOAM (9) ((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45))	CLAY LOAM
 mask = ((clay >= 27) & (clay < 40) & (sand > 20) & (sand <= 45))
 texture_class[mask] = 9
 #SANDY CLAY (10) (clay >= 35 && sand > 45)	SANDY CLAY
 mask = (clay >= 35) & (sand > 45)
 texture_class[mask] = 10
 #SILTY CLAY (11) (clay >= 40 && silt >= 40)	SILTY CLAY
 mask = (clay >= 40) & (silt >= 40)
 texture_class[mask] = 11
 #CLAY (12) (clay >= 40 && sand <= 45 && silt < 40)	CLAY
 mask = (clay >= 40) & (sand <= 45) & (silt < 40)
 texture_class[mask] = 12

 return texture_class
