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

def Theta_1500t_Saxton2006(S,C,OM):

 return -0.024*S + 0.487*C + 0.006*OM + 0.005*S*OM - 0.013*C*OM + 0.068*S*C + 0.031

def Theta_1500_Saxton2006(S,C,OM):

 Theta1500t = Theta_1500t_Saxton2006(S,C,OM)

 return Theta1500t + (0.14*Theta1500t - 0.02)

def Theta_33t_Saxton2006(S,C,OM):
 
 return -0.251*S + 0.195*C + 0.011*OM + 0.006*S*OM - 0.027*C*OM + 0.452*S*C + 0.299

def Theta_33_Saxton2006(S,C,OM):

 Theta33t = Theta_33t_Saxton2006(S,C,OM)

 return Theta33t + 1.283*Theta33t**2 - 0.374*Theta33t - 0.015

def Theta_S33t_Saxton2006(S,C,OM):

 return 0.278*S + 0.034*C + 0.022*OM - 0.018*S*OM - 0.027*C*OM - 0.584*S*C + 0.078

def Theta_S33_Saxton2006(S,C,OM):

 Theta_S33t = Theta_S33t_Saxton2006(S,C,OM)

 return Theta_S33t + 0.636*Theta_S33t - 0.107

def ThetaS_Saxton2006(S,C,OM):

 Theta33 = Theta_33_Saxton2006(S,C,OM) 
 ThetaS33 = Theta_S33_Saxton2006(S,C,OM)

 return Theta33 + ThetaS33 - 0.097*S + 0.043

def FAO_Soil_Texture(S,C,ST):

 ##Coarse  Medium   Fine    CM     CF     MF    CMF
 ids = np.arange(1,8)
 msand = np.array([83,37,17,60,50,27,46])
 mclay = np.array([9,30,67,20,38,48,35])
 msilt = np.array([8,33,17,20,12,25,19])

 #Calculate the closest texture class
 m = (S.mask == 0) & (C.mask == 0) & (ST.mask == 0)
 dsand,dclay,dsilt = [],[],[]
 for i in xrange(7):
  dsand.append(S[m] - msand[i])
  dclay.append(C[m] - mclay[i])
  dsilt.append(ST[m] - msilt[i])
 dsand = np.array(dsand)
 dclay = np.array(dclay)
 dsilt = np.array(dsilt)
 
 #Calculate the euclidean distance
 ed = (dsand**2 + dsilt**2 + dclay**2)**0.5

 #Determine the id
 tmp = np.argmin(ed,axis=0) + 1
 tclass = np.zeros(S.shape)
 tclass[:] = -9999.0
 tclass[m] = tmp
 tclass = np.ma.masked_array(tclass,tclass==-9999)

 return tclass
