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

def Lambda_Saxton_2006(S,C,OM):

 Theta1500 = Theta_1500_Saxton2006(S,C,OM)
 Theta33 = Theta_33_Saxton2006(S,C,OM)
 return (np.log(Theta33) - np.log(Theta1500))/(np.log(1500) - np.log(33))

def Ksat_Saxton2006(S,C,OM):

 ThetaS = ThetaS_Saxton2006(S,C,OM)
 Theta33 = Theta_33_Saxton2006(S,C,OM)
 Lambda = Lambda_Saxton_2006(S,C,OM)
 return 1930*(ThetaS - Theta33)**(3-Lambda)

def Psisat_Saxton2006(S,C,OM):

  vwcr = 0
  b = 1/Lambda_Saxton_2006(S,C,OM)
  vwc33 = Theta_33_Saxton2006(S,C,OM)
  vwc0 = ThetaS_Saxton2006(S,C,OM)
  vwc = vwc33
  psi = 33

  return psi*((vwc-vwcr)/(vwc0-vwcr))**b

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

def Run_Tests():

 clay = np.array([88.0,80.0,65.0,40.0,20.0,10.0,60.0,30.0,10.0,10.0,50.0,25.0])/100
 sand = np.array([5.0,5.0,10.0,20.0,15.0,5.0,25.0,35.0,35.0,45.0,40.0,50.0])/100
 om = 2.5
 ksat = np.array([108.,96.7,50.3,15.5,16.1,22.,11.3,4.3,5.7,3.7,1.4,1.1])
 thetas = np.array([46.,46.,45.,46.,48.,48.,43.,48.,51.,52.,44.,50.])/100
 theta33 = np.array([10.,12.,18.,28.,31.,30.,27.,36.,38.,41.,36.,42.])/100
 theta1500 = np.array([5.,5.,8.,14.,11.,6.,17.,22.,22.,27.,25.,30.])/100

 print "Comparing ksat"
 print ksat
 print Ksat_Saxton2006(clay,sand,om)
 print np.allclose(ksat,Ksat_Saxton2006(clay,sand,om),atol=5e-01)
 print "Comparing thetas"
 print np.allclose(thetas,ThetaS_Saxton2006(clay,sand,om),atol=1e-02)
 print thetas
 print ThetaS_Saxton2006(clay,sand,om)
 print "Comparing theta33"
 print np.allclose(theta33,Theta_33_Saxton2006(clay,sand,om),atol=1e-02)
 print theta33
 print Theta_33_Saxton2006(clay,sand,om)
 print "Comparing theta1500"
 print np.allclose(theta1500,Theta_1500_Saxton2006(clay,sand,om),atol=1e-02)
 print theta1500
 print Theta_1500_Saxton2006(clay,sand,om)
 print "Psisat"
 print 10*Psisat_Saxton2006(clay,sand,om)/100 #m

 return
