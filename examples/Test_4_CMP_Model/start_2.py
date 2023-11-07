import math

########################################################## Stack parameter
T_ba        = 100  
T_teos      = 350

########################################################## Subroutine
def d_max_P3(H,alpha,beta,width,space):
    
    return H*math.pow(width,alpha)*math.pow(space,beta)

def function_fitting(width,space,k1,k2,k3):
    density = width/(width+space)

    tmp = k1*math.exp(-math.pow(width,k2))*math.pow(math.exp(density*(1-density)),k3)

    return tmp

class P3_Model:
    def __init__(
        self, density,
        H_Cu, H_Ox,
        Blanket_Cu, Blanket_BA, Blanket_TEOS, Blanket_BD,
        D_max, K_U, K_D,
        delta, ddl
    ):
        self.Cu_density = density
        
        self.H_Cu = H_Cu
        self.H_Ox = H_Ox
        self.SH   = H_Cu - H_Ox

        self.Blanket_Cu     = Blanket_Cu
        self.Blanket_BA     = Blanket_BA
        self.Blanket_TEOS   = Blanket_TEOS
        self.Blanket_BD     = Blanket_BD

        #!待解释
        self.D_max = D_max*(1-self.Cu_density)/self.Cu_density
        #self.D_max = D_max*(1-self.Cu_density)/self.Cu_density 
        self.D_p   = D_max
        self.K_U   = K_U    # 疑虑在于 K_U 和 K_D 在形貌反转上可能是否需要进行变换（先前简单的考虑并没有考虑该点）
        self.K_D   = K_D
    
        self.time       = 0
        self.step       = delta
        self.time_line  = ddl

    def check(self):
        if self.H_Ox <= T_ba:
            self.Blanket_Ox = self.Blanket_BA
        elif ((self.H_Ox > T_ba) and (self.H_Ox <= (T_ba + T_teos))):
            self.Blanket_Ox = self.Blanket_TEOS
        else:
            self.Blanket_Ox = self.Blanket_BD


    def RR_model_update(self):
        if(self.SH <= (-self.D_p)):
            self.RR_ox = 0
            self.RR_Cu = (self.K_D*self.Blanket_Cu)/self.Cu_density
        elif( (self.SH > (-self.D_p)) and (self.SH < self.D_max)):
            self.RR_ox = (self.K_U*self.Blanket_Ox)*(1+self.SH/self.D_p)
            self.RR_Cu = (self.K_D*self.Blanket_Cu)*( 1- (1-self.Cu_density)/self.Cu_density*self.SH/self.D_p)
        elif(self.SH >= (self.D_max)):
            self.RR_ox = (self.K_U*self.Blanket_Ox)/(1-self.Cu_density)
            self.RR_Cu = 0
        else:
            print("Something wrong in RR_model_update() !")

    def topography_update(self):
        self.check()
        self.RR_model_update()

        self.H_Cu += self.RR_Cu * self.step
        self.H_Ox += self.RR_ox * self.step
        self.SH   = self.H_Cu - self.H_Ox
        
        self.time += self.step

    def topography_evolution(self):
        while(self.time < self.time_line):
            self.topography_update()

        return self.time, (self.H_Cu-self.H_Ox), self.H_Ox
    
########################################################## ODE
"""
f[0] = width
f[1] = space
f[2] = H_Cu
f[3] = H_Ox
f[4] = r_Cu_P3
f[5] = r_BA_P3
f[6] = r_TEOS_P3
f[7] = r_BD_P3
"""

def derivative(f):