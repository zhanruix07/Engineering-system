'''
Modelling of Engineering Systems
Coursework Question 8
'''
# Created on 13/12/2022
# Zhan RuiXin 710082148

import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate


#Define pre step size.
h=0.01
print("h=")
print(h)
#To make the error less than 1%, it was tested and found to be possible at 0.01 per step.

#Define the end time, at 50 seconds the system is in a stable state.
t = np.arange(0, 50, h) 

#Define input force is a unit inpluse
def impulse(t):
    if t < 1:
        return 1
    else:
        return 0

#Define the matrix Z which is calculated from equation 6 and equation 7 (Can find in pdf) 
def Matrix(t,Z):
    Z1 = Z[1]
    Z2 = Z[3] - Z[1] - Z[0]
    Z3 = Z[3]
    Z4 = impulse(t)/2 + Z[1] - Z[3]
    return np.array([Z1, Z2, Z3, Z4])
 # t is one unit of time
 

#Define Euler method
def Euler(t):
    Z = np.zeros((len(t), 4))
    '''
    Define the number of rows of the matrix 
    as the number of time instances.
    The same applies to the Heun method matrix
    '''
#The iterative process of Euler's method
    i=0
    while(i<len(t)-1):
     Z[i]= Z[i]+ h * Matrix(t[i], Z[i])
     Z[i+1]=Z[i]
     i=i+1
    return Z[:,2]
#Output the Euler method calculation results
Result_Euler=Euler(t)


#Define Heun method    
def Heun(t):
    Z = np.zeros((len(t), 4))
#Defining the constant used in Heun method.(From slide in week 8)
    a_1 = 1/2
    a_2 = 1/2
    p_1 = 1
    q_11 = 1
#The iterative process of Heun method    
    i=0
    while(i<len(t)-1):
     k_1 = Matrix(t[i], Z[i]) 
     k_2 = Matrix(t[i] + p_1 * h, Z[i] + q_11 * k_1 * h)
     Z[i] = Z[i] + (a_1 * k_1 + a_2 * k_2) * h
     Z[i+1]=Z[i]
     i=i+1
    return Z[:,2]
#Output the Heun method calculation results
Result_Heun=Heun(t)


#Define Odeint method
F = lambda Z,t: [Z[1],Z[3] - Z[1] - Z[0],Z[3],impulse(t)/2 + Z[1] - Z[3]]
Odeint= integrate.odeint(F,[0,0,0,0],t)  
#Output calculation results 
Result_Odeinet=Odeint[:,2]


#Plot the figure of Euler method
plt.plot(t,Result_Euler,label='Euler',color='g')
plt.legend()
plt.grid()
plt.title('Result of Euler method')
plt.show()
plt.savefig('./Ruler.jpg')

#Plot the figure of Heun method
plt.plot(t,Result_Heun,label='Heun',color='b')
plt.legend()
plt.grid()
plt.title('Result of Heun method')
plt.show()

#Plot the figure of Ode method
plt.plot(t,Result_Odeinet,label='Odeint',color='r')
plt.legend()
plt.grid()
plt.title('Result of Odeint method')
plt.show()

#Plot the figure of three method together
plt.plot(t,Result_Euler,label='Euler',color='g')
plt.plot(t,Result_Heun,label='Heun',color='b')
plt.plot(t,Result_Odeinet,label='Odeint',color='r')
plt.legend()
plt.grid()
plt.title('Result together')
plt.show()

#Error analysis
'''
   The exact solution calculated by odeint is compared
   with the predictions calculated by Euler's and Heun's 
   methods to find the percentage of error.
   The result show in the figure.
'''

Error_euler=np.zeros(len(t))  
Error_Heun=np.zeros(len(t))

 # Error calculation
def E_t(y,s):
   result=abs((y-s)/s*100)
   return result

for i in range(len(t)-1): 
 #Euler error
 Error_euler[i]=E_t(Result_Euler[i],Result_Odeinet[i])
 #Heun error
 Error_Heun[i]=E_t(Result_Heun[i],Result_Odeinet[i])

#Plot the figure of the result of erroy calculation
plt.plot(t,Error_euler,label='Euler')
plt.plot(t,Error_Heun,label='Heun')
plt.ylim((0, 5))
plt.yticks(np.arange(0, 5, 1))
plt.legend()
plt.title('Error in this step size')
plt.xlabel('$t$')
plt.ylabel('$Error\%$')
plt.grid()
plt.show()
plt.savefig('./error.jpg')


 # Error value calculation
E_Euler = abs(1/2-abs(Result_Euler[-1]))
E_Euler_percentage=E_Euler/0.5*100
print("The Euler method error is",E_Euler)
print("The Euler method error percentage is",E_Euler_percentage)


E_Heun = abs(1/2-abs(Result_Heun[-1]))
E_Heun_percentage=E_Heun/0.5*100
print("The Heun method error is",E_Heun)
print("The Heun method error percentage is",E_Heun_percentage)


