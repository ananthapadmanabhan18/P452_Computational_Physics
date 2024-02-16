############################################################################################################################
#LCG with parameters already included
def lcg1(p,n):
    seed=p
    a=1103515245
    c=12345
    m=32768
    xlist=[]
    ylist=[]
    p=seed
    for k in range(1,n+1):
        xlist.append(k)
        (q)=((a*p)+c)%(m)
        ylist.append(q)
        p=q
    return(ylist)




#LCG without parameters already included-u have to give the parameters
def lcg2(i,a,c,m,n):
    xlist=[]
    ylist=[]
    p=i
    for k in range(1,500):
        xlist.append(k)
        (q)=((a*p)+c)%(m)
        ylist.append(q)
        p=q
    return(ylist)









class rng():
    def __init__(self,seed, a = 1103515245, c = 12345 ,m = 32768):
        # initiation of data input
        self.term = seed
        self.a = a
        self.c = c
        self.m = m
    def gen(self):
        # generates a random number
        self.term = (((self.a * self.term) + self.c) % self.m)
        return self.term
    def genlist(self,length):
        # returns a list of 'n' random numbers in the range (0,1) where 'n' is 'length'.
        RNs = []
        for i in range(length):
            self.term = (((self.a * self.term) + self.c) % self.m)
            RNs.append(self.term / self.m)
        return RNs    
##########################################################################################################################











#To print the matrix
def printmatrix(A):
    for i in range(0,len(A)):
        print("\n")
        for j in range(0,len(A[i])):
            print(A[i][j],",",end='')
    print("\n")    









#checking if the matrix is symmetric or not
def symcheck(A):
    for i in range(0,len(A)):
        for j in range(0,len(A[i])):  
            if A[i][j]!=A[j][i]:
                return False
            else:
                continue    
    if not False:
        return True
    else:
        return False         


def inv_solve_LU(A):
    n=len(A)
    Is = [0 for i in range(n)]
    for j in range(n):
         Is[j] = [[0] for k in range(n)]
         Is[j][j][0] = 1
    invA = [[] for l in range(n)]
    for i in range(n):
        tA = copy_matrix(A)
        inv = eqn_solve_lu(tA,Is[i]) 
        for j in range(n):
            invA[j].append(inv[j][0])
    return invA
def get_det(A):
    n = len(A)
    if n != len(A[0]):
        print('Not a square matrix')
        return None
    for curr in range(n):
        if A[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,n):

                if abs(A[row][curr]) > abs(A[max_row][curr]):
                    max_row = row
            if max_row == curr:
                print('The matrix is singular!')
                return None
            A[curr],A[max_row] = A[max_row], A[curr]
        for i in range(curr + 1,n):
            if A[i][curr] != 0:
                lead_term = A[i][curr]/A[curr][curr]
                for j in range(curr,len(A[i])):
                    A[i][j] = A[i][j] - (A[curr][j] * lead_term)
    prdct = 1
    for i in range(n):
        prdct *= A[i][i]
    return prdct

def mat_mul(A,B):
    result=[[0 for i in range(len(B[0]))] for i in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result            


def copy_matrix(A):
    newA = []
    for i in range(len(A)):
        newA.append(A[i][:])
    return newA



def dot_prod(A,B):
    U=get_transpose(A)
    P=mat_mul(U,B)
    return P


def get_transpose(matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    matrix_T = []
    for j in range(columns):
        row = []
        for i in range(rows):
           row.append(matrix[i][j])
        matrix_T.append(row)
    return matrix_T


def get_inv_mat_GJ(A):

    if len(A) != len(A[0]): #if not square matrix, exit
        raise ValueError('Matrix is not square')

    n = len(A) #the dimension of matrix will be n*n
    I = []
    for row in range(n):
        I.append(list())
        for col in range(n):
            if col == row:
                I[row].append(1)
            else:
                I[row].append(0)
    
    
    for curr in range(n): #curr takes the index of each column we are processing
        #row swap if zero
        if A[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,n):

                if abs(A[row][curr]) > abs(A[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                return None
            A[curr],A[max_row] = A[max_row], A[curr]
            I[curr],I[max_row] = I[max_row], I[curr]
        #making the pivot element 1
        if A[curr][curr] != 1:
            pivot = A[curr][curr]
            for i in range(n):
                A[curr][i] = A[curr][i] / pivot
                I[curr][i] = I[curr][i] / pivot

        #making others zero
        for i in range(0,n):
            if i == curr: #skipping the pivot point
                continue
            if A[i][curr] != 0:
                lead = A[i][curr]
                for j in range(0,n):
                    A[i][j] = A[i][j] - (A[curr][j] * lead)
                    I[i][j] = I[i][j] - (I[curr][j] * lead)

    return I

def solve_GJ(A,B):
    #solves linear equations using Gauss-Jordon method
    #takes the A and B matrix as input from the form A.X = B where X is the unknown matrix
    #returns solved X 

    #constructing augmented matrix 
    augmat = A
    for row in range(len(augmat)):
        augmat[row].append(B[row][0])
    
    for curr in range(len(augmat)): #curr takes the index of each column we are processing
        #row swap if zero
        if augmat[curr][curr] == 0:
            max_row = curr
            for row in range(curr + 1,len(augmat)):
                if abs(augmat[row][curr]) > abs(augmat[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                return None
            augmat[curr],augmat[max_row] = augmat[max_row], augmat[curr]

        #making the pivot element 1
        if augmat[curr][curr] != 1:
            pivot_term = augmat[curr][curr]
            for i in range(len(augmat[curr])):
                augmat[curr][i] = augmat[curr][i] / pivot_term

        #making others zero
        for i in range(0,len(augmat)):
            if i == curr: #skipping the pivot point
                continue
            if augmat[i][curr] != 0:
                lead_term = augmat[i][curr]
                for j in range(curr,len(augmat[i])): #elements before the curr column are zero in curr row, so no need to calculate
                    augmat[i][j] = augmat[i][j] - (augmat[curr][j] * lead_term)
    solution = []
    for i in range(len(augmat)):
        solution.append([augmat[i][-1]]) #Taking last elements into a list to form column matrix
    return solution

    


#Gauss Jordan


def getone(sd,pp):
    for i in range(len(sd[0])):
        if sd[pp][pp] != 1:
            q00 = sd[pp][pp]

            for j in range(len(sd[0])):
                sd[pp][j] = sd[pp][j] / q00

def getzero(sd,r, c):
    for i in range(len(sd[0])):
        if sd[r][c] != 0:
            q04 = sd[r][c]

            for j in range(len(sd[0])):
                sd[r][j] = sd[r][j] - ((q04) * sd[c][j])



def gauss_jordan(sd):
    for i in range(len(sd)):
        getone(sd,i)

        for j in range(len(sd)):
            if i != j:
                getzero(sd,j, i)
    return sd            


def get_inv(A):
    A=gauss_jordan(A)
    inv=[[0 for i in range(int(len(A[1])/2))] for i in range(int(len(A[1])/2))]
    for i in range(len(inv)):
        for j in range(len(inv)):
            inv[i][j]=A[i][int((len(A)/2))+j+1]
    return inv   






#Cholesky_Decomposition
def Chol_dec(A):  
    from math import sqrt
    if len(A) != len(A[0]):
        return None
    n = len(A)

    for row in range(n):
        for col in range(row,n): 
            if row == col:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] ** 2)
                A[row][row] = sqrt(A[row][row] - sum)
            else:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] * A[i][col])
                A[row][col] = (A[row][col] - sum) / A[row][row]
                A[col][row] = A[row][col]   
    for row in range(n):
        for col in range(row + 1,n):
            A[row][col] = 0   
    return A

def sup_chol_dec(A):
    L=Chol_dec(A)
    LT=get_transpose(L)
    LP=[[0 for i in range(len(L))] for i in range(len(L))]
    for i in range(len(L)):
        for j in range(len(L)):
            LP[i][j]=L[i][j]+LT[i][j]
    for i in range(len(L)):
        LP[i][i]=LP[i][i]/2
    return LP    

def for_back_chol_dec(A,B):
    Y = []
    n = len(B)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += (A[i][j] * Y[j])
        Y.append((B[i][0]-sum)/A[i][i])

    X = Y
    for i in range(n-1,-1,-1):
        sum = 0
        for j in range(i + 1, n):
            sum += (A[i][j] * X[j])
        X[i] = (Y[i] - sum) / A[i][i]

    for i in range(n):
        X[i] = [X[i]]

    return Y


def lu_dec(A):
    

    n = len(A) 
    if n != len(A[0]): 
        print('Not square!')
        return None

    for j in range(n):
        for i in range(1,n):
            if i <= j:
                sum = 0
                for k in range(0,i):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = A[i][j] - sum
            else:
                sum = 0
                for k in range(0,j):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = (A[i][j] - sum) / A[j][j]
    return A






def for_back_LU(A,B):
    Y = []
    n = len(B)
    for i in range(n):
        s = 0
        for j in range(i):
            s += (A[i][j] * Y[j])
        Y.append(B[i][0]-s)
    
    X = Y[:]
    for i in range(n-1,-1,-1):
        s = 0
        for j in range(i + 1, n):
            s+= (A[i][j] * X[j])
        X[i] = (Y[i] - s) / A[i][i]

    for i in range(n):
        X[i] = [X[i]]

    return X

#Iteration methods-Solving linear equation-Gauss seidel method
def diagonal_dom_check(m) :
    n=len(m)
    for i in range(0, n) :        
        sum = 0
        for j in range(0, n) :
            sum = sum + abs(m[i][j])    
        sum = sum - abs(m[i][i])
        if (abs(m[i][i]) < sum) :
            return False 
    return True

def make_DD(A,B):
    n = len(A)
    for i in range(n):
        sum=0
        for j in range(n):
            if j != i:
                sum += abs(A[i][j])
        if abs(A[i][i])>sum:
            continue
        else:
            c = i + 1
            flag = 0
            while c<n:
                sum = 0
                for j in range(n):
                    if j != i:
                        sum += abs(A[c][j])
                if abs(A[c][i])>sum:
                    flag = 1
                    break
                else:
                    c+=1
            if flag==0:
                return None,None
            else:
                A[i],A[c]=A[c],A[i]
                B[i],B[c]=B[c],B[i]
    return A,B


def eqn_solve_lu(A,B):
    L=lu_dec(A)
    L1=for_back_LU(L,B)
    return L1




def gauss_seidel(A,B,tolerance):
    n = len(A)
    X = list([0] for i in range(n))
    for steps in range(200):
        flag = 1
        for i in range(n):
            sum = 0
            for j in range(i):
               sum += (A[i][j] * X[j][0])
            for j in range(i+1,n):
                sum += (A[i][j] * X[j][0])
            temp = (B[i][0] - sum) / (A[i][i])
            if (abs(temp) - abs(X[i][0])) > tolerance:
               flag = 0
            X[i][0] = temp
        if flag == 1:
            return X,steps + 1
    print('Eqn not solved after 200 steps')
    return None,200

def GSM_with_DD_check(A,B,T):
    if diagonal_dom_check(A)==True:
        return gauss_seidel(A,B,T)
    else:
        print(" The given Matrix is not Diagonallly Dominant")    
    





#Iteration methods-Solving linear equation-Gauss-Jacobi method
def jacobi_method(A,B,guess,tol):
    XK=guess
    XK1=[[0] for i in range(len(A))]
    count=0
    sumaijxj=0
    flag=0
    while flag!=1:
        for i in range(len(XK1)):
            sumaijxj=0
            for j in range(len(A[i])):
                if i!=j:
                    sumaijxj+=A[i][j]*XK[j][0]    
            XK1[i][0]=(1/A[i][i])*(B[i][0]-sumaijxj)
        for i in range(len(A)):
            if abs(XK[i][0]-XK1[i][0])<tol:
                flag=1    
        count+=1
        for i in range(len(XK1)):
            XK[i][0]=XK1[i][0]
    return XK,count

def GJacobi_with_DD_Check(A,B,guess,T):
    if diagonal_dom_check(A)==True:
        return jacobi_method(A,B,guess,T)
    else:
        print(" The given Matrix was not Diagonallly Dominant") 
        A,B=make_DD(A,B)
        return jacobi_method(A,B,guess,T)



#Roots of non linear equation-BIsection method
def bracket(a0,b0,f):
    n=0
    while f(a0)*f(b0)>=0:
        if abs(f(a0))>abs(f(b0)):
            b0=b0+1.5*(b0-a0)
        else:
            a0=a0-1.5*(b0-a0)       
    return(a0,b0)

def bijection(a0,b0,f,T):
    a0,b0=bracket(a0,b0,f)
    epsilon=T
    delta=0.001
    count=0
    while (abs(b0-a0))>epsilon:
        c0=(a0+b0)/2
        if f(a0)*f(c0)>0:
            a0=c0
        else:
            b0 = c0 
        count+=1       
    return c0,count 




#Roots of non linear equation-Newton raphson
def newton_raphson(x0,f,fd,T):
    count=0
    xn=x0
    xn1=xn-(f(xn)/fd(xn))
    while True:
        xn=xn1
        xn1=xn-(f(xn)/fd(xn)) 
        count+=1
        if abs(xn1-xn)<T:
            break          
    return xn1,count+1  



#Roots of non linear equation-Regula falsi
def regula_falsi(a0,b0,f,T):
    epsilon=T
    a0,b0=bracket(a0,b0,f)
    for i in range(0,1):
        c0=b0-(((b0-a0)*f(b0))/(f(b0)-f(a0)))
        cold=c0  
        if f(a0)*f(c0)<0:
            b0=c0       
        else:    
            a0=c0 

     
    count=1

    while True:
        cold=c0
        c0=b0-(((b0-a0)*f(b0))/(f(b0)-f(a0)))
        if f(a0)*f(c0)<0:
            b0=c0       
        else:    
           a0=c0 
        if abs(cold-c0)<epsilon:
            break 
        count+=1
    return c0,count



import math as m

def px_deflate(px, root):

    #Here, the passed 'root' is a verified root from main body with a certain tolerance, therefore checking not done.
    if len(px) == 1:
        print('P(x) doesnt contain any x: deflation exited!')
        return None
    if px[0] != 1:
        lead = px[0]
        for i in range(len(px)):
            px[i] = px[i] / lead
    for i in range(1,len(px)):
        px[i] = (px[i-1] * root) + px[i]
    return px[:-1]



def px_derivative(px):
    n = len(px) - 1
    i = 0
    while n != 0:
        px[i] = px[i] * (n)
        i += 1
        n -= 1
    return px[:-1]

def px_value(px,x):
    print('finding',px)
    sum = 0
    i = 0
    n = len(px) - 1
    while n >= 0:
        sum+= (px[i] * (x**n))
        i += 1
        n -= 1
    return sum




def root_laguire(px,guess,tolerance):
    roots = []
    steps = 1
    N = len(px) - 1
    n = N
    while steps <= N:
        if abs(px_value(px,guess)) < tolerance:
            print('{}th root {} found in 0 steps.'.format(steps,guess))
            steps += 1
            roots.append(guess)
            px = px_deflate(px,guess)
            n -= 1
            continue
        # print('Finding the {}th root for {}'.format(steps,px))
        dpx = px_derivative(px[:])
        ddpx = px_derivative(dpx[:])
        i = 1
        theguess = guess
        while True:
            g = px_value(dpx,theguess)/px_value(px,theguess)
            h = (g**2) - (px_value(ddpx,theguess)/px_value(px,theguess))
            if g < 0:
                a = (n / (g - m.sqrt((n-1)*((n*h)-(g*2)))))
            else:
                a = (n / (g + m.sqrt((n-1)*((n*h)-(g*2)))))
            # print('a is',a)
            newguess = theguess - a
            # print(px_value(px,theguess),'and',px_value(px, newguess))
            
            if i < 26 :
                # if abs(a) < tolerance and px_value(px,newguess) < tolerance:
                if px_value(px,newguess) < tolerance:
                    print('{}th root is {} found in {} steps.'.format(steps,newguess,i))
                    steps += 1
                    roots.append(newguess)
                    px = px_deflate(px,newguess)
                    n -= 1
                    break
                # else:
                #     print('Guess discarded.\n')

            else:
                print('The guess for {} th root was not found in 25 steps'.format(steps))
                return None
            theguess = newguess
            i += 1
    return roots
            


# interpolation
def l_interpolation(xlist,ylist,x):
    s=0
    for i in range(0,len(xlist)):
        p=1
        for k in range(0,len(xlist)):
            if i==k:
                p=p
            else:    
               p=p*(((x-xlist[k])/(xlist[i]-xlist[k])))
        p=p*ylist[i]       
        s=s+p
    return s


#Least square fitting
def line_fit(X,Y,degree):
    n=len(X)
    k=degree
    A=[[0 for i in range(k+1)] for i in range(k+1)]
    B=[[0] for i in range(k+1)]

    for i in range(0,k+1):
        for j in range(0,k+1):
            sum=0
            for t in range(n):
                sum+=pow(X[t],i+j)
            A[i][j]=sum
    A[0][0]=n
    for i in range(0,k+1):
            sum=0
            for t in range(n):
                sum+=pow(X[t],i)*Y[t]
            B[i][0]=sum    
    return A,B


def line_fit_coeff(x,y,degree):
    A,B=line_fit(x,y,degree)
    coeff=eqn_solve_lu(A,B)
    return coeff


def linear_fitting(x,y):
    import math
    Sx, Sy,Sxx,Syy,Sxy,S=0,0,0,0,0,len(x)
    for i in range(len(x)):
        Sx=Sx+x[i]
        Sy=Sy+y[i]
        Sxx+=(x[i])**2
        Syy+=(y[i])**2
        Sxy+=(x[i])*(y[i])
    d=S*Sxx-(Sx)**2
    w0=(Sxx*Sy-Sx*Sxy)/d
    wc=(Sxy*S-Sx*Sy)/d
    r=Sxy/(math.sqrt((Sxx*2) * (Syy*2)))
    return w0,wc






def pearson_r(x,y):
    #for sigma =1
    n=len(x)
    (sx,sy,sxy,sx2,sy2)=(0,0,0,0,0)
    for i in range(n):
        sx+=x[i]
        sy+=y[i]
        sxy+=x[i]*y[i]
        sx2+=x[i]**2
        sy2+=y[i]**2
    R=(n*sxy - sx*sy)/(((n*sx2 - sx**2)**0.5)*((n*sy2 - sy**2)**0.5))
    return abs(R)    




# file = open('Assignments\Assignment-3\Q2-Input.txt')
# lines = file.readlines()
# lineindex = 1
# while lineindex != len(lines):
#     A = []                                                                            #reading from file 
#     while lines[lineindex][0] != '#':
#         linelist = lines[lineindex].split(',')
#         for i in range(len(linelist)):
#             linelist[i]= float(linelist[i])
#         lineindex += 1
#         A.append(linelist)
#     lineindex += 1
#     B = []
#     while lines[lineindex][0] != '#':
#         linelist = lines[lineindex].split(',')
#         for i in range(len(linelist)):
#             linelist[i]= float(linelist[i])
#         lineindex += 1
#         B.append(linelist)
#     lineindex += 1    



#integration##########################
def midpoint_int(f,a,b,N):
    n=N
    h=(b-a)/N
    xi=0
    sum=0
    for i in range(n):
        xi=((a+(i*h))+(a+((i+1)*h)))/2
        sum+=h*f(xi)
    return sum

def trap_int(f,a,b,N):
    n=N
    h=(b-a)/N
    xi=0
    sum=0
    xi=a
    for i in range(1,n):
        xi1=a+i*h
        sum+=0.5*h*(f(xi)+f(xi1))
        xi=xi1
    return sum

def simpson_int(f,a,b,N):
    def w(n,nmax):
        if n==nmax or n==0:
            return 1
        else:    
            if n%2==0:
                return 2
            else:
                return 4
    h=(b-a)/(N)
    x0=a
    sn=0
    for i in range(N+1):
        sn+=(h/3)*w(i,N)*(f(x0))
        x0=x0+h
    return sn 



def monte_carlo(a,b,f,N,seed):
    randno=lcg1(seed,N)
    F=0
    for i in range(len(randno)):
        randno[i]=((b-a)*(randno[i]/32768))+a
        F+=((b-a)*f(randno[i]))/N   
    return F   


def monte_carlo_error(a,b,f,N,seed):
    randno=lcg1(seed,N)
    F=0
    F1=0
    for i in range(len(randno)):
        randno[i]=((b-a)*(randno[i]/32768))+a
        F+=f(randno[i])
        F1+=pow(f(randno[i]),2)  
    return (F1/N)-pow((F/N),2)    


def RK4_1eqn(x0,y0,h,f,limit):
    xlist=[]
    ylist=[]
    while x0<limit:
        xlist.append(x0)
        k1=h*(f(y0,x0))
        k2=h*f(y0+(k1/2),x0+h/2)
        k3=h*f(y0+(k2/2),x0+h/2)
        k4=h*(f(y0+k3,x0+h))
        ylist.append(y0)
        y0=y0+((1/6)*(k1+2*k2+2*k3+k4))
        x0=x0+h
    return xlist,ylist



def gen_RK4_coupled(fnlist,x0,y0s,limit,h):
    limit -= h/2 
    n = len(y0s) 
    k1 = [0 for i in range(n)]
    k2 = [0 for i in range(n)]
    k3 = [0 for i in range(n)]
    k4 = [0 for i in range(n)]
    tys= [0 for i in range(n)] 
    datY = []
    for i in range(n):
        datY.append([y0s[i]])
    datT = [x0]
    while x0 < limit:
        for i in range(n):
            k1[i] = h*fnlist[i](y0s,x0)    
        for i in range(n):
            tys[i] = y0s[i] + (k1[i] / 2)
        for i in range(n):
            k2[i] = h*fnlist[i](tys, (x0 + (h/2)))
        for i in range(n):
            tys[i] = y0s[i] + (k2[i] / 2)
        for i in range(n):
            k3[i] = h*fnlist[i](tys, (x0 + (h/2)))   
        for i in range(n):
            tys[i] = y0s[i] + k3[i]
        for i in range(n):
            k4[i] = h*fnlist[i](tys, (x0 + h))
        for i in range(n):
            y0s[i] += ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + (k4[i])) / 6)
        x0 += h
        for i in range(n):
            datY[i].append(y0s[i])
        datT.append(x0)
    return datT, datY


def shooting_method(fns,x0,y0,x1,y1,guess1,tol,h):    
    X,Y = gen_RK4_coupled(fns,x0,[y0,guess1],x1,h)
    ye1 = Y[0][-1]
    if abs(ye1 - y1) < tol:
        return X,Y
    if ye1 < y1:
        guess1side = -1
    else :
        guess1side = 1
    guess2 = guess1 + 2   
    X,Y =gen_RK4_coupled(fns,x0,[y0,guess2],x1,h)
    ye2 = Y[0][-1]
    if ye2 < y1:
        guess2side = -1
    else :
        guess2side = 1
    while guess1side * guess2side != -1:

        if abs(y1-ye2) > abs(y1-ye1):
            guess2 = guess1 - abs(guess2-guess1)
        else:
            guess2 += abs(guess2-guess1)
        X,Y = gen_RK4_coupled(fns,x0,[y0,guess2],x1,h)
        ye2 = Y[0][-1]
        if ye2 < y1:
            guess2side = -1
        else :
            guess2side = 1
    i = 0
    while True:
        newguess = guess1 + (((guess2 - guess1)/(ye1 - ye2))*(y1 - ye2))
        i += 1
        X,Y = gen_RK4_coupled(fns,x0,[y0,newguess],x1,h)
        yvalnew = Y[0][-1]
        if abs(yvalnew - y1) <tol:
            break
        if yvalnew < y1:
            guess1 = newguess
            ye1 = yvalnew
        else:
            guess2 = newguess
            ye2 = yvalnew
    return X,Y


def pde_explicit(f,nx,nt,lx,lt,N):
    hx=lx/nx
    ht=lt/nt
    a = ht/(pow(hx,2))
    V0,V1 = [0],[0]
    for i in range(nx+1):
        V1.append(0)
    for i in range(nx+1):
        V0.append(f(i*hx))
    for j in range(N):
        for i in range(nx+1):
            if i==0:
                V1[i]=(1-2*a)*V0[i] + (a*V0[i+1])
            elif i==nx:
                V1[i]=(a*V0[i-1]) + (1-2*a)*V0[i]
            else:
                V1[i]=(a*V0[i-1]) + (1-2*a)*V0[i] + (a*V0[i+1])
        for i in range(nx+1):
            V0[i] = V1[i]        
    if N==0:
        return V0
    else:
        return V1









def eigen_val_power_iter(A,x0,tol):
    ev0=mat_mul(A,x0)
    den=(dot_prod(mat_mul(A,x0),x0))[0][0]
    ev1=mat_mul(A,ev0)
    num=dot_prod(ev1,x0)[0][0]
    l=num/den
    den=num
    ev0=ev1
    ev1=mat_mul(A,ev0)
    num=dot_prod(ev1,x0)[0][0]
    l1=num/den
    den=num
    ev0=ev1    
    i=2
    while abs(l1-l)>tol: 
        l=l1 
        ev1=mat_mul(A,ev0)
        num=dot_prod(ev1,x0)[0][0]
        l1=num/den
        den=num
        i+=1
        ev0=ev1
    sum=0
    for k in range(len(ev1)):
        sum+=ev1[k][0]**2
    for j in range(len(ev1)):
        ev1[j][0]=ev1[j][0]/(sum**(0.5))
    return l1,ev1,i


def radioactivity(Nlist,tlist,step,tf,seed1):
    n = len(Nlist)
    Nlist.append(0)
    print(n)
    result = Nlist[:]
    for i in range(len(result)):
        result[i] = [result[i]]
    for i in range(n):
        tlist[i] = ((m.log(2)/tlist[i]) * step)
    print(tlist)
    ran = rng(seed = seed1)
    t = step
    times = [0]
    while t <= tf:
        for i in range(n):
            trials = Nlist[i]
            for j in range(trials):
                if ran.gen() <= tlist[i]:
                    # print("decayed")
                    Nlist[i] -= 1
                    Nlist[i+1] += 1
        times.append(t)
        for i in range(len(result)):
            result[i].append(Nlist[i])
        t += step
    return times,result    

































def partition_problem(N,dt,limit,seed):
    N0 = N
    dat=[N,0]
    Ns = dat[:]
    for i in range(len(dat)):
        dat[i] = [dat[i]]
    ran = rng(seed)
    t = dt
    ts = [0]
    while t <=limit:
        cr = Ns[0]
        p = cr/N0
        if ran.gen() <p:
            Ns[0] -= 1
            Ns[1] += 1
        else:
            Ns[0] += 1
            Ns[1] -= 1
        ts.append(t)
        for i in range(len(dat)):
            dat[i].append(Ns[i])
        t += dt
    return ts,dat




