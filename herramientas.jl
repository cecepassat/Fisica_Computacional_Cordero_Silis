
## Comandos para definir la biblioteca

__precompile__()
module herramientas

## Programas

##### Método de Newton
export Metodo_Newton
function Metodo_Newton(f,derivadaf,inicial)
    
    x=inicial;
    for i in 1:100
        x=x-f(x)/derivadaf(x)
    end
    #La función regresa el último valor de x obtenido
    return x
end

##### Newton Funciones Arbitrarias
export programa
function metodo(f,derivada,inicial)
      #x=Sym("x")
   
    
    list=zeros(Complex64,30);
    j=inicial;
    for i in 1:30
            #@show [j]
        j=j-f(j)/derivada(j);
        list[i]=j;
    end
    
    return j
end

##### Derivada Numérica
export DerivadaNumerica
function DerivadaNumerica(f,x0,h)
    return df=(f(x0+h)-f(x0))/h
end


##### Derivada Simétrica
export Derivada_Simetrica
function Derivada_Simetrica(f,x0,h)
    (f(x0+h)-f(x0-h))/2h
end

##### Integración método del Rectángulo
export Rectangulo
function Rectangulo(f, a, b, n)
    S=0
    A = linspace.(a,b,n)      ##Se define la partición del intervalo X
    for j in 1:n-1            ##Iteramos la suma usando el método del rectángulo
        S = S + (A[j+1]-A[j])*f((A[j+1]+A[j])/2)
    end   
    S    
end 

##### Integración método del Trapecio
export Trapecio
function Trapecio(f, a, b, n)
    S=0
    A = linspace.(a,b,n)
    for j in 1:n-1
        S = S + (A[j+1]-A[j])*(f(A[j+1])+f(A[j]))/2     ##Iteramos la suma con el método del trapecio
    end   
    S    
end

##### Integración Simpson
export Simpson
function Simpson(f, a, b, n)
    S=0
    A = linspace.(a,b,n)
    for j in 1:n-1
        p = (A[j+1] + A[j])/2
        S = S + (A[j+1]-A[j])/6 * (f(A[j]) + 4f(p) + f(A[j+1]))     ##Iteramos la suma con el método de simpson
    end   
    S    
end

##### Euler Simétrico
export Metodo_Euler
function Metodo_Euler(f,x0,t0,tf,h)
    x=x0
    t=t0
    ## Definimos las listas a utilizar "listt" y "listx", así como las funciones a las que están asociadas x,t
    listt=[]
    listx=[]
    n = round((tf-t0)/h) ## Aplicamos la fórmula de Euler
    for i in 1:n-1
        x += h*f(x,t)
        t += h
        push!(listt,t) ## Limitamos los valores que se tomarán para no saturarse de resultados
        push!(listx,x)
    end
    return listt, listx
end


##### Euler más de una dimensión
export Metodo_Euler2
function Metodo_Euler2(f,x0,u0,listt)
    n = length(listt)   
    listx = zeros(n)            
    listu = zeros(n)             
    listx[1] = x0                
    
    listu[1] = u0                
    h = (listt[n]-listt[1])/n   
    for i in 1:n-1
        listx[i+1] = listx[i] + h*listu[i]                       
        listu[i+1] = listu[i] + h*f(listx[i],listu[i],listt[i])  
    end
    return listx
end



##### Usando ahora el método de Newton definimos el método implícito de Euler (en función del de Newton) agregando ahora una dependencia de una lista (listt)
export Imp_Euler
function Imp_Euler(f,derivadaf,listt,x0)
## Definimos listx como un arreglo de raíces del mismo orden (entiéndase cardinalidad) que la lista de t's.
    listx=zeros(length(listt))
    listx[1]=x0
    h=listt[2]-listt[1]     ## Es la distancia entre las primeras entradas
    for i in 1:length(listt)-1     ## Guardamos las entradas de listx siguientes
        
## La siguiente función es a la que tenemos que encontrar su raíz (usando el método implícito) es decir, el "elemento" x_{i+1}.
    g(z)=z-f(z,listt[i+1])*h-listx[i]
    derivadag(z)= 1-derivadaf(z,listt[i+1])*h     ## Consideremos también su derivada
        
        
## Aplicamos el método de Newton a g, usando la condición inicial x_{i} y tomando cada raíz como el siguiente elemento(x_{i+1}).
    listx[i+1]=Metodo_Newton(g,derivadag,listx[i])    
    end
    
    return listx     ## Obtenemos los valores para x
    
end

##### Euler Punto Medio

## Definimos la función para el método del punto medio
export Punto_M
function Punto_M(f,listt,inicial)
    h=listt[2]-listt[1]     ## Nuevamente tenemos la distancia entre las primeros dos entradas
    
    listx=zeros(length(listt))     ## Misma longitud de las listas
    
    listx[1]=inicial
    
## Apliquemos las iteraciones dada la condición de recurrencia
    for i in 1:length(listt)-1
       listx[i+1]= listx[i]+h*f(listx[i]+(h/2)*f(listx[i],listt[i]),listt[i]+h/2)
    end 
    
    return listx     ## Obtenemos nuestros resultados
    
end

##### Runge-Kutta orden 4
export RK_04
function RK_O4(f,listt,x0)
    listx=zeros(length(listt))
    listx[1]=x0
    h=listt[2]-listt[1]
    
    for i in 1:length(listt)-1
        it1=f(listx[i],listt[i])
        it2=f(listx[i]+(h/2)*it1,listt[i]+h/2)
        it3=f(listx[i]+(h/2)*it2,listt[i]+h/2)
        it4=f(listx[i]+h*it3,listt[i+1])     ## Aplicamos la relación de recurrencia para cada iteración hasta orden 4
        
        listx[i+1]=listx[i]+(h/6)*(it1+2*it2+2*it3+it4)     ## Lista donde se guardan las x's obtenidas
        
    end
      
    return listx     ## Obtenemos las x's del paso anterior
    
end

end
