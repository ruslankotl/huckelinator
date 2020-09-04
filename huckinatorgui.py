# -*- coding: utf-8 -*-
"""

Fit GUI to the Huckinator
Created on Mon Aug 31 11:10:37 2020

@author: Руслан Алексеевич
"""

import tkinter as tk
import tkinter.simpledialog
import tkinter.ttk

import numpy as np

import math



class Application(tk.Frame):
    def __init__(self, master = None):
        #initialise the program
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    '''Interface functions'''

    def create_widgets(self):
        #main menu
        self.lin = tk.Button(self)
        self.lin['text'] = 'Linear polyenes'
        self.lin['command'] = self.linear
        self.lin.grid(row=0, column=0)

        self.cyc=tk.Button(self, text='Cyclic polyenes',
                         command=self.cyclic)
        self.cyc.grid(row=1, column=0)

        self.poly=tk.Button(self, text='Selected polyhedra',
                         command=self.polyhedron)
        self.poly.grid(row=2, column=0)

        self.quit=tk.Button(self, text = 'QUIT', fg = 'red',
                         command=self.master.destroy)
        self.quit.grid(row=3, column = 0)


    def linear(self):

        length = tkinter.simpledialog.askinteger("Linear polyene", "Number of atoms",
                                 parent=self.master,
                                 minvalue=1)


        try:
            lm = self.linear_matrix(length)
            self.result(lm)
        except TypeError:
            pass







    def cyclic(self):

        length = tkinter.simpledialog.askinteger("Cyclic polyene", "Number of atoms",
                                 parent=self.master,
                                 minvalue=1)

        try:
            cm = self.cyclic_matrix(length)
            self.result_deg(cm)
        except TypeError:
            pass

    def polyhedron(self):

        newWindow=tk.Toplevel(self.master)
        label=tk.Label(newWindow, text='Choose the polyhedron')
        label.grid(row=0, column=0, columnspan=3)

        tetr = tk.Button(newWindow, text='Tetrahedron', command=lambda: self.result_deg(self.tetrahedron()))
        cube = tk.Button(newWindow, text='Cube', command=lambda: self.result_deg(self.cube()))
        octa = tk.Button(newWindow, text='Octahedron', command=lambda: self.result_deg(self.octahedron()))
        dode = tk.Button(newWindow, text='Dodecahedron', command=lambda: self.result_deg(self.dodecahedron()))
        icos = tk.Button(newWindow, text='Icosahedron', command=lambda: self.result_deg(self.icosahedron()))
        full = tk.Button(newWindow, text='C60-fullerene', command=lambda: self.result_deg(self.fullerene()))

        tetr.grid(row=1, column=0)
        cube.grid(row=1, column=1)
        octa.grid(row=1, column=2)
        dode.grid(row=2, column=0)
        icos.grid(row=2, column=1)
        full.grid(row=2, column=2)


        back=tk.Button(newWindow, text = 'Back', fg = '#ff00ff',
                         command=newWindow.destroy)
        back.grid(row=3, column=0, columnspan=3)


    '''Working functions'''


    def linear_matrix(self, n):
        #specify Huckel matrix for a linear polyene of length n

        H = np.zeros((n,n)) #initialise n*n matrix
        for i in range(n-1):
            H[i,i+1] = 1
        for i in range(n-1):
            H[i+1,i] = 1

        return H

    def cyclic_matrix(self, n):

        H = np.zeros((n,n))
        for i in range(n):
            H[i,(i+1)%n] = 1
        for i in range(n):
            H[(i+1)%n,i] = 1
        return H

    def get_evals(self, m):

        evals, evecs = np.linalg.eig(m)

        #print('Eigenvalues check out: '+ str(check_evals(evals)) + '\n')

        result = sorted(evals.real)

        #for each in result:
            #print(round(each, 4))

        return result

    def degeneracy_count(self, evals):
        out_list=[]
        count=0
        for i in range(len(evals)):
            if round(evals[i],4)==round(evals[(i+1)%(len(evals))],4):
                count+=1
            else:
                out_list.append([round(evals[i],3),count+1])
                count=0
        #print (out_list)
        return (out_list)


    def result(self, matrix):
        evals = self.get_evals(matrix)

        resultWindow = tk.Toplevel(self.master)
        result_box=tk.Text(resultWindow, width=25, height=10)
        result_box.grid(row=0, column=0)
        for value in evals:
            out_line=('alpha {:+4.3f} beta'.format(value))
            result_box.insert(tk.END, out_line+'\n')

        back=tk.Button(resultWindow, text = 'Back', fg = '#ff00ff',
                         command=resultWindow.destroy)
        back.grid(row=1, column=0)



    def result_deg(self, matrix):
        evals = self.get_evals(matrix)
        evals_output = self.degeneracy_count(evals)

        resultWindow = tk.Toplevel(self.master)
        result_box=tk.Text(resultWindow)
        result_box.grid(row=0, column=0)
        for value in evals_output:
            out_line='a{:+4.3f} b \t Degeneracy:{:1}'.format(value[0],value[1])
            result_box.insert(tk.END, out_line+'\n')

        back=tk.Button(resultWindow, text = 'Back', fg = '#ff00ff',
                         command=resultWindow.destroy)
        back.grid(row=1, column=0)

    '''polyhedra huckel matrices'''
    def tetrahedron(self):
        T=[[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]
        return T

    def cube(self):
        C = self.cyclic_matrix(8) #define a cube as an 8-membered ring, add extra connections
        C[0,3]=C[3,0]=1
        C[1,6]=C[6,1]=1
        C[2,5]=C[5,2]=1
        C[4,7]=C[7,4]=1
        return C

    def octahedron(self):
        O = self.cyclic_matrix(6) #define octahedron as a 6-membered ring, add missed connections
        O[0,2]=O[0,3]=O[1,4]=O[1,5]=O[2,0]=O[3,0]=1
        O[4,1]=O[5,1]=O[3,5]=O[5,3]=O[2,4]=O[4,2]=1
        return O


    def dodecahedron(self): #generating the matrix for dodecahedron
        #specifying coordinates
        phi = (math.sqrt(5)+1)/2
        coordinates = []
        for dummy in range(8):
            coordinates.append(((-1)**(dummy//4),(-1)**(dummy//2),(-1)**(dummy)))
        for dummy in range(4):
            coordinates.append((0,phi*(-1)**(dummy//2),(1/phi)*(-1)**(dummy)))
        for dummy in range(4):
            coordinates.append(((1/phi)*(-1)**(dummy//2),0,phi*(-1)**(dummy)))
        for dummy in range(4):
            coordinates.append((phi*(-1)**(dummy),(1/phi)*(-1)**(dummy//2),0))

        #create the matrix
        H = np.zeros((20,20))
        for r1 in coordinates:
            for r2 in coordinates:
                distance2 = 0
                distance2 =((r1[0])-(r2[0]))**2+((r1[1])-(r2[1]))**2+((r1[2])-(r2[2]))**2
                if round(distance2,6)==round((2/phi)**2,6):
                    H[coordinates.index(r1),coordinates.index(r2)] = 1
        return H



    def icosahedron(self):
        phi = (math.sqrt(5)+1)/2
        coordinates = []    #couldn't set up cyclic permutation
        for dummy in range(4):
            coordinates.append((0,(-1)**(dummy//2),(phi)*(-1)**(dummy)))
        for dummy in range(4):
            coordinates.append(((-1)**(dummy//2),(phi)*(-1)**(dummy),0))
        for dummy in range(4):
            coordinates.append(((phi)*(-1)**(dummy),0,(-1)**(dummy//2)))
        H = np.zeros((12,12))
        for r1 in coordinates:
            for r2 in coordinates:
                distance2 = 0
                distance2 =((r1[0])-(r2[0]))**2+((r1[1])-(r2[1]))**2+((r1[2])-(r2[2]))**2
                if round(distance2,6)==round((2)**2,6):
                    H[coordinates.index(r1),coordinates.index(r2)] = 1
        return H

    def fullerene(self): #treating as a truncated icosahedron
        phi = (math.sqrt(5)+1)/2
        coordinates = []
        #sets of 0, 1, phi
        for dummy in range(4):
            coordinates.append((0,(-1)**(dummy//2),(3*phi)*(-1)**(dummy)))
        for dummy in range(4):
            coordinates.append(((-1)**(dummy//2),(3*phi)*(-1)**(dummy),0))
        for dummy in range(4):
            coordinates.append(((3*phi)*(-1)**(dummy),0,(-1)**(dummy//2)))
        #sets of 1, 2+phi, 2*phi
        for dummy in range(8):
            coordinates.append(((1)*(-1)**(dummy//4),(2+phi)*(-1)**(dummy//2),((2*phi)*(-1)**(dummy))))
        for dummy in range(8):
            coordinates.append(((2+phi)*(-1)**(dummy//4),(2*phi)*(-1)**(dummy//2),((1)*(-1)**(dummy))))
        for dummy in range(8):
            coordinates.append(((2*phi)*(-1)**(dummy//4),(1)*(-1)**(dummy//2),((2+phi)*(-1)**(dummy))))
         #sets of phi, 2. 2*phi+1
        for dummy in range(8):
            coordinates.append(((phi)*(-1)**(dummy//4),(2)*(-1)**(dummy//2),((2*phi+1)*(-1)**(dummy))))
        for dummy in range(8):
            coordinates.append(((2*phi+1)*(-1)**(dummy//4),(phi)*(-1)**(dummy//2),((2)*(-1)**(dummy))))
        for dummy in range(8):
            coordinates.append(((2)*(-1)**(dummy//4),(2*phi+1)*(-1)**(dummy//2),((phi)*(-1)**(dummy))))
        #create matrix
        H = np.zeros((60,60))
        for r1 in coordinates:
            for r2 in coordinates:
                distance2 = 0

                distance2 =((r1[0])-(r2[0]))**2+((r1[1])-(r2[1]))**2+((r1[2])-(r2[2]))**2
                # check if distance equals 2 - a bond length in given coordinates
                if round(distance2,6)==round((2)**2,6):
                    H[coordinates.index(r1),coordinates.index(r2)] = 1

        return H


root = tk.Tk()
root.title('Huckinator')
app = Application(master=root)

root.mainloop()

#%%
