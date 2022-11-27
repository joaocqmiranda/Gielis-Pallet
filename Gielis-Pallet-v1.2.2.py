import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backend_tools import ToolBase
import tkinter
import pandas as pd
from tkinter import *
from PIL import ImageTk, Image

ws = Tk()
ws.title("Gielis Pallet")
# ws.geometry('850x500')

# Create an object of tkinter ImageTk

image = Image.open(r"C:\Users\joaoc\AppData\Local\Programs\Python\Python39\equacoes.png")
new_image = image.resize((450, 180))
photo = ImageTk.PhotoImage(new_image)

# Create an image label
img_label = Label(image=photo)
img_label.image = photo

img_label.grid(rowspan=50, columnspan=5, row=21, column=0)


# especificações de Tabela
lab_Graf = Label(ws, text="Especificações de Tabela").grid(row=0, column=1)

Tipo = IntVar()
lab_Graf = Label(ws, text="Tipo de Gráfico").grid(row=1, column=0)
op_2D = Radiobutton(ws, text="2D", indicatoron=0, variable=Tipo, value=2,
                    command=lambda: [lab_reso2.grid_forget(), reso2.grid_forget()])
op_2D.grid(row=1, column=2)
op_3D = Radiobutton(ws, text="3D", indicatoron=0, variable=Tipo, value=3,
                    command=lambda: [lab_reso2.grid(row=5, column=2), reso2.grid(row=5, column=3)])
op_3D.grid(row=1, column=3)

VisiEixos = IntVar()
lab_Exv = Label(ws, text="Eixos Visiveis?").grid(row=2, column=0)
op_evyes = Radiobutton(ws, text="Sim", indicatoron=0, variable=VisiEixos, value=1)
op_evyes.grid(row=2, column=2)
op_evno = Radiobutton(ws, text="Não", indicatoron=0, variable=VisiEixos, value=0)
op_evno.grid(row=2, column=3)

lab_Ng = Label(ws, text="Nº de Gráficos").grid(row=3, column=0)
Ng = Entry(ws, width=5)
Ng.grid(row=3, column=1)
Ng.insert(0, "4")

lab_Mat = Label(ws, text="Formato: Nº Linhas").grid(row=4, column=0)
Matriz1 = Entry(ws, width=5)
Matriz1.grid(row=4, column=1)
Matriz1.insert(0, "2")
lab_Mat = Label(ws, text="Nº Colunas").grid(row=4, column=2)
Matriz2 = Entry(ws, width=5)
Matriz2.grid(row=4, column=3)
Matriz2.insert(0, "2")

# Parametros da Grelha
Label(ws, text="Param/Legenda").grid(row=0, column=6)
Label(ws, text="Sx     ").grid(row=1, column=6)
Label(ws, text="Sy     ").grid(row=2, column=6)
Label(ws, text="Sz     ").grid(row=3, column=6)
Label(ws, text="\u03A62i [º]").grid(row=4, column=6)
Label(ws, text="\u03A62f [º]").grid(row=5, column=6)
Label(ws, text="\u03B11 [º]").grid(row=6, column=6)
Label(ws, text="\u03B12 [º]").grid(row=7, column=6)
Label(ws, text="a      ").grid(row=8, column=6)
Label(ws, text="b      ").grid(row=9, column=6)
Label(ws, text="m1      ").grid(row=10, column=6)
Label(ws, text="m2     ").grid(row=11, column=6)
Label(ws, text="n1     ").grid(row=12, column=6)
Label(ws, text="n2     ").grid(row=13, column=6)
Label(ws, text="n3     ").grid(row=14, column=6)
Label(ws, text="\u03A61i [º]").grid(row=15, column=6)
Label(ws, text="\u03A61f [º]").grid(row=16, column=6)
Label(ws, text="\u03B13 [º]").grid(row=17, column=6)
Label(ws, text="\u03B14 [º]").grid(row=18, column=6)
Label(ws, text="c      ").grid(row=19, column=6)
Label(ws, text="d      ").grid(row=20, column=6)
Label(ws, text="m3      ").grid(row=21, column=6)
Label(ws, text="m4     ").grid(row=22, column=6)
Label(ws, text="n4     ").grid(row=23, column=6)
Label(ws, text="n5     ").grid(row=24, column=6)
Label(ws, text="n6     ").grid(row=25, column=6)


def grelha():
    agt = Label(ws, text="Aviso! Nº Graf> Tabela")
    agt2 = Label(ws, text="Nº de Gráficos certos")
    if int(Matriz2.get()) * int(Matriz1.get()) < int(Ng.get()):
        agt2.grid_forget()
        agt.grid(row=3, column=4)
    else:
        agt.grid_forget()
        agt2.grid(row=3, column=4)

    Matriz = [int(Matriz1.get()), int(Matriz2.get())]

    Sx = [1] * int(Ng.get())
    Sy = [1] * int(Ng.get())
    Sz = [1] * int(Ng.get())
    phi2i = [1] * int(Ng.get())
    phi2f = [1] * int(Ng.get())
    alpha1 = [1] * int(Ng.get())
    alpha2 = [1] * int(Ng.get())
    a = [1] * int(Ng.get())
    b = [1] * int(Ng.get())
    m = [1] * int(Ng.get())
    m2 = [1] * int(Ng.get())
    n1 = [1] * int(Ng.get())
    n2 = [1] * int(Ng.get())
    n3 = [1] * int(Ng.get())
    Phi1i = [1] * int(Ng.get())
    Phi1f = [1] * int(Ng.get())
    Alpha1 = [1] * int(Ng.get())
    Alpha2 = [1] * int(Ng.get())
    A = [1] * int(Ng.get())
    B = [1] * int(Ng.get())
    M = [1] * int(Ng.get())
    M2 = [1] * int(Ng.get())
    N1 = [1] * int(Ng.get())
    N2 = [1] * int(Ng.get())
    N3 = [1] * int(Ng.get())

    tLeg = [1] * int(Ng.get())
    tSx = [1] * int(Ng.get())
    tSy = [1] * int(Ng.get())
    tSz = [1] * int(Ng.get())
    tphi2i = [1] * int(Ng.get())
    tphi2f = [1] * int(Ng.get())
    talpha1 = [1] * int(Ng.get())
    talpha2 = [1] * int(Ng.get())
    ta = [1] * int(Ng.get())
    tb = [1] * int(Ng.get())
    tm = [1] * int(Ng.get())
    tm2 = [1] * int(Ng.get())
    tn1 = [1] * int(Ng.get())
    tn2 = [1] * int(Ng.get())
    tn3 = [1] * int(Ng.get())
    tPhi1i = [1] * int(Ng.get())
    tPhi1f = [1] * int(Ng.get())
    tAlpha1 = [1] * int(Ng.get())
    tAlpha2 = [1] * int(Ng.get())
    tA = [1] * int(Ng.get())
    tB = [1] * int(Ng.get())
    tM = [1] * int(Ng.get())
    tM2 = [1] * int(Ng.get())
    tN1 = [1] * int(Ng.get())
    tN2 = [1] * int(Ng.get())
    tN3 = [1] * int(Ng.get())

    for i in range(int(Ng.get())):
        tLeg[i] = Entry(ws, width=10)
        tLeg[i].grid(row=0, column=i + 7)
        tLeg[i].insert(0, " ")
        tSx[i] = Entry(ws, width=10)
        tSx[i].grid(row=1, column=i + 7)
        tSx[i].insert(0, "1")
        tSy[i] = Entry(ws, width=10)
        tSy[i].grid(row=2, column=i + 7)
        tSy[i].insert(0, "1")
        tSz[i] = Entry(ws, width=10)
        tSz[i].grid(row=3, column=i + 7)
        tSz[i].insert(0, "1")

        tphi2i[i] = Entry(ws, width=10)
        tphi2i[i].grid(row=4, column=i + 7)
        tphi2i[i].insert(0, "-180")
        tphi2f[i] = Entry(ws, width=10)
        tphi2f[i].grid(row=5, column=i + 7)
        tphi2f[i].insert(0, "180")

        talpha1[i] = Entry(ws, width=10)
        talpha1[i].grid(row=6, column=i + 7)
        talpha1[i].insert(0, "0")
        talpha2[i] = Entry(ws, width=10)
        talpha2[i].grid(row=7, column=i + 7)
        talpha2[i].insert(0, "0")
        ta[i] = Entry(ws, width=10)
        ta[i].grid(row=8, column=i + 7)
        ta[i].insert(0, "1")
        tb[i] = Entry(ws, width=10)
        tb[i].grid(row=9, column=i + 7)
        tb[i].insert(0, "1")
        tm[i] = Entry(ws, width=10)
        tm[i].grid(row=10, column=i + 7)
        tm[i].insert(0, "1")
        tm2[i] = Entry(ws, width=10)
        tm2[i].grid(row=11, column=i + 7)
        tm2[i].insert(0, "1")
        tn1[i] = Entry(ws, width=10)
        tn1[i].grid(row=12, column=i + 7)
        tn1[i].insert(0, "1")
        tn2[i] = Entry(ws, width=10)
        tn2[i].grid(row=13, column=i + 7)
        tn2[i].insert(0, "1")
        tn3[i] = Entry(ws, width=10)
        tn3[i].grid(row=14, column=i + 7)
        tn3[i].insert(0, "1")
        tPhi1i[i] = Entry(ws, width=10)
        tPhi1i[i].grid(row=15, column=i + 7)
        tPhi1i[i].insert(0, "-90")
        tPhi1f[i] = Entry(ws, width=10)
        tPhi1f[i].grid(row=16, column=i + 7)
        tPhi1f[i].insert(0, "90")
        tAlpha1[i] = Entry(ws, width=10)
        tAlpha1[i].grid(row=17, column=i + 7)
        tAlpha1[i].insert(0, "0")
        tAlpha2[i] = Entry(ws, width=10)
        tAlpha2[i].grid(row=18, column=i + 7)
        tAlpha2[i].insert(0, "0")
        tA[i] = Entry(ws, width=10)
        tA[i].grid(row=19, column=i + 7)
        tA[i].insert(0, "1")
        tB[i] = Entry(ws, width=10)
        tB[i].grid(row=20, column=i + 7)
        tB[i].insert(0, "1")
        tM[i] = Entry(ws, width=10)
        tM[i].grid(row=21, column=i + 7)
        tM[i].insert(0, "1")
        tM2[i] = Entry(ws, width=10)
        tM2[i].grid(row=22, column=i + 7)
        tM2[i].insert(0, "1")
        tN1[i] = Entry(ws, width=10)
        tN1[i].grid(row=23, column=i + 7)
        tN1[i].insert(0, "1")
        tN2[i] = Entry(ws, width=10)
        tN2[i].grid(row=24, column=i + 7)
        tN2[i].insert(0, "1")
        tN3[i] = Entry(ws, width=10)
        tN3[i].grid(row=25, column=i + 7)
        tN3[i].insert(0, "1")

    def grafico():
        for i in range(int(Ng.get())):

            Sx[i] = (float(tSx[i].get()))
            Sy[i] = (float(tSy[i].get()))
            Sz[i] = (float(tSz[i].get()))
            phi2i[i] = (float(tphi2i[i].get()))*(np.pi)/180
            phi2f[i] = (float(tphi2f[i].get()))*(np.pi)/180
            alpha1[i] = (float(talpha1[i].get()))*(np.pi)/180
            alpha2[i] = (float(talpha2[i].get()))*(np.pi)/180
            a[i] = (float(ta[i].get()))
            b[i] = (float(tb[i].get()))
            m[i] = (float(tm[i].get()))
            m2[i] = (float(tm2[i].get()))
            n1[i] = (float(tn1[i].get()))
            n2[i] = (float(tn2[i].get()))
            n3[i] = (float(tn3[i].get()))
            Phi1i[i] = (float(tPhi1i[i].get())) * (np.pi) / 180
            Phi1f[i] = (float(tPhi1f[i].get())) * (np.pi) / 180
            Alpha1[i] = (float(tAlpha1[i].get()))*(np.pi)/180
            Alpha2[i] = (float(tAlpha2[i].get()))*(np.pi)/180
            A[i] = (float(tA[i].get()))
            B[i] = (float(tB[i].get()))
            M[i] = (float(tM[i].get()))
            M2[i] = (float(tM2[i].get()))
            N1[i] = (float(tN1[i].get()))
            N2[i] = (float(tN2[i].get()))
            N3[i] = (float(tN3[i].get()))

        if len(a) > ((Matriz[0]) * (Matriz[1])):
            tkinter.messagebox.showwarning(title="ALERTA", message="EXIXTTEM FIGURAS POR MOSTRAR")

        if int(Tipo.get()) == 3:
            angle = np.linspace(Phi1i[i], Phi1f[i], int(reso.get()))  # Max rec. 141 em cada e ser numero impar por causa da simetria
            angle2 = np.linspace(phi2i[i], phi2f[i], int(reso2.get()))
            theta, phi = np.meshgrid(angle2, angle)

            if A != 0:
                lab_ia = Label(ws, text="IRREGULAR")
                lab_ia.grid(row=4, column=5)
                fig = plt.figure()
                i = -1
                for i1 in range(Matriz[0]):
                    for i2 in range(Matriz[1]):
                        i = i + 1
                        if i <= len(a) - 1:
                            def r1(fi):
                                r11 = 1 / (
                                        (((abs(np.cos(m[i] * (fi - alpha1[i]) / 4) / (a[i]))) ** (n2[i])) + (
                                                (abs(np.sin(m2[i] * (fi - alpha2[i]) / 4) / (b[i]))) ** (
                                            n3[i]))) ** (
                                                1 / n1[i]))
                                if int(Logaritr.get()) == 1:
                                    r11 = r11 * (((angle + 3 * np.pi) / (2 * np.pi)) ** float(ir.get()))
                                return r11

                            def r2(fi):
                                return 1 / (
                                        (((abs(np.cos(M[i] * (fi - Alpha1[i]) / 4) / (A[i]))) ** (N2[i])) + (
                                                (abs(np.sin(M2[i] * (fi - Alpha2[i]) / 4) / (B[i]))) ** (
                                            N3[i]))) ** (
                                                1 / N1[i]))

                            X = r1(theta) * np.cos(theta) * r2(phi) * np.cos(phi) * Sx[i]
                            Y = r1(theta) * np.sin(theta) * r2(phi) * np.cos(phi) * Sy[i]
                            Z = r2(phi) * np.sin(phi) * Sz[i]
                            if int(toroidal.get()) == 1:
                                X = np.cos(theta) * (r1(theta) + r2(phi) * np.cos(phi)) * Sx[i]
                                Y = np.sin(theta) * (r1(theta) + r2(phi) * np.cos(phi)) * Sy[i]
                            if int(Logaritz.get()) == 1 or int(Logaritz.get()) == 3:
                                Z = Z * (((angle + 3 * np.pi) / (2 * np.pi)) ** float(ia.get()))

                            if int(Logaritz.get()) == 2 or int(Logaritz.get()) == 3:
                                Z = Z - r2(np.pi / 2) * np.sin(np.pi / 2) + r2(np.pi / 2) * np.sin(np.pi / 2) * Sz[
                                    i] * (((angle + 3 * np.pi) / (2 * np.pi)) ** float(id.get()))

                            ax = fig.add_subplot(Matriz[0], Matriz[1], i + 1,
                                                 projection='3d')  # fig.gca(projection = '3d')
                            maxs = max(Sx[i], Sy[i], Sz[
                                i]) * float(scale.get())  # max([max(X),b,max(Y),abs(min(Y)),max(Z),abs(min(Z))]) Não funciona por ser 3d
                            ax.set_xlim3d(-maxs, maxs)
                            ax.set_ylim3d(-maxs, maxs)
                            ax.set_zlim3d(-maxs, maxs)
                            ax.plot_surface(X, Y, Z, color=str(cl3.get()), rstride=1, cstride=1)
                            ax.set_title(str((tLeg[i].get())))

                            if int(VisiEixos.get()) == 0:
                                ax.axis('off')

            else:
                lab_ia.grid_forget()

                fig = plt.figure()
                i = -1
                for i1 in range(Matriz[0]):
                    for i2 in range(Matriz[1]):
                        i = i + 1
                        if i <= len(a) - 1:
                            def r(fi):
                                return 1 / ((((abs(np.cos(m[i] * fi / 4) / a[i])) ** (n2[i])) + (
                                        (abs(np.sin(m[i] * fi / 4) / b[i])) ** (n3[i]))) ** (1 / n1[i]))

                            X = r(theta) * np.cos(theta) * r(phi) * np.cos(phi) * Sx[i]
                            Y = r(theta) * np.sin(theta) * r(phi) * np.cos(phi) * Sy[i]
                            Z = r(phi) * np.sin(phi) * Sy[i]
                            ax = fig.add_subplot(Matriz[0], Matriz[1], i + 1,
                                                 projection='3d')  # fig.gca(projection = '3d')
                            maxs = max(Sx[i], Sy[i], Sz[
                                i])* float(scale.get())  # max([max(X),b,max(Y),abs(min(Y)),max(Z),abs(min(Z))]) Não funciona por ser 3d
                            ax.set_xlim3d(-maxs, maxs)
                            ax.set_ylim3d(-maxs, maxs)
                            ax.set_zlim3d(-maxs, maxs)
                            ax.plot_surface(X, Y, Z, color=str(cl3.get()), rstride=1, cstride=1)
                            ax.set_title(str((tLeg[i].get())))
                            if int(VisiEixos.get()) == 0:
                                ax.axis('off')

        if int(Tipo.get()) == 2:
            theta = np.linspace(phi2i[i], phi2f[i], int(reso.get()))
            fig, axs = plt.subplots(Matriz[0], Matriz[1])  # , subplot_kw=dict(projection='polar')
            i = -1
            for i1 in range(Matriz[0]):
                for i2 in range(Matriz[1]):
                    i = i + 1
                    if i <= len(a) - 1:
                        def r(fi):
                            return 1 / ((((abs(np.cos(m[i] * fi / 4) / a[i])) ** (n2[i])) + (
                                    (abs(np.sin(m[i] * fi / 4) / b[i])) ** (n3[i]))) ** (1 / n1[i]))

                        X = r(theta) * np.cos(theta)
                        Y = r(theta) * np.sin(theta)
                        axs[i1, i2].plot(X, Y, 'tab:'+str(cl2.get()))
                        axs[i1, i2].set_title(str((tLeg[i].get())))
                        maxs = max(Sx[i], Sy[i])* float(scale.get())
                        # axs[i1, i2].set(xlim=(min(X)-0.5, max(X)+0.5), ylim=(min(Y)-0.5, max(Y)+0.5))
                        axs[i1, i2].set(xlim=(-maxs, maxs), ylim=(-maxs, maxs))  # (xlim=(-3.6, 3.6), ylim=(-3.6, 3.6))
                        # axs[i1, i2].gca().set_aspect('equal', adjustable='box')#NOvo testar
                        if VisiEixos.get() == 0:
                            axs[i1, i2].axis('off')

        # plt.savefig('demo.png', transparent=True)
        matplotlib.rcParams["savefig.directory"] = str(patht.get())
        if int(Transp.get())==0:
            matplotlib.rcParams["savefig.transparent"] = "True"
        else:
            matplotlib.rcParams["savefig.transparent"] = "False"
        matplotlib.rcParams["savefig.format"] = str(formt.get())
        matplotlib.rcParams["savefig.dpi"] = str(dpi.get())
        # Display the mesh

        plt.show()

    Radiobutton(ws, text="Criar Gráficos", indicatoron=0, command=lambda: grafico()).grid(row=3, column=5)


Radiobutton(ws, text="Criar Grelha", indicatoron=0, command=lambda: grelha()).grid(row=3, column=3)

# Verifica se o numero de formas cabe na tabela
lab_reso = Label(ws, text="Resolução r1").grid(row=5, column=0)
reso = Entry(ws, width=10)
reso.grid(row=5, column=1)
reso.insert(0, "141")
lab_reso2 = Label(ws, text="Resolução r2")
lab_reso2.grid(row=5, column=2)
lab_reso2.grid_forget()
reso2 = Entry(ws, width=10)
reso2.grid(row=5, column=3)
reso2.insert(0, "141")
reso2.grid_forget()

# Parametros de Forma

lab_Graf = Label(ws, text="Parametros de Forma").grid(row=7, column=1)

toroidal = IntVar()
lab_Tor = Label(ws, text="Toroidal").grid(row=8, column=0)
op_tyes = Radiobutton(ws, text="Sim", indicatoron=0, variable=toroidal, value=1).grid(row=8, column=2)
op_tno = Radiobutton(ws, text="Não", indicatoron=0, variable=toroidal, value=0).grid(row=8, column=3)

Logaritr = IntVar()
lab_Exr = Label(ws, text="Exponencial radial?").grid(row=9, column=0)
op_eryes = Radiobutton(ws, text="Sim", indicatoron=0, variable=Logaritr, value=1,
                       command=lambda: [lab_ir.grid(row=11, column=1), ir.grid(row=12, column=1)]).grid(row=9, column=2)
op_erno = Radiobutton(ws, text="Não", indicatoron=0, variable=Logaritr, value=0,
                      command=lambda: [lab_ir.grid_forget(), ir.grid_forget()]).grid(row=9, column=3)

Logaritz = IntVar()
lab_Exz = Label(ws, text="Exponencial em z?").grid(row=10, column=0)
op_ezno = Radiobutton(ws, text="Não", indicatoron=0, variable=Logaritz, value=0,
                      command=lambda: [lab_ia.grid_forget(), ia.grid_forget(),
                                       lab_id.grid_forget(), id.grid_forget()]).grid(row=10, column=2)
op_ezamp = Radiobutton(ws, text="Amplitude", indicatoron=0, variable=Logaritz, value=1,
                       command=lambda: [lab_ia.grid(row=11, column=2), ia.grid(row=12, column=2),
                                        lab_id.grid_forget(), id.grid_forget()]).grid(row=10, column=3)
op_ezdesl = Radiobutton(ws, text="Deslocamento", indicatoron=0, variable=Logaritz, value=2,
                        command=lambda: [lab_id.grid(row=11, column=3), id.grid(row=12, column=3),
                                         lab_ia.grid_forget(), ia.grid_forget()]).grid(row=10, column=4)
op_exmist = Radiobutton(ws, text="Misto", indicatoron=0, variable=Logaritz, value=3,
                        command=lambda: [lab_id.grid(row=11, column=3), id.grid(row=12, column=3),
                                         lab_ia.grid(row=11, column=2), ia.grid(row=12, column=2)]).grid(row=10,
                                                                                                         column=5)
# Intensidades

lab_ir = Label(ws, text="Intensidade Radial")
lab_ir.grid(row=11, column=1)
lab_ir.grid_forget()
ir = Entry(ws, width=5)
ir.grid(row=12, column=1)
ir.insert(0, "1")
ir.grid_forget()

lab_ia = Label(ws, text="Intensidade Amplitude")
lab_ia.grid(row=11, column=2)
lab_ia.grid_forget()
ia = Entry(ws, width=5)
ia.grid(row=12, column=2)
ia.insert(0, "1")
ia.grid_forget()

lab_id = Label(ws, text="Intensidade Deslocamento")
lab_id.grid(row=11, column=3)
lab_id.grid_forget()
id = Entry(ws, width=5)
id.grid(row=12, column=3)
id.insert(0, "1")
id.grid_forget()

lab_tscale = Label(ws, text="Escala")
lab_tscale.grid(row=13, column=0)
scale = Entry(ws, width=5)
scale.grid(row=13, column=1)
scale.insert(0, "1")

Label(ws, text="Parametros da Imagem").grid(row=14, column=1)

lab_tpatht = Label(ws, text="Path")
lab_tpatht.grid(row=15, column=0)
patht = Entry(ws)
patht.grid(row=15, columnspan=2, sticky='e')
patht.insert(0, "C:Desktop/")

Transp = IntVar()
lab_Transp = Label(ws, text="Fundo").grid(row=16, column=0)
op_tpyes = Radiobutton(ws, text="Transparente", indicatoron=0, variable=Transp, value=0)
op_tpyes.grid(row=16, column=1)
op_tpno = Radiobutton(ws, text="Branco", indicatoron=0, variable=Transp, value=1)
op_tpno.grid(row=16, column=2)

lab_tformt = Label(ws, text="Formato")
lab_tformt.grid(row=17, column=0)
formt = Entry(ws, width=5)
formt.grid(row=17, column=1)
formt.insert(0, "svg")

lab_tdpi = Label(ws, text="DPI")
lab_tdpi.grid(row=18, column=0)
dpi = Entry(ws, width=5)
dpi.grid(row=18, column=1)
dpi.insert(0, "2400")

lab_cl2 = Label(ws, text="Cor Linha 2D")
lab_cl2.grid(row=19, column=0)
cl2 = Entry(ws, width=5)
cl2.grid(row=19, column=1)
cl2.insert(0, "grey")

lab_cl3 = Label(ws, text="Cor da Sup. 3D")
lab_cl3.grid(row=20, column=0)
cl3 = Entry(ws, width=5)
cl3.grid(row=20, columnspan=3)
cl3.insert(0, "w")

ws.mainloop()
