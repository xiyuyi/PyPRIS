import os
import pandas as pd
import numpy as np
import copy
from csv import writer
from sklearn.cluster import DBSCAN
import pickle
from bokeh.plotting import figure, output_file, show, output_notebook
import bokeh.io
from bokeh.layouts import gridplot


def pris_show(locs, candidate_locs, lb, outputhtml, saveoption, locs_g=None):
    locs=np.asarray(locs)
    blur = lb.b.reshape(128, 128)
    recb = lb.recb.reshape(128, 128)
    w=128
    h=128
    pw=400
    ph=400
    p1 = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur')
    p1.x_range.range_padding = p1.y_range.range_padding = 0
    p1.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")



    p2 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='reconstructed blur')
    p2.x_range.range_padding = p2.y_range.range_padding = 0
    p2.image(image=[recb], x=0, y=0, dw=w, dh=h, palette="Greys256")



    p3 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur + fits (red dots)')
    p3.x_range.range_padding = p3.y_range.range_padding = 0
    p3.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    p3.scatter(locs[:, 1], locs[:, 0], marker='dot', size=25, color='red')



    p4 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='residual')
    p4.x_range.range_padding = p4.y_range.range_padding = 0
    p4.image(image=[blur - recb], x=0, y=0, dw=w, dh=h, palette="Greys256")


    p5 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur + ground truth (green crosses)')
    p5.x_range.range_padding = p3.y_range.range_padding = 0
    p5.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    if locs_g is not None:
        p5.scatter(locs_g[:, 1], locs_g[:, 0], marker='cross', size=15, color='green', line_width=2)


    p6 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur + final fits (red dots) + ground truth (green crosses)')
    p6.x_range.range_padding = p4.y_range.range_padding = 0
    p6.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    p6.scatter(locs[:, 1], locs[:, 0], marker='dot', size=25, color='red')
    if locs_g is not None:
        p6.scatter(locs_g[:, 1], locs_g[:, 0], marker='cross', size=15, color='green', line_width=2)

    if saveoption:
        output_file(outputhtml, title="PRIS0to4 candidate pools")
    grid = gridplot([[p1, p2, p4], [p5, p3, p6]], plot_width=pw, plot_height=ph)
    show(grid)
    return

def bokeh_imshow(im):
    p1 = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                    title='blur')
    p1.x_range.range_padding = p1.y_range.range_padding = 0
    p1.image(image=[im], x=0, y=0, dw=500, dh=500, palette="Greys256")
    pw=450
    ph=450
    grid = gridplot([[p1]], plot_width=500, plot_height=500)
    show(grid)

def pris_show_pris0to4(locs, cands0, cands1, cands2, cands3, cands4, lb, outputhtml, saveoption):
    locs=np.asarray(locs)
    blur = lb.b.reshape(128, 128)
    recb = lb.recb.reshape(128, 128)
    w=128
    h=128
    pw=500
    ph=500
    x0 = cands0
    x1 = cands1
    x2 = cands2
    x3 = cands3
    x4 = cands4
    x5 = locs
    "show the support of pris0 and pris1"
    p1 = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris0 and pris1')
    p1.x_range.range_padding = p1.y_range.range_padding = 0
    p1.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        p1.scatter(x0[:, 2], x0[:, 1], marker='dot', size=25, color='red')
        p1.scatter(x1[:, 2], x1[:, 1], marker='dot', size=25, color='darkorange')
    except Warning:
        print("no pypris object detected.")
        pass



    p2 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris1 and pris2')
    p2.x_range.range_padding = p1.y_range.range_padding = 0
    p2.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        p2.scatter(x1[:, 2], x1[:, 1], marker='dot', size=25, color='darkorange')
        p2.scatter(x2[:, 2], x2[:, 1], marker='dot', size=25, color='green')
    except Warning:
        print("no pypris object detected.")
        pass



    p3 = figure(x_range=p2.x_range, y_range=p2.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris2 and pris3')
    p3.x_range.range_padding = p2.y_range.range_padding = 0
    p3.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        p3.scatter(x2[:, 2], x2[:, 1], marker='dot', size=25, color='green')
        p3.scatter(x3[:, 2], x3[:, 1], marker='dot', size=25, color='dodgerblue')
    except Warning:
        print("no pypris object detected.")
        pass



    p4 = figure(x_range=p2.x_range, y_range=p2.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris3 and pris4')
    p4.x_range.range_padding = p2.y_range.range_padding = 0
    p4.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        p4.scatter(x3[:, 2], x3[:, 1], marker='dot', size=25, color='dodgerblue')
        p4.scatter(x4[:, 2], x4[:, 1], marker='dot', size=25, color='blue')
    except Warning:
        print("no pypris object detected.")
        pass


    p5 = figure(x_range=p2.x_range, y_range=p2.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris4 and final fitting location')
    p5.x_range.range_padding = p2.y_range.range_padding = 0
    p5.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        p5.scatter(x4[:, 2], x4[:, 1], marker='dot', size=25, color='blue')
        p5.scatter(x5[:, 1], x5[:, 0], marker='dot', size=25, color='purple')
    except Warning:
        print("no pypris object detected.")
        pass

    p6 = figure(x_range=p2.x_range, y_range=p2.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='support of pris1/2/3/4, and final fitting locations.')
    p6.x_range.range_padding = p2.y_range.range_padding = 0
    p6.image(image=[blur], x=0, y=0, dw=w, dh=h, palette="Greys256")
    try:
        #p6.scatter(x0[:, 2], x0[:, 1], marker='dot', size=25, color='red')
        p6.scatter(x1[:, 2], x1[:, 1], marker='dot', size=25, color='darkorange')
        p6.scatter(x2[:, 2], x2[:, 1], marker='dot', size=25, color='green')
        p6.scatter(x3[:, 2], x3[:, 1], marker='dot', size=25, color='dodgerblue')
        p6.scatter(x4[:, 2], x4[:, 1], marker='dot', size=25, color='blue')
        p6.scatter(x5[:, 1], x5[:, 0], marker='dot', size=25, color='purple')
    except Warning:
        print("no pypris object detected.")
        pass




    if saveoption:
        output_file(outputhtml, title="result (PRIS4)")
    grid = gridplot([[p1, p2, p3], [p4, p5, p6]], plot_width=pw, plot_height=ph)
    show(grid)
    return

