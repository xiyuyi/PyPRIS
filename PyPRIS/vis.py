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


def pris_show(locs, lb, outputhtml, saveoption):
    locs=np.asarray(locs)
    blur = lb.b.reshape(128, 128)
    recb = lb.recb.reshape(128, 128)
    w=128
    h=128
    pw=450
    ph=450
    p1 = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur')
    p1.x_range.range_padding = p1.y_range.range_padding = 0
    p1.image(image=[blur], x=0, y=0, dw=w, dh=h)



    p2 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='reconstructed blur')
    p2.x_range.range_padding = p2.y_range.range_padding = 0
    p2.image(image=[recb], x=0, y=0, dw=w, dh=h)



    p3 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='blur + fits')
    p3.x_range.range_padding = p3.y_range.range_padding = 0
    p3.image(image=[blur], x=0, y=0, dw=w, dh=h)
    p3.scatter(locs[:, 1], locs[:, 0], marker='dot', size=25, color='red')



    p4 = figure(x_range=p1.x_range, y_range=p1.y_range,
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                title='residual')
    p4.x_range.range_padding = p4.y_range.range_padding = 0
    p4.image(image=[blur - recb], x=0, y=0, dw=w, dh=h)


    # p5 = figure(x_range=p1.x_range, y_range=p1.y_range,
    #             tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
    #             title='residual')
    # p5.x_range.range_padding = p4.y_range.range_padding = 0
    # p5.image(image=[blur - recb], x=0, y=0, dw=w, dh=h)
    #
    #
    # p6 = figure(x_range=p1.x_range, y_range=p1.y_range,
    #             tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
    #             title='residual')
    # p6.x_range.range_padding = p4.y_range.range_padding = 0
    # p6.image(image=[blur - recb], x=0, y=0, dw=w, dh=h)

    if saveoption:
        output_file(outputhtml, title="pris result")
    grid = gridplot([[p1, p2], [p3, p4]], plot_width=pw, plot_height=ph)
    show(grid)
    return

def bokeh_imshow(im):
    p1 = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                    title='blur')
    p1.x_range.range_padding = p1.y_range.range_padding = 0
    p1.image(image=[im], x=0, y=0, dw=500, dh=500)
    pw=450
    ph=450
    grid = gridplot([[p1]], plot_width=500, plot_height=500)
    show(grid)


