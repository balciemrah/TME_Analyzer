# -*- coding: utf-8 -*-
"""
Created on Monday Oct 23rd 2023
fMO.openFile(self)
fMO.openSeg(self)
fMO.openFolder(self)
fMO.LoadWorkspace(self)
fMO.LoadWorkspaceFolder(self)
fMO.ChangeChannelSetup(self)
fMO.cropFile(self)
fMO.ResetCurrent(self)
fMO.RemoveCurrent(self)
fMO.Reset(self)
"""

import tkinter
import os
from tkinter import ttk as ttkinter
# import tkinter.messagebox
import tkinter.filedialog
import skimage.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.backends.backend_tkagg as Tk_Agg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NT2Tk
import skimage.measure
from math import log
import pickle
from imageMenu import DestroyTK, popupmsg
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


# def DestroyTK(tk_param):
#     try:
#         tk_param.destroy()
#     except tkinter.TclError:
#         pass
#     except ttkinter.TclError:
#         pass

# def popupmsg(msg, self_control=True):
#     def Quit(*a):
#         popupnew2 = tkinter.Tk()
#         popupnew2.wm_title("!")
#         label2 = tkinter.Label(popupnew2, text="This will stop the current" +
#                                " operation following this iteration!\n" +
#                                "Are you sure?")
#         label2.pack(side="top", fill="x", pady=10)
#         B1 = tkinter.Button(popupnew2, text="Go Ahead",
#                             command=lambda: [DestroyTK(popupnew2),
#                                              DestroyTK(popupnew)])
#         B1.pack()
#         B2 = tkinter.Button(popupnew2, text="Go Back",
#                             command=lambda:[DestroyTK(popupnew2)])
#         B2.pack()
#         popupnew2.mainloop()
#     popupnew = tkinter.Tk()
#     popupnew.wm_title("!")
#     label = tkinter.Label(popupnew, text=msg)
#     label.pack(side="top", fill="x", pady=10)
#     if self_control:
#         B1 = tkinter.Button(popupnew, text="Okay", command=lambda:[DestroyTK(popupnew)])
#         B1.pack()
#     else:
#         popupnew.protocol("WM_DELETE_WINDOW", Quit)
#     popupnew.update()
#     return popupnew, label

def openFile(self):
    ftypes = [('Tiff images', '.tif .tiff')]
    filenames = tkinter.filedialog.askopenfilenames(
            parent=self.master, filetypes=ftypes)
    [popup_int, label_int] = popupmsg("...", False)
    for n_im, filename in enumerate(filenames):
        label_int['text'] = ("Loading image " +
                                str(n_im + 1) + " of " +
                                str(len(filenames)) +
                                " images.\n Please hold.")
        popup_int.update()
        im_raw = np.squeeze(np.float32(np.array(
            skimage.io.imread(filename))))
        channel_variable = []
        color_variable = []
        self.activeImage = len(self.FileDictionary)
        self.FileDictionary[len(self.FileDictionary)] = filename
        n_channels = min(im_raw.shape)
        Color_pointers = self.Color_info.copy()[0:n_channels]
        Channel_pointers = self.Markers.copy()[0:n_channels]
        if n_channels == im_raw.shape[0]:
            im_raw = im_raw.transpose(1, 2, 0)
        elif n_channels == im_raw.shape[1]:
            im_raw = im_raw.transpose(0, 2, 1)
        for i in range(np.size(Channel_pointers), n_channels):
            Channel_pointers.append("garbage")
        for i in range(np.size(Color_pointers), n_channels):
            Color_pointers.append("black")
        for i in range(n_channels):
            channel_variable.append([])
            color_variable.append([])
        self.activeImagePointer = []
        self.Channel_pointers.append(Channel_pointers)
        self.Color_pointers.append(Color_pointers)
        self.channel_variable.append(channel_variable)
        self.color_variable.append(color_variable)
        self.im_raw.append(im_raw)
        self.n_channels.append(n_channels)
        self.im_analyzed.append([])
        self.analyze_index.append([])
        self.foreground_threshold.append([])
        self.analysis_params.append(self.default_analysis_params)
        self.Cell_props.append([])
        self.Tissue_props.append({})
        self.small_images.append([])
    DestroyTK(popup_int)
    self.remake_side_window()

def openSeg(self):
    ftypes = [('Tiff images', '.tif .tiff')]
    filename = tkinter.filedialog.askopenfilename(
            parent=self.master, filetypes=ftypes)
    if True:
        analysis_params = self.analysis_params[self.active_image].copy()
        analysis_params.pop("Foreground")
        analysis_params["Foreground"] = {"thres": []}
        im_raw = np.squeeze(np.float32(np.array(
            skimage.io.imread(filename))))
        im_analyzed = self.im_analyzed[self.activeImage]
        analyze_index = self.analyze_index[self.activeImage]
        if np.size(im_raw.shape) > 2:
            n_channels = min(im_raw.shape)
            if n_channels == im_raw.shape[0]:
                im_raw = im_raw.transpose(1, 2, 0)
            elif n_channels == im_raw.shape[1]:
                im_raw = im_raw.transpose(0, 2, 1)
            Fore_mask = im_raw[:, :, 0] < 2
            Tumor_mask = im_raw[:, :, 0] == 0
            Stroma_mask = im_raw[:, :, 0] == 1
        else:
            Fore_mask = im_raw < 2
            Tumor_mask = im_raw == 0
            Stroma_mask = im_raw == 1
        if "Foreground" in analyze_index:
            im_analyzed[analyze_index.index("Foreground")] = Fore_mask
        else:
            im_analyzed.append(Fore_mask)
            analyze_index.append("Foreground")
        if "Tumor" in analyze_index:
            im_analyzed[analyze_index.index("Tumor")] = Tumor_mask
        else:
            im_analyzed.append(Tumor_mask)
            analyze_index.append("Tumor")
        if "Stroma" in analyze_index:
            im_analyzed[analyze_index.index("Stroma")] = Stroma_mask
        else:
            im_analyzed.append(Stroma_mask)
            analyze_index.append("Stroma")
        if np.size(im_raw.shape) > 2:
            if "DAPI" in analysis_params["Segments"]:
                analysis_params["Segments"].pop("DAPI")
            analysis_params["Segments"]["DAPI"] = {"thres": []}
            self.single_cells = skimage.measure.label(
                    im_raw[:, :, 1].astype('int32'))
            self.voronoi_image = skimage.measure.label(
                    im_raw[:, :, 1].astype('int32') +
                    im_raw[:, :, 2].astype('int32'))
            im_analyzed.append(im_raw[:, :, 1] > 0)
            analyze_index.append("DAPI")
            im_analyzed.append(self.single_cells)
            analyze_index.append("Nuclei")
            im_analyzed.append(self.voronoi_image)
            analyze_index.append("Cells")
            self.im_analyzed[self.activeImage] = im_analyzed
            self.analyze_index[self.activeImage] = analyze_index
            self.Get_cell_props()
        else:
            self.im_analyzed[self.activeImage] = im_analyzed
            self.analyze_index[self.activeImage] = analyze_index
        self.analysis_params[self.active_image] = analysis_params
        self.remake_side_window()

def openFolder(self):
    def Folder_stitch_selected(verts, *a):
        ROI_path = matplotlib.path.Path(verts)
        points = np.transpose((self.xy_locs[:, 0].ravel(),
                                self.xy_locs[:, 1].ravel()))
        mask = ROI_path.contains_points(points)
        mask = ~mask
        mask = np.nonzero(mask)[0]
        file_list = self.file_list
        filedir = self.filedir
        file_list_temp = self.file_list_temp
        for i in mask:
            file_list.pop(file_list.index(file_list_temp[i]))
        DestroyTK(self.popup)
        openFolder_real(filedir, file_list)
        openFolder_sub(1)

    def openFolder_sub(redo=0, *a):
        if redo == 0:
            ftypes = [('Tiff images', '.tif .tiff')]
            filename = tkinter.filedialog.askopenfilename(
                    parent=self.master, filetypes=ftypes)
            filedir = filename[:-1*filename[::-1].find('/')]
            file_list = os.listdir(filedir)
            self.filedir = filedir
        else:
            filedir = self.filedir
            file_list = os.listdir(filedir)
        self.file_list = file_list
        if self.combine_large == "PreDetermined stitch":
            filenames = []
            for i, filename in enumerate(file_list):
                if filename.find('_component_data.tif') >= 0:
                    filenames.append(filename)
            filenames_orig = filenames
            images_to_combine = []
            filenames = []
            master_file_list = []
            for i in range(len(filenames_orig)):
                filename = filenames_orig[i]
                filename = filename[:-1*filename[::-1].find('[')-1]
                filenames.append(filename)
                if filename != '':
                    if filenames[i] in master_file_list:
                        images_to_combine[master_file_list.index(
                                filenames[i])].append(filenames_orig[i])
                    else:
                        master_file_list.append(filenames[i])
                        images_to_combine.append([])
                        images_to_combine[master_file_list.index(
                                filenames[i])].append(filenames_orig[i])
            for i, file_list in enumerate(images_to_combine):
                if np.size(file_list) > 1:
                    xy_locs = []
                    for i, sub in enumerate(file_list):
                        xy_locs.append([
                                int(sub[sub.find('[')+1:sub.find(
                                        ',', sub.find('['))]),
                                int(sub[sub.find(',', sub.find(
                                        '['))+1:sub.find(']')])])
                    xy_locs = np.array(xy_locs)
                    xy_locs = (xy_locs - xy_locs.min(0))*self.magnification
            self.xy_locs = xy_locs
            popup = tkinter.Tk()
            popup.wm_title("Stitchin region selection")
            label = ttkinter.Label(
                    popup, text="Phenotype selection for file :\n" +
                    master_file_list[-1])
            label.pack(side="bottom", fill="x", pady=10)
            internal_windows = tkinter.Frame(popup, width=440, height=440)
            internal_windows.pack(side=tkinter.LEFT, anchor="w")
            f2 = plt.Figure(figsize=(6, 6), dpi=100)
            f2.patch.set_visible(False)
            f2.subplots_adjust(left=0.1, bottom=0.1, right=0.92, top=0.92,
                                wspace=0, hspace=0)
            ax2 = f2.gca()
            ax2 = f2.gca()
            image_canvas2 = Tk_Agg.FigureCanvasTkAgg(
                    f2, master=internal_windows)
            image_canvas2.draw()
            image_canvas2.get_tk_widget().pack(
                    side=tkinter.TOP, expand=True, anchor='n')
            image_canvas2._tkcanvas.pack(
                    side=tkinter.BOTTOM, fill=tkinter.BOTH, expand=True)
            image_toolbar2 = NT2Tk(
                    image_canvas2, internal_windows)
            image_toolbar2.update()
            ax2.clear()
            ax2.set_position([0.1, 0.1, 0.65, 0.65])
            ax2.scatter(xy_locs[:, 0], xy_locs[:, 1], c='C0')
            ax2.set_xlim((np.min(xy_locs[:, 0])-1000,
                            np.max(xy_locs[:, 0])+1000))
            ax2.set_ylim((np.min(xy_locs[:, 1])-1000,
                            np.max(xy_locs[:, 1])+1000))
            self.popup = popup
            ROIPolygon = matplotlib.widgets.PolygonSelector(
                    ax2, Folder_stitch_selected,
                    lineprops=dict(color='r', linestyle='-',
                                    linewidth=1, alpha=0.5),
                    markerprops=dict(marker='o', markersize=3,
                                        mec='r', mfc='r', alpha=0.5))
            ROIPolygon.active = True
            ROIPolygon.visible = True
            self.file_list_temp = file_list
            popup.mainloop()
        else:
            openFolder_real(filedir, file_list)

    def openFolder_real(filedir, file_list, *a):
        combine_large = self.combine_large in ["Combine images",
                                                "Combine and downsample"]
        downsample = self.combine_large in ["Stitch and downsample",
                                            "Combine and downsample"]
        stitch_large = self.combine_large in ["Stitch images",
                                                "Stitch and downsample",
                                                "PreDetermined stitch"]
        self.loaded_images = []

        for i, filename in enumerate(file_list):
            if filename.find('_component_data.tif') >= 0:
                im_raw = np.squeeze(np.float32(np.array(
                    skimage.io.imread(filedir + filename))))
                channel_variable = []
                color_variable = []
                w = []
                c = []
                self.activeImage = len(self.FileDictionary)
                self.FileDictionary[len(self.FileDictionary)] = filename
                n_channels = min(im_raw.shape)
                Color_pointers = self.Color_info.copy()[0:n_channels]
                Channel_pointers = self.Markers.copy()[0:n_channels]
                if n_channels == im_raw.shape[0]:
                    im_raw = im_raw.transpose(1, 2, 0)
                elif n_channels == im_raw.shape[1]:
                    im_raw = im_raw.transpose(0, 2, 1)
                for i in range(np.size(Channel_pointers), n_channels):
                    Channel_pointers.append("garbage")
                for i in range(np.size(Color_pointers), n_channels):
                    Color_pointers.append("black")
                for i in range(n_channels):
                    channel_variable.append([])
                    color_variable.append([])

                self.loaded_images.append(self.activeImage)
                self.activeImagePointer = []
                self.Channel_pointers.append(Channel_pointers)
                self.Color_pointers.append(Color_pointers)
                self.channel_variable.append(channel_variable)
                self.color_variable.append(color_variable)
                self.im_raw.append(im_raw)
                self.n_channels.append(n_channels)
                self.im_analyzed.append([])
                self.analyze_index.append([])
                self.foreground_threshold.append([])
                self.analysis_params.append(self.default_analysis_params)
                self.Cell_props.append([])
                self.Tissue_props.append({})
                self.small_images.append([])
                filename = filename[:filename.find(
                        '_component_data.tif')] + '_binary_seg_maps.tif'
                if filename in file_list:
                    analysis_params = self.analysis_params[
                            self.activeImage].copy()
                    im_raw = np.squeeze(np.float32(np.array(
                        skimage.io.imread(filedir + filename))))
                    im_analyzed = self.im_analyzed[self.activeImage]
                    analyze_index = self.analyze_index[self.activeImage]
                    analysis_params.pop("Foreground")
                    analysis_params["Foreground"] = {"thres": []}
                    if np.size(im_raw.shape) > 2:
                        n_channels = min(im_raw.shape)
                        if n_channels == im_raw.shape[0]:
                            im_raw = im_raw.transpose(1, 2, 0)
                        elif n_channels == im_raw.shape[1]:
                            im_raw = im_raw.transpose(0, 2, 1)
                        Fore_mask = im_raw[:, :, 0] < 2
                        Tumor_mask = im_raw[:, :, 0] == 0
                        Stroma_mask = im_raw[:, :, 0] == 1
                    else:
                        Fore_mask = im_raw < 2
                        Tumor_mask = im_raw == 0
                        Stroma_mask = im_raw == 1
                    if "Foreground" in analyze_index:
                        im_analyzed[analyze_index.index(
                                "Foreground")] = Fore_mask
                    else:
                        im_analyzed.append(Fore_mask)
                        analyze_index.append("Foreground")
                    if "Tumor" in analyze_index:
                        im_analyzed[analyze_index.index(
                                "Tumor")] = Tumor_mask
                    else:
                        im_analyzed.append(Tumor_mask)
                        analyze_index.append("Tumor")
                    if "Stroma" in analyze_index:
                        im_analyzed[analyze_index.index(
                                "Stroma")] = Stroma_mask
                    else:
                        im_analyzed.append(Stroma_mask)
                        analyze_index.append("Stroma")
                    if np.size(im_raw.shape) > 2:
                        if "DAPI" in analysis_params["Segments"]:
                            analysis_params["Segments"].pop("DAPI")
                        analysis_params["Segments"]["DAPI"] = {
                                "thres": []}
                        self.single_cells = skimage.measure.label(
                                im_raw[:, :, 1].astype('int32'))
                        self.voronoi_image = skimage.measure.label(
                                im_raw[:, :, 1].astype('int32') +
                                im_raw[:, :, 2].astype('int32'))
                        im_analyzed.append(im_raw[:, :, 1] > 0)
                        analyze_index.append("DAPI")
                        im_analyzed.append(self.single_cells)
                        analyze_index.append("Nuclei")
                        im_analyzed.append(self.voronoi_image)
                        analyze_index.append("Cells")
                        self.im_analyzed[self.activeImage] = im_analyzed
                        self.analyze_index[
                                self.activeImage] = analyze_index
                        self.Get_cell_props()
                    else:
                        self.im_analyzed[self.activeImage] = im_analyzed
                        self.analyze_index[
                                self.activeImage] = analyze_index
                    self.analysis_params[
                            self.activeImage] = analysis_params.copy()
        if (np.size(self.loaded_images) > 1) & combine_large:
            nx_size = np.uint32(np.floor(np.sqrt(np.size(
                    self.loaded_images))))
            ny_size = nx_size
            n_image = 0
            while np.size(self.loaded_images) > (nx_size * ny_size):
                ny_size = ny_size + 1
            im_temp = self.im_raw[self.loaded_images[n_image]]
            i = 2
            while (i < nx_size) & downsample:
                im_temp = ((skimage.transform.pyramid_reduce(
                        im_temp/im_temp.max()) * im_temp.max()))
                i = i * 2
            im_raw_combined = im_temp
            for i in range(1, nx_size):
                n_image = n_image + 1
                if n_image < np.size(self.loaded_images):
                    im_temp = self.im_raw[self.loaded_images[n_image]]
                    i = 2
                    while (i < nx_size) & downsample:
                        im_temp = ((skimage.transform.pyramid_reduce(
                                im_temp/im_temp.max()) * im_temp.max()))
                        i = i * 2
                else:
                    im_temp = np.zeros_like(im_temp)
                im_raw_combined = np.concatenate((
                        im_raw_combined, im_temp), axis=1)
            for i in range(1, ny_size):
                n_image = n_image + 1
                if n_image < np.size(self.loaded_images):
                    im_temp = self.im_raw[self.loaded_images[n_image]]
                    i = 2
                    while (i < nx_size) & downsample:
                        im_temp = ((skimage.transform.pyramid_reduce(
                                im_temp/im_temp.max()) * im_temp.max()))
                        i = i * 2
                else:
                    im_temp = np.zeros_like(im_temp)
                im_combined = im_temp
                for i in range(1, nx_size):
                    n_image = n_image + 1
                    if n_image < np.size(self.loaded_images):
                        im_temp = self.im_raw[self.loaded_images[n_image]]
                        i = 2
                        while (i < nx_size) & downsample:
                            im_temp = (skimage.transform.pyramid_reduce(
                                    im_temp/im_temp.max()) * im_temp.max())
                            i = i * 2
                    else:
                        im_temp = np.zeros_like(im_temp)
                    im_combined = np.concatenate((im_combined, im_temp),
                                                    axis=1)
                im_raw_combined = np.concatenate((im_raw_combined,
                                                    im_combined), axis=0)
            pheno_list = []
            for i in self.loaded_images:
                pheno_list = np.unique(np.concatenate((
                        pheno_list, self.analyze_index[i])))
            im_analyzed = []
            analyze_index = []
            for pheno in pheno_list:
                im_analyzed_temp = np.zeros_like(im_raw_combined[:, :, 0],
                                                    dtype="int64")
                for n_image in self.loaded_images:
                    if pheno in self.analyze_index[n_image]:
                        im_temp2 = self.im_analyzed[n_image][
                                self.analyze_index[n_image].index(pheno)]
                        i = 2
                        while (i < nx_size) & downsample:
                            im_temp2 = (skimage.transform.pyramid_reduce(
                                    im_temp2/im_temp2.max()) *
                                        im_temp2.max())
                            i = i * 2
                        n_y = np.uint32(n_image/nx_size)
                        n_x = n_image - nx_size * n_y
                        if pheno in ["Nuclei", "Cells"]:
                            im_temp2[im_temp2 > 0] = np.max(
                                    im_analyzed_temp)
                        im_analyzed_temp[n_y*im_temp2.shape[0]:(
                                n_y+1)*im_temp2.shape[0],
                                            n_x*im_temp2.shape[1]:(
                                n_x+1)*im_temp2.shape[1]] = im_temp2
                im_analyzed.append(im_analyzed_temp)
                analyze_index.append(pheno)
            im_raw = im_raw_combined
            self.activeImage = len(self.FileDictionary)
            self.FileDictionary[len(
                    self.FileDictionary)] = filedir + 'combined'
            DestroyTK(self.rightWindow)
            self.rightWindow = tkinter.Frame(self.rightWindow_master,
                                                width=100)
            self.rightWindow.pack(side=tkinter.RIGHT)
            channel_variable = []
            color_variable = []
            w = []
            c = []
            n_channels = min(im_raw.shape)
            Color_pointers = self.Color_info.copy()[0:n_channels]
            Channel_pointers = self.Markers.copy()[0:n_channels]
            # if n_channels == 3:
            #     Color_pointers[0] = 'red'
            #     Color_pointers[1] = 'green'
            #     Color_pointers[2] = 'blue'
            if n_channels == im_raw.shape[0]:
                im_raw = im_raw.transpose(1, 2, 0)
            elif n_channels == im_raw.shape[1]:
                im_raw = im_raw.transpose(0, 2, 1)
            for i in range(np.size(Channel_pointers), n_channels):
                Channel_pointers.append("garbage")
            for i in range(np.size(Color_pointers), n_channels):
                Color_pointers.append("black")
            for i in range(n_channels):
                channel_variable.append([])
                color_variable.append([])
                internal_windows = tkinter.Frame(self.rightWindow,
                                                    width=100, height=20)
                internal_windows.pack(side=tkinter.TOP)
                channel_variable[i] = tkinter.StringVar(
                        internal_windows)
                channel_variable[i].set(Channel_pointers[i])
                channel_variable[i].trace("w", self.cc_changed)
                w.append(tkinter.OptionMenu(internal_windows,
                                            channel_variable[i],
                                            *self.Markers))
                w[i].config(width=10)
                w[i].pack(side=tkinter.LEFT)
                color_variable[i] = tkinter.StringVar(internal_windows)
                color_variable[i].set(Color_pointers[i])
                color_variable[i].trace("w", self.cc_changed)
                c.append(tkinter.OptionMenu(internal_windows,
                                            color_variable[i],
                                            *self.Color_info))
                c[i].config(width=10)
                c[i].pack(side=tkinter.LEFT)
            internal_windows = tkinter.Frame(self.rightWindow,
                                                width=100, height=20)
            internal_windows.pack(side=tkinter.TOP)
            im_display_button = tkinter.Button(
                    master=internal_windows, text='Display',
                    command=self.display_composite_image)
            im_display_button.pack(side=tkinter.TOP)
            internal_windows = tkinter.Frame(self.rightWindow,
                                                width=100, height=20)
            internal_windows.pack(side=tkinter.TOP)
            active_image = tkinter.StringVar(internal_windows)
            active_image.set('(' + str(self.activeImage) + ', '
                                + self.FileDictionary[self.activeImage] + ')')
            active_image.trace("w", self.image_changed)
            w = tkinter.OptionMenu(internal_windows, active_image,
                                    *self.FileDictionary.items())
            w.config(width=10)
            w.pack(side=tkinter.TOP)

            self.activeImagePointer = active_image
#                self.activeImagePointer = []
            self.Channel_pointers.append(Channel_pointers)
            self.Color_pointers.append(Color_pointers)
            self.channel_variable.append(channel_variable)
            self.color_variable.append(color_variable)
            self.im_raw.append(im_raw)
            self.n_channels.append(n_channels)
            self.im_analyzed.append(im_analyzed)
            self.analyze_index.append(analyze_index)
            self.foreground_threshold.append([])
            self.analysis_params.append(self.default_analysis_params)
            self.Cell_props.append([])
            self.Tissue_props.append({})
            self.small_images.append([])

        if (np.size(self.loaded_images) > 1) & stitch_large:
            images_to_combine = []
            filenames = []
            master_file_list = []
            for i in range(np.size(self.loaded_images)):
                filename = self.FileDictionary[self.loaded_images[i]]
                filename = filename[:-1*filename[::-1].find('[')-1]
                filenames.append(filename)
                if filename != '':
                    if filenames[i] in master_file_list:
                        images_to_combine[master_file_list.index(filenames[
                                i])].append(self.loaded_images[i])
                    else:
                        images_to_combine.append([self.loaded_images[i]])
                        master_file_list.append(filenames[i])
            for i, file_list in enumerate(images_to_combine):
                if np.size(file_list) > 1:
                    xy_locs = []
                    nx_size = np.uint32(np.floor(
                            np.sqrt(np.size(file_list))))
                    i = 2
                    n_downsample = 1
                    while (i < nx_size) & downsample:
                        i = i * 2
                        n_downsample = n_downsample/2
                    for i, file_number in enumerate(file_list):
                        sub = self.FileDictionary[file_number]
                        xy_locs.append([int(sub[sub.find('[')+1:sub.find(
                                                ',', sub.find('['))]),
                                        int(sub[sub.find(',', sub.find(
                                                '['))+1:sub.find(']')])])
                    xy_locs = np.array(xy_locs)
                    xy_locs = (xy_locs - xy_locs.min(0))*self.magnification
                    im_temp = np.float32(np.array(self.im_raw[
                            file_number]))
                    im_temp = im_temp.transpose(
                            np.argsort(im_temp.shape)[2],
                            np.argsort(im_temp.shape)[1],
                            np.argsort(im_temp.shape)[0])
                    im_raw = np.zeros([
                            np.uint32(n_downsample * (
                                    100 + im_temp.shape[0] +
                                    xy_locs.max(0)[0])),
                            np.uint32(n_downsample * (
                                    100 + im_temp.shape[1] +
                                    xy_locs.max(0)[1])),
                            im_temp.shape[2]], 'float32')

                    for i, file_number in enumerate(file_list):
                        im_temp = np.float32(np.array(self.im_raw[
                                file_number]))
                        im_temp = im_temp.transpose(
                                np.argsort(im_temp.shape)[2],
                                np.argsort(im_temp.shape)[1],
                                np.argsort(im_temp.shape)[0])
                        j = 2
                        while (j < nx_size) & downsample:
                            j = j * 2
                            im_temp = ((
                                    skimage.transform.pyramid_reduce(
                                            im_temp/im_temp.max()) *
                                    im_temp.max()))
                        im_raw[np.uint32(n_downsample * xy_locs[
                                i, 0]):np.uint32(n_downsample * xy_locs[
                                        i, 0])+im_temp.shape[0],
                                np.uint32(n_downsample * xy_locs[
                                        i, 1]):np.uint32(
                                            n_downsample * xy_locs[
                                                    i, 1]) +
                                im_temp.shape[1], :] = im_temp
                    im_raw = im_raw.transpose(np.argsort(im_raw.shape)[2],
                                                np.argsort(im_raw.shape)[1],
                                                np.argsort(im_raw.shape)[0])

                    self.activeImage = len(self.FileDictionary)
                    self.FileDictionary[len(self.FileDictionary)] = sub[
                            :-1 * sub[::-1].find('[') - 1] + 'stitched'
                    DestroyTK(self.rightWindow)
                    self.rightWindow = tkinter.Frame(
                            self.rightWindow_master, width=100)
                    self.rightWindow.pack(side=tkinter.RIGHT)
                    channel_variable = []
                    color_variable = []
                    w = []
                    c = []
                    n_channels = min(im_raw.shape)
                    Color_pointers = self.Color_info.copy()[0:n_channels]
                    Channel_pointers = self.Markers.copy()[0:n_channels]
                    if n_channels == im_raw.shape[0]:
                        im_raw = im_raw.transpose(1, 2, 0)
                    elif n_channels == im_raw.shape[1]:
                        im_raw = im_raw.transpose(0, 2, 1)
                    for i in range(np.size(Channel_pointers), n_channels):
                        Channel_pointers.append("garbage")
                    for i in range(np.size(Color_pointers), n_channels):
                        Color_pointers.append("black")
                    for i in range(n_channels):
                        internal_windows = tkinter.Frame(self.rightWindow,
                                                            width=100,
                                                            height=20)
                        internal_windows.pack(side=tkinter.TOP)
                        channel_variable.append(tkinter.StringVar(
                                internal_windows))
                        channel_variable[i].set(Channel_pointers[i])
                        channel_variable[i].trace("w", self.cc_changed)
                        w.append(tkinter.OptionMenu(internal_windows,
                                                    channel_variable[i],
                                                    *self.Markers))
                        w[i].config(width=10)
                        w[i].pack(side=tkinter.LEFT)
                        color_variable.append(tkinter.StringVar(
                                internal_windows))
                        color_variable[i].set(Color_pointers[i])
                        color_variable[i].trace("w", self.cc_changed)
                        c.append(tkinter.OptionMenu(internal_windows,
                                                    color_variable[i],
                                                    *self.Color_info))
                        c[i].config(width=10)
                        c[i].pack(side=tkinter.LEFT)
                    internal_windows = tkinter.Frame(self.rightWindow,
                                                        width=100, height=20)
                    internal_windows.pack(side=tkinter.TOP)
                    im_display_button = tkinter.Button(
                            master=internal_windows, text='Display',
                            command=self.display_composite_image)
                    im_display_button.pack(side=tkinter.TOP)
                    internal_windows = tkinter.Frame(self.rightWindow,
                                                        width=100, height=20)
                    internal_windows.pack(side=tkinter.TOP)
                    active_image = tkinter.StringVar(internal_windows)
                    active_image.set('(' + str(self.activeImage) + ', '
                                        + self.FileDictionary[
                                                self.activeImage] + ')')
                    active_image.trace("w", self.image_changed)
                    w = tkinter.OptionMenu(internal_windows, active_image,
                                            *self.FileDictionary.items())
                    w.config(width=10)
                    w.pack(side=tkinter.TOP)

                    self.activeImagePointer = active_image
                    self.Channel_pointers.append(Channel_pointers)
                    self.Color_pointers.append(Color_pointers)
                    self.channel_variable.append(channel_variable)
                    self.color_variable.append(color_variable)
                    self.im_raw.append(im_raw)
                    self.n_channels.append(n_channels)
                    self.im_analyzed.append([])
                    self.analyze_index.append([])
                    self.foreground_threshold.append([])
                    self.analysis_params.append(
                            self.default_analysis_params)
                    self.Cell_props.append([])
                    self.Tissue_props.append({})
                    self.small_images.append([])
        self.remake_side_window()

    def openFolder_test(*a):
        stitch_large = self.combine_large in ["Stitch images",
                                                "Stitch and downsample",
                                                "PreDetermined stitch"]
        if stitch_large:
            def get_tick_values(*args):
                self.magnification = 1
                for i in range(dropdown_options.index(
                        combine_option.get())):
                    self.magnification = self.magnification * 2
            popup = tkinter.Tk()
            dropdown_options = ["10x",
                                "20x",
                                "40x"]
            combine_option = tkinter.StringVar(popup)
            combine_option.set(dropdown_options[0])
            combine_option.trace("w", get_tick_values)
            label = tkinter.Label(popup, text="The magnification is:")
            label.pack()
            c = tkinter.OptionMenu(popup, combine_option,
                                    *dropdown_options)
            c.pack()
            folderButton = ttkinter.Button(
                    popup, text="Load Folder", command=lambda: [
                            get_tick_values(), DestroyTK(popup),
                            openFolder_sub()])
            folderButton.pack()
            popup.mainloop()

        else:
            openFolder_sub()

    def get_tick_values(*args):
        self.combine_large = combine_option.get()
    popup = tkinter.Tk()
    dropdown_options = ["Do not combine images",
                        "Combine images",
                        "Combine and downsample",
                        "Stitch images",
                        "Stitch and downsample",
                        "PreDetermined stitch"]
    combine_option = tkinter.StringVar(popup)
    combine_option.set(dropdown_options[0])
    combine_option.trace("w", get_tick_values)
    c = tkinter.OptionMenu(popup, combine_option, *dropdown_options)
    c.pack()
    folderButton = ttkinter.Button(
            popup, text="Load Folder", command=lambda: [
                    get_tick_values(), DestroyTK(popup),
                    openFolder_test()])
    folderButton.pack()
    popup.mainloop()

def LoadWorkspace(self):
    filenames = tkinter.filedialog.askopenfilenames(parent=self.master)
    [popup_int, label_int] = popupmsg("...", False)
    for n_im, filename in enumerate(filenames):
        label_int['text'] = ("Loading workspace for " +
                                str(n_im + 1) + " of " +
                                str(len(filenames)) +
                                " images.\n Please hold.")
        popup_int.update()
        try:
            pickle_data = pickle.load(open(filename, "rb"))
            self.activeImage = len(self.FileDictionary)
            if isinstance(pickle_data["im_raw"], str):
                im_raw = np.squeeze(np.float32(np.array(
                    skimage.io.imread(pickle_data["im_raw"]))))
                self.im_raw.append(im_raw)
            else:
                self.im_raw.append(pickle_data["im_raw"])

            self.FileDictionary[len(self.FileDictionary)] = pickle_data[
                    "filename"]
            self.Channel_pointers.append(pickle_data["Channel_pointers"])
            self.Color_pointers.append(pickle_data["Color_pointers"])
            self.im_analyzed.append(pickle_data["im_analyzed"])
            self.analyze_index.append(pickle_data["analyze_index"])
            self.Cell_props.append(pickle_data["Cell_props"])
            self.analysis_params.append(pickle_data["analysis_params"])
            self.n_channels.append(min(self.im_raw[
                    self.activeImage].shape))
            if "Tissue_props" in pickle_data:
                self.Tissue_props.append(pickle_data["Tissue_props"])
            else:
                self.Tissue_props.append({})
            self.foreground_threshold.append(pickle_data[
                    "foreground_threshold"])
            self.small_images.append([])
            Color_pointers = self.Color_pointers[self.activeImage]
            channel_variable = []
            color_variable = []
            n_channels = self.n_channels[self.activeImage]
            for i in range(n_channels):
                channel_variable.append([])
                color_variable.append([])
            for i in range(n_channels, len(Color_pointers)):
                color_variable.append([])
            self.channel_variable.append(channel_variable)
            self.color_variable.append(color_variable)
            self.activeImagePointer = []
            self.remake_side_window()
        except NameError:
            pass
        except AttributeError:
            pass
    DestroyTK(popup_int)
    self.remake_side_window()

def LoadWorkspaceFolder(self):
    filename = tkinter.filedialog.askopenfilename(parent=self.master)
    filedir = filename[:-1*filename[::-1].find('/')]
    file_list = os.listdir(filedir)
    [popup_int, label_int] = popupmsg("...", False)
    for n_im, filename in enumerate(file_list):
        print(filename)
        label_int['text'] = ("Loading workspace for " +
                                str(n_im + 1) + " of " +
                                str(len(file_list)) +
                                " images.\n Please hold.")
        popup_int.update()
        try:
            filename = filedir + filename
            pickle_data = pickle.load(open(filename, "rb"))
            self.activeImage = len(self.FileDictionary)
            self.im_raw.append(pickle_data["im_raw"])
            self.FileDictionary[len(self.FileDictionary)] = pickle_data[
                    "filename"]
            self.Channel_pointers.append(pickle_data["Channel_pointers"])
            self.Color_pointers.append(pickle_data["Color_pointers"])
            self.im_analyzed.append(pickle_data["im_analyzed"])
            self.analyze_index.append(pickle_data["analyze_index"])
            self.Cell_props.append(pickle_data["Cell_props"])
            self.analysis_params.append(pickle_data["analysis_params"])
            if "Tissue_props" in pickle_data:
                self.Tissue_props.append(pickle_data["Tissue_props"])
            else:
                self.Tissue_props.append({})
            self.n_channels.append(min(self.im_raw[
                    self.activeImage].shape))
            self.foreground_threshold.append(pickle_data[
                    "foreground_threshold"])
            self.small_images.append([])
            print(pickle_data["analysis_params"])
            Color_pointers = self.Color_pointers[self.activeImage]
            channel_variable = []
            color_variable = []
            n_channels = self.n_channels[self.activeImage]
            for i in range(n_channels):
                channel_variable.append([])
                color_variable.append([])
            for i in range(n_channels, len(Color_pointers)):
                color_variable.append([])
            self.channel_variable.append(channel_variable)
            self.color_variable.append(color_variable)
            self.activeImagePointer = []
        except NameError:
            pass
        except AttributeError:
            pass
    DestroyTK(popup_int)
    self.remake_side_window()

def ChangeChannelSetup(self):
    def get_tick_values(*args):
        self.Markers = lookup_markers[dropdown_options.index(
                combine_option.get())]
        self.Color_info = lookup_colors[dropdown_options.index(
                combine_option.get())]
    lookup_markers = self.all_lookup_markers
    lookup_colors = self.all_lookup_colors
    popup = tkinter.Tk()
    dropdown_options = self.all_marker_dropdown
    combine_option = tkinter.StringVar(popup)
    combine_option.set(dropdown_options[0])
    combine_option.trace("w", get_tick_values)
    c = tkinter.OptionMenu(popup, combine_option, *dropdown_options)
    c.pack()
    folderButton = ttkinter.Button(popup, text="Confirm",
                                    command=lambda: [get_tick_values(),
                                                    DestroyTK(popup)])
    folderButton.pack()
    popup.mainloop()

def cropFile(self):
    def SelectCrop(*a):
        DestroyTK(self.popup_statusBar)
        t = ("Select the crop region with a rectangle.")
        self.popup_statusBar = tkinter.Label(self.popup, text=t, bd=1,
                                                relief=tkinter.SUNKEN,
                                                anchor=tkinter.W)
        self.popup_statusBar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
        self.CropRectangle = matplotlib.widgets.RectangleSelector(
                self.ax, onCropSelect,
                lineprops=dict(color='w', linestyle='-',
                                linewidth=2, alpha=0.5),
                marker_props=dict(marker='o', markersize=7,
                                    mec='w', mfc='w', alpha=0.5))
        self.CropRectrangle.set_active(True)
        if self.showCropmessage == 1:
            self.showCropmessage = 0
            popupmsg(t)

    def onCropSelect(eclick, erelease, *a):
        vertices = [eclick.xdata, erelease.xdata, eclick.ydata,
                    erelease.ydata]
        vertices = [np.max([np.int32(np.round(x)), 0]) for x in vertices]
        self.CropRectangle.set_active(False)
        mask = np.zeros_like(self.reduced_image[:, :, 0])
        vertices[:2] = [np.min([x, mask.shape[1]]) for x in vertices[:2]]
        vertices[2:] = [np.min([x, mask.shape[0]]) for x in vertices[2:]]
        mask[vertices[2]:vertices[3], vertices[0]:vertices[1]] = 1
        self.activeCrop = mask
        ax = self.ax
        if len(self.activeCrop) != 0:
            im_temp = self.imCropImage
            mask_copy = np.float32(self.activeCrop)
            mask_copy = (1 - mask_copy)/2
            for i in range(3):
                im_temp[:, :, i] = mask_copy + im_temp[:, :, i]
            ax.clear()
            ax.imshow(im_temp, aspect='equal')
            ax.autoscale(False)
            ax.axis('off')
            for i in range(3):
                im_temp[:, :, i] = im_temp[:, :, i] - mask_copy
                self.ax_canvas.draw()

    def addCrop(*a):
        mask = self.activeCrop
        im_raw = self.im_raw[self.activeImage]
        while round(log(mask.size/im_raw[:, :, 0].size, 2)) < 0:
            mask = skimage.transform.pyramid_expand(mask)
        a = im_raw[:, :, 0][np.ix_(mask.any(1), mask.any(0))]
        a = np.zeros([a.shape[0], a.shape[1], im_raw.shape[2]],
                        dtype=im_raw.dtype)
        for i in range(im_raw.shape[2]):
            a[:, :, i] = im_raw[:, :, i][np.ix_(mask.any(1), mask.any(0))]
        im_raw = a
        self.FileDictionary[len(
                self.FileDictionary)] = self.FileDictionary[
                                            self.activeImage] + "- cropped"
        self.Channel_pointers.append(
                self.Channel_pointers[self.activeImage])
        self.Color_pointers.append(self.Color_pointers[self.activeImage])
        self.channel_variable.append(
                self.channel_variable[self.activeImage])
        self.color_variable.append(self.color_variable[self.activeImage])
        self.im_raw.append(im_raw)
        self.n_channels.append(self.n_channels[self.activeImage])
        self.im_analyzed.append([])
        self.analyze_index.append([])
        self.foreground_threshold.append([])
        self.analysis_params.append(self.default_analysis_params)
        self.Cell_props.append([])
        self.Tissue_props.append({})
        self.remake_side_window()
        self.small_images.append([])

    def SaveCrop(*a):
        ftypes = [('Tiff images', '.tif .tiff')]
        filename = tkinter.filedialog.asksaveasfilename(parent=self.master,
                                                        filetypes=ftypes)
        if filename:
            mask = self.activeCrop
            im_raw = self.im_raw[self.activeImage]
            while round(log(mask.size/im_raw[:, :, 0].size, 2)) < 0:
                mask = skimage.transform.pyramid_expand(mask)
            a = im_raw[:, :, 0][np.ix_(mask.any(1), mask.any(0))]
            a = np.zeros([a.shape[0], a.shape[1], im_raw.shape[2]],
                            dtype=im_raw.dtype)
            for i in range(im_raw.shape[2]):
                a[:, :, i] = im_raw[:, :, i][np.ix_(mask.any(1),
                                                    mask.any(0))]
            im_raw = a
            if not(filename[-4:] == '.tif') | (filename[-5:] == '.tiff'):
                filename += '.tif'
            skimage.io.imsave(filename, im_raw.transpose(2, 0, 1))

    def QuitCrop(*a):
        popup2 = tkinter.Tk()
        popup2.wm_title("Are you sure?")
        label = ttkinter.Label(popup2,
                                text="You are about to quit image cropping")
        label.pack(side="top", fill="x", pady=10)
        B1 = ttkinter.Button(
                popup2, text="Okay",
                command=lambda: [DestroyTK(self.popup), DestroyTK(popup2)])
        B1.pack()
        B2 = ttkinter.Button(popup2, text="Go Back",
                                command=lambda:[DestroyTK(popup2)])
        B2.pack()
        popup2.mainloop()
    im_raw = self.im_raw[self.activeImage]
    image_size = np.size(im_raw)
    if image_size > 100000000:
        im_raw_norm = im_raw.max()
        im_raw = im_raw/im_raw_norm

        while np.size(im_raw) > 100000000:
            im_raw = skimage.transform.pyramid_reduce(im_raw)
        im_raw = im_raw*im_raw_norm
    self.small_images[self.activeImage] = im_raw
    popup = tkinter.Tk()
    self.reduced_image = im_raw
    self.popup = popup
    Color_pointers = self.Color_pointers[self.activeImage]
    color_variable = self.color_variable[self.activeImage]
    Channel_pointers = self.Channel_pointers[self.activeImage]
    channel_variable = self.channel_variable[self.activeImage]
    n_channels = self.n_channels[self.activeImage]
    im_2_display = np.zeros((im_raw.shape[0],
                                im_raw.shape[1], 3), dtype=np.float32)
    for i in range(n_channels):
        Color_pointers_temp = color_variable[i].get()
        Color_pointers[i] = color_variable[i].get()
        Channel_pointers[i] = channel_variable[i].get()
        for j in range(3):
            im_temp = im_2_display[:, :, j]
            im2_add = (im_raw[:, :, i]/im_raw[:, :, i].max())
            im2_add = im2_add * self.LUT[Color_pointers_temp][j]
            im_temp = im_temp + im2_add
            im_2_display[:, :, j] = im_temp
    popup.wm_title("Crop Image")
    label = ttkinter.Label(popup, text="Image crop for file :" +
                            self.FileDictionary[self.activeImage])
    toolbar = tkinter.Frame(popup)
    selectButton = ttkinter.Button(toolbar, text="Select",
                                    command=SelectCrop)
    selectButton.pack(side=tkinter.LEFT, padx=2, pady=2)
    addButton = ttkinter.Button(toolbar, text="Crop",
                                command=addCrop)
    addButton.pack(side=tkinter.LEFT, padx=2, pady=2)
    saveButton = ttkinter.Button(toolbar, text="Save Image",
                                    command=SaveCrop)
    saveButton.pack(side=tkinter.LEFT, padx=2, pady=2)
    cancelButton = ttkinter.Button(toolbar, text="Quit",
                                    command=QuitCrop)
    cancelButton.pack(side=tkinter.LEFT, padx=2, pady=2)
    toolbar.pack(side=tkinter.TOP, fill=tkinter.X)
    label.pack(side="bottom", fill="x", pady=10)
    f = plt.Figure(figsize=(10, 7), dpi=100)
    f.patch.set_visible(False)
    f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    ax = f.gca()
    image_canvas = Tk_Agg.FigureCanvasTkAgg(f, master=popup)
    image_canvas.draw()
    image_canvas.get_tk_widget().pack(side=tkinter.TOP,
                                        expand=True, anchor='n')
    image_toolbar = NT2Tk(image_canvas, popup)
    image_toolbar.update()
    image_canvas._tkcanvas.pack(side=tkinter.BOTTOM, fill=tkinter.BOTH,
                                expand=True)
    ax.clear()
    ax.imshow(im_2_display, aspect='equal')
    self.imCropImage = im_2_display
    ax.autoscale(False)
    ax.axis('off')
    self.ax = ax
    image_canvas.draw()
    self.ax_canvas = image_canvas
    self.popup_statusBar = tkinter.Label(popup, text="Cropping Image",
                                            bd=1, relief=tkinter.SUNKEN,
                                            anchor=tkinter.W)
    self.popup_statusBar.pack(side=tkinter.BOTTOM, fill=tkinter.X)
    if len(self.activeCrop) != 0:
        im_temp = im_2_display
        mask_copy = np.float32(self.activeCrop)
        mask_copy = (1 - mask_copy)/2
        for i in range(3):
            im_temp[:, :, i] = mask_copy + im_temp[:, :, i]
        ax.clear()
        ax.imshow(im_temp, aspect='equal')
        ax.autoscale(False)
        ax.axis('off')
        for i in range(3):
            im_temp[:, :, i] = im_temp[:, :, i] - mask_copy
            self.ax_canvas.draw()
    popup.mainloop()

def ResetCurrent(self):
    def ResetSure(*a):
        n_channels = self.n_channels[self.activeImage]
        self.channel_variable[self.activeImage] = self.channel_variable[
            self.activeImage][0:n_channels]
        self.color_variable[self.activeImage] = self.color_variable[
            self.activeImage][0:n_channels]
        self.Color_pointers[self.activeImage] = self.Color_pointers[
            self.activeImage][0:n_channels]
        self.Channel_pointers[self.activeImage] = self.Channel_pointers[
            self.activeImage][0:n_channels]
        self.im_analyzed[self.activeImage] = []
        self.analyze_index[self.activeImage] = []
        self.foreground_threshold[self.activeImage] = []
        self.analysis_params[self.activeImage] = self.default_analysis_params.copy()
        self.Cell_props[self.activeImage] = []
        self.Tissue_props[self.activeImage] = {}
        self.small_images[self.activeImage] = []
        self.remake_side_window()
    popup2 = tkinter.Tk()
    popup2.wm_title("Are you sure?")
    label = tkinter.Label(
            popup2, text="You are about to reset your current anaylsis")
    label.pack(side="top", fill="x", pady=10)
    B1 = tkinter.Button(popup2, text="Go Ahead",
                        command=lambda: [DestroyTK(popup2), ResetSure()])
    B1.pack()
    B2 = tkinter.Button(popup2, text="Go Back", command=lambda:[DestroyTK(popup2)])
    B2.pack()
    popup2.mainloop()

def RemoveCurrent(self):
    def ResetSure(*a):
        self.n_channels.pop(self.activeImage)
        self.channel_variable.pop(self.activeImage)
        self.color_variable.pop(self.activeImage)
        self.Color_pointers.pop(self.activeImage)
        self.Channel_pointers.pop(self.activeImage)
        self.im_analyzed.pop(self.activeImage)
        self.analyze_index.pop(self.activeImage)
        self.foreground_threshold.pop(self.activeImage)
        self.analysis_params.pop(self.activeImage)
        self.Cell_props.pop(self.activeImage)
        self.Tissue_props.pop(self.activeImage)
        self.small_images.pop(self.activeImage)
        self.im_raw.pop(self.activeImage)
        FileDictionary_init = self.FileDictionary.copy()
        self.FileDictionary = {}
        for i in FileDictionary_init.keys():
            if (i == self.activeImage):
                continue
            self.FileDictionary[len(self.FileDictionary)] = FileDictionary_init[i]
        self.activeImage = 0
        self.remake_side_window()
    popup2 = tkinter.Tk()
    popup2.wm_title("Are you sure?")
    label = tkinter.Label(
            popup2, text="You are about to remove the image and its analysis")
    label.pack(side="top", fill="x", pady=10)
    B1 = tkinter.Button(popup2, text="Go Ahead",
                        command=lambda: [DestroyTK(popup2), ResetSure()])
    B1.pack()
    B2 = tkinter.Button(popup2, text="Go Back", command=lambda:[DestroyTK(popup2)])
    B2.pack()
    popup2.mainloop()

def Reset(self):
    def ResetSure(*a):
        self.im_raw = []
        self.channel_variable = []
        self.color_variable = []
        self.Color_pointers = []
        self.Channel_pointers = []
        self.FileDictionary = {}
        self.ThresholdValues = []
        self.activeImage = []
        self.n_channels = []
        self.activeROI = []
        self.im_analyzed = []
        self.analyze_index = []
        self.foreground_threshold = []
        self.activeFore = []
        self.NucLimits = []
        self.HoleLimits = []
        self.ForeLimits = []
        self.Cell_props = []
        self.Tissue_props = []
        self.overall_data_to_export = {}
        self.analysis_params = []
        self.fill_ch = -1
        self.small_images = []
        self.remake_side_window()
    popup2 = tkinter.Tk()
    popup2.wm_title("Are you sure?")
    label = tkinter.Label(
            popup2, text="You are about to reset your whole anaylsis")
    label.pack(side="top", fill="x", pady=10)
    B1 = tkinter.Button(popup2, text="Go Ahead",
                        command=lambda: [DestroyTK(popup2), ResetSure()])
    B1.pack()
    B2 = tkinter.Button(popup2, text="Go Back", command=lambda:[DestroyTK(popup2)])
    B2.pack()
    popup2.mainloop()