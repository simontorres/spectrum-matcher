from ccdproc import CCDData
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import (Slider, Button, SpanSelector)
from pipeline.wcs import WCS
import re
import os
import numpy as np
from astropy.modeling import models

wcs = WCS()

matplotlib.rcParams['toolbar'] = 'None'

class MatchLamps(object):

    def __init__(self, reference_dir):
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.ax5 = None
        self.ax6 = None
        self.ax7 = None
        self.__zero_xlim = None
        self.__zero_ylim = None
        self.slider_1 = None
        self.reference_wavelength = None
        self.stretch = None
        self.ref_dir = reference_dir
        self.ref_wav = None
        self.ref_intens = None
        self.new_lamp = None
        self.new_ccd = None
        self.ref_plot = None
        self.ref_ccd = None
        self.reference_pixel = None
        self.ref_pix_line = None
        self.magnification = None
        self.second_order_par = None
        self.model = None

    def __call__(self, reference, new_file):
        n_columns = 15
        n_rows = 40
        self.reference_file = reference
        self.ref_ccd = CCDData.read(os.path.join(self.ref_dir,
                                                 self.reference_file),
                                    unit='adu')
        self.new_ccd = CCDData.read(new_file, unit='adu')
        self.ref_wav, self.ref_intens = wcs.read(self.ref_ccd)
        self.fig = plt.figure()
        self.ax1 = plt.subplot2grid((n_rows, n_columns), (0, 0), rowspan=27, colspan=n_columns)
        plt.setp(self.ax1.get_xticklabels(), visible=False)
        self.ax2 = plt.subplot2grid((n_rows, n_columns), (27, 0), rowspan=4, colspan=n_columns, sharex=self.ax1)
        self.ax6 = plt.subplot2grid((n_rows, n_columns), (n_rows - 7, 0), colspan=n_columns)
        self.ax5 = plt.subplot2grid((n_rows, n_columns), (n_rows - 6, 0), colspan=n_columns)
        self.ax3 = plt.subplot2grid((n_rows, n_columns), (n_rows - 5, 0), colspan=n_columns)
        self.ax4 = plt.subplot2grid((n_rows, n_columns), (n_rows - 4, 0), colspan=n_columns)
        self.ax7 = plt.subplot2grid((n_rows, n_columns), (n_rows - 3, 0), colspan=n_columns)
        self.ax10 = plt.subplot2grid((n_rows, n_columns), (n_rows - 2, 0), colspan=n_columns)
        self.ax8 = plt.subplot2grid((n_rows, n_columns), (n_rows - 1, n_columns - 2), colspan=2)
        self.ax9 = plt.subplot2grid((n_rows, n_columns), (n_rows - 1, n_columns - 6), colspan=2)
        self.ax11 = plt.subplot2grid((n_rows, n_columns), (n_rows - 1, n_columns - 4), colspan=2)
        self._plot_reference()
        self.reference_wavelength = Slider(self.ax3, 'Reference Wavelength',
                                           self.ref_wav[0],
                                           self.ref_wav[-1],
                                           valinit=self.ref_wav[0])
        self.reference_wavelength.on_changed(self._update_plots)
        self.stretch = Slider(self.ax4, 'Linear Parameter', 0, 2, valinit=self.ref_wav[1] - self.ref_wav[0])
        self.stretch.on_changed(self._update_plots)
        self.second_order_par = Slider(self.ax7, 'Second Order', -1e-4, 1e-4, valinit=0, valfmt="%.2E")
        self.second_order_par.on_changed(self._update_plots)
        self.reference_pixel = Slider(self.ax5, 'Reference Pixel', 0, 4060, valinit=0)
        self.reference_pixel.on_changed(self._update_plots)
        self.magnification = Slider(self.ax6, 'Magnification', 1, 100, valinit=1)
        self.magnification.on_changed(self._update_plots)
        self.third_order = Slider(self.ax10, 'Third Order', -1e-8, 1e-8, valinit=0, valfmt="%.2E")
        self.third_order.on_changed(self._update_plots)


        self.reference_range = SpanSelector(self.ax1, self._plot_reference, 'horizontal', rectprops=dict(alpha=0.5, facecolor='red'))
        self.reset_button = Button(self.ax8, 'Reset Values')
        self.reset_button.on_clicked(self.__reset_values)
        self.reset_zoom = Button(self.ax11, 'Reset Zoom')
        self.reset_zoom.on_clicked(self.__reset_zoom)

        self.save_button = Button(self.ax9, "Save Solution")
        self.save_button.on_clicked(self.__save_solution)
        self._update_plots(value=0)
        plt.subplots_adjust(left=0.125,
                            right=0.930,
                            top=0.985,
                            bottom=0.015,
                            hspace=0.17,
                            wspace=0.11)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        plt.show()

    def __save_solution(self, event):
        self.new_ccd = wcs.write_gsp_wcs(self.new_ccd, self.model)
        self.new_ccd.write('new_file.fits', overwrite=True)

    def __reset_zoom(self, event):
        if self.__zero_xlim is not None and self.__zero_ylim is not None:
            self.ax1.set_xlim(self.__zero_xlim)
            self.ax1.set_ylim(self.__zero_ylim)

    def __reset_values(self, event):
        self.reference_pixel.reset()
        self.reference_wavelength.reset()
        self.magnification.reset()
        self.stretch.reset()
        self.second_order_par.reset()
        self.third_order.reset()



    def _update_range(self, xmin, xmax):
        index_min, index_max = np.searchsorted(self.ref_wav, (xmin, xmax))
        index_max = min(len(self.ref_wav) - 1, index_max)

        this_wav = self.ref_wav[index_min:index_max]
        this_inte = self.ref_intens[index_min:index_max]
        self.ref_plot.set_data(this_wav, this_inte)
        self.ax1.set_xlim(this_wav[0], this_wav[-1])
        self.ax1.set_ylim(this_inte.min(), this_inte.max())
        self.fig.canvas.draw()

    def _plot_residuals(self):
        pass


    def _plot_reference(self, xmin=None, xmax=None):
        self.ref_plot, = self.ax1.plot(self.ref_wav, self.ref_intens, color='b')
        if xmin is not None and xmax is not None:
            index_min, index_max = np.searchsorted(self.ref_wav, (xmin, xmax))
            index_max = min(len(self.ref_wav) - 1, index_max)
            this_wav = self.ref_wav[index_min:index_max]
            this_inte = self.ref_intens[index_min:index_max]
            self.ref_plot.set_data(this_wav, this_inte)
            self.ax1.set_xlim(this_wav[0], this_wav[-1])
            self.ax1.set_ylim(this_inte.min(), this_inte.max())

        elif self.__zero_xlim is None and self.__zero_ylim is None:
            self.__zero_xlim = (self.ref_wav[0], self.ref_wav[-1])
            self.__zero_ylim = (self.ref_intens.min(), self.ref_intens.max())
            self.ax1.set_xlim(self.__zero_xlim)
            self.ax1.set_ylim(self.__zero_ylim)
        # else:
        #     self.ax1.set_xlim(self.__zero_xlim)
        #     self.ax1.set_ylim(self.__zero_ylim)

            # self.ax1.set_xlim(xmin, xmax)
            #
            # self.ax1.set_ylim()

    def _update_plots(self, value):
        try:
            self.new_lamp.remove()
            self.ax1.relim()
        except ValueError:
            pass
        except AttributeError:
            pass

        self.model = models.Chebyshev1D(degree=3)
        self.model.c0.value = self.reference_wavelength.val - (self.reference_pixel.val * self.stretch.val)
        self.model.c1.value = self.stretch.val
        self.model.c2.value = self.second_order_par.val
        self.model.c3.value = self.third_order.val

        try:
            self.ref_plot.remove()
            self.ax1.relim()
            # self._plot_reference(xmin=cheb(0), xmax=cheb(len(self.new_ccd.data)))
            self._plot_reference()
        except ValueError:
            pass
        except AttributeError:
            pass

        self.new_lamp, = self.ax1.plot(self.model(range(len(self.new_ccd.data))),
                                       self.new_ccd.data * self.magnification.val, color='r')

        try:
            self.ref_pix_line.remove()
            self.ax1.relim()
        except ValueError:
            pass
        except AttributeError:
            pass

        self.ref_pix_line = self.ax1.axvline(self.model(self.reference_pixel.val), color='k')


ref_dir = "/data/simon/documentation/soar/general_documentation/" \
                       "WCS_problem/noao_ref"

ref_file = 'cuar-3000-7050.fits'

path_william = '/user/simon/data/soar/work/william_data/2016-03-07/comp_files'
filw = os.path.join(
            path_william,
            "ll_ext_cfzsto_0141.SO2016A-006-0141.SO2016A-006_0307.fits")

match_lamps = MatchLamps(reference_dir=ref_dir)

match_lamps(reference=ref_file, new_file=filw)
#
# path_ref = "/data/simon/documentation/soar/general_documentation/WCS_problem/noao_ref"
#
# path_lib = "/data/simon/data/soar/comp_lamp_lib/work/comparison_lamp_library/not_completed/"
# path_lib_completed = "/data/simon/data/soar/comp_lamp_lib/work/comparison_lamp_library/completed/"
#
# path_william = '/user/simon/data/soar/work/william_data/2016-03-07/comp_files'
# # file_1 = os.path.join(path_ref,"fear-3000-7050.fits")
# file_1 = os.path.join(path_ref,"cuar-3000-7050.fits")
# file_2 = os.path.join(path_william,
#                       "ll_ext_cfzsto_0141.SO2016A-006-0141.SO2016A-006_0307.fits")
#
#
# file_3 = os.path.join(path_lib_completed,
#                       "ll_ext_cfzsto_0152-0156_goodman_comp_600Blue_CuHeAr.fits")
# ccd1 = CCDData.read(file_1, unit=u.adu)
# wav_int_1 = wcs.read(ccd=ccd1)
# ccd2 = CCDData.read(file_2, unit=u.adu)
# ccd3 = CCDData.read(file_3, unit=u.adu)
#
# fig, (ax1, ax2) = plt.subplots(2, 1)
#
#
#
# #
# ax2.set_title(os.path.basename(file_2))
# ax2.plot(ccd2.data)
# ax2.set_xlim((0, ccd2.data.size))
#
# if True:
#     ax1.set_title(os.path.basename(file_3))
#     ax1.plot(ccd3.data)
#     pix_keys = ccd3.header['GSP_P*']
#     for key in pix_keys:
#         if re.match(r'GSP_P\d{3}', key) is not None:
#             x_val = pix_keys[key]
#             y_val = ccd3.data[int(round(x_val))]
#             text = ccd3.header[re.sub('GSP_P', 'GSP_A', key)]
#             ax1.text(x_val, y_val, text, rotation=270)
#     ax1.set_xlim((0, ccd3.data.size))
# else:
#     ax1.set_title('Reference from NOAO')
#     ax1.plot(wav_int_1[0], wav_int_1[1])
#     ax1.set_xlim((3500, 6160))
#
# pix_keys = ccd2.header['GSP_P*']
# for key in pix_keys:
#     if re.match(r'GSP_P\d{3}', key) is not None:
#         x_val = pix_keys[key]
#         y_val = ccd2.data[int(round(x_val))]
#         ax2.text(x_val, y_val, x_val, rotation=270)
#         # print(key, pix_keys[key])
# plt.show()