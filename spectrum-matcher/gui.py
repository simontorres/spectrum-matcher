from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
from matplotlib.widgets import (Slider, SpanSelector)

class MatchGui(object):

    def __init__(self, rows=40, columns=1, subplots=7):
        self.fig = plt.figure()

        self.fig.canvas.set_window_title('Spectrum Matcher')

        data_width = int(0.5 * rows)

        residual_width = int(0.3 * data_width)

        self.ax_data = plt.subplot2grid((rows, columns),
                                        (0, 0),
                                        rowspan=data_width)

        self.ax_residuals = plt.subplot2grid((rows, columns),
                                             (data_width + 1, 0),
                                             rowspan=residual_width)

        self.ax_slider_magnification = plt.subplot2grid((rows, columns),
                                                        (rows - 5, 0))

        self.magnification = Slider(self.ax_slider_magnification,
                                    'Magnification',
                                    1, 100,
                                    valinit=1)

        self.magnification.on_changed(self._update_plots)

        self.ax_slider_ref_pix = plt.subplot2grid((rows, columns),
                                                  (rows - 4, 0))

        self.reference_pixel = Slider(self.ax_slider_ref_pix,
                                      'Reference Pixel',
                                      0,
                                      4600,
                                      valinit=2300)

        self.reference_pixel.on_changed(self._update_plots)

        self.ax_slider_ref_wave = plt.subplot2grid((rows, columns),
                                                   (rows - 3, 0))

        self.reference_wave = Slider(self.ax_slider_ref_wave,
                                     'Reference Wavelength',
                                     3000,
                                     10000,
                                     valinit=5000)

        self.reference_wave.on_changed(self._update_plots)

        self.ax_slider_linear_par = plt.subplot2grid((rows, columns),
                                                     (rows - 2, 0))
        self.linear_par = Slider(self.ax_slider_linear_par,
                                 'Linear Parameter',
                                 0,
                                 2,
                                 valinit=0.65)

        self.linear_par.on_changed(self._update_plots)

        self.ax_slider_second_order = plt.subplot2grid((rows, columns),
                                                       (rows - 1, 0))

        self.second_order = Slider(self.ax_slider_second_order,
                                   'Second Order',
                                   0,
                                   1e-4,
                                   valinit=1e-7)

        self.second_order.on_changed(self._update_plots)

    def __call__(self, *args, **kwargs):
        manager = plt.get_current_fig_manager()
        if plt.get_backend() == u'GTK3Agg':
            manager.window.maximize()
        elif plt.get_backend() == u'Qt5Agg':
            manager.window.showMaximized()
        plt.show()

    def _update_plots(self, value):
        pass


if __name__ == '__main__':
    match_gui = MatchGui()
    match_gui()
