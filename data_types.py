from collections import namedtuple

Range = namedtuple('Range', ['start', 'end'])

ColumnIndexes = namedtuple('ColumnIndexes', ['wavelength', 'flux', 'error'])

PointData = namedtuple('PointData', ['wavelength', 'flux', 'error'])

RangesData = namedtuple('RangesData', ['wavelength', 'flux', 'error'])

FigureData = namedtuple('FigureData', ['spectrum_file_name', 'wavelength_from', 'wavelength_to', 'z', 'snr', 'snr_mean_in_ehvo'])

FigureDataOriginal = namedtuple('FigureDataOriginal', ['FigureData','bf', 'cf', 'power_law_data_x', 'power_law_data_y'])

DataNormalized = namedtuple('DataNormalized', ['flux_normalized', 'error_normalized'])