#! /usr/bin/env python
import ConfigParser


Config = ConfigParser.ConfigParser()

# create config file
cfg = open('SE_parameters.cfg', 'w')

# add section for calibration
Config.add_section('Calibration')
Config.set('Calibration', '-c', 'config.sex')
Config.set('Calibration', '-THRESH_TYPE', 'RELATIVE')
Config.set('Calibration', '-DETECT_MINAREA', '5')
Config.set('Calibration', '-DETECT_THRESH', '1.0')
Config.set('Calibration', '-ANALYSIS_THRESH', '1.0')
Config.set('Calibration', '-WEIGHT_TYPE', 'NONE')
Config.set('Calibration', '-PHOT_APERTURES', '30')
Config.set('Calibration', '-BACK_SIZE', '64')
Config.set('Calibration', '-MAG_ZEROPOINT', '0.0')

# section for final catalog
Config.add_section('Catalog')
Config.set('Catalog', '-c', 'config.sex')
Config.set('Catalog', '-THRESH_TYPE', 'RELATIVE')
Config.set('Catalog', '-DETECT_MINAREA', '5')
Config.set('Catalog', '-DETECT_THRESH', '2.2')
Config.set('Catalog', '-ANALYSIS_THRESH', '2.2')
Config.set('Catalog', '-WEIGHT_TYPE', 'NONE')
Config.set('Catalog', '-PHOT_APERTURES', '30')
Config.set('Catalog', '-BACK_SIZE', '64')
Config.set('Catalog', '-MAG_ZEROPOINT', '0.0')

Config.write(cfg)
cfg.close()

