import apec_plot_functions as func

log = True
n = [1, 2, 3, 4, 5, 6]
colors = ['tab:red', 'tab:orange', 'tab:olive', 'tab:green', 'tab:blue', 'tab:purple']
temp = ['0.11 keV', '0.187 keV', '0.318 keV', '0.54 keV', '0.919 keV', '1.56 keV']

data_file_format = "./data_files/28jun25_herve_v2/28jun25bvapecdata_"
model_file_format = "./data_files/28jun25_herve_v2/28jun25bvapecmodel_"

#subplot figures
func.isothermal_model_subplots(data_file_format, n, colors, set = 1, space = 0, temp=temp, log=log, save=True)

func.isothermal_model_subplots(model_file_format, n, colors, set = 1, space = 1, temp=temp, log=log, save=True)

#composite figures
func.isothermal_model_composite(
    file_format=data_file_format, n=n, colors=colors, space = 0, set=1, log=log, temp = temp, save=True)

func.isothermal_model_composite(
    file_format=model_file_format, n=n, colors=colors, space = 1, set=1, log=log, temp = temp, save=True)

#combined spectrum figures
func.isothermal_model_combined(file_format=data_file_format, n=n, space = 0, set=1, log=log, save=True)

func.isothermal_model_combined(file_format=model_file_format, n=n, space = 1, set=1, log=log, save=True)

#stand alone figures
func.isothermal_model_loop(file_format=data_file_format, n=n, colors=colors, set=1, space=0, log = log, temp=temp, save=True)

func.isothermal_model_loop(file_format=model_file_format, n=n, colors=colors, set=1, space=1, log = log, temp=temp, save=True)
