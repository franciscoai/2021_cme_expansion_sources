def save_img(path,name)
    import matplotlib.pyplot as plt
    import astropy.units as u
    import sunpy.data.sample
    import sunpy.map
    from astropy.io import fits
    hmi_map= sunpy.map.Map(path)
    hmi_map.plot_settings['cmap'] = "hmimag"
    hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
    fig = plt.figure(figsize=(6, 5))
    ax2 = fig.add_subplot(1, 1, 1, projection=hmi_map)
    hmi_map.plot(axes=ax2)
    plt.savefig(name + '.png')
    return



