
class WrappedPSFModel(crowdsource.psf.SimplePSF):
    """
    wrapper for photutils GriddedPSFModel
    """
    def __init__(self, psfgridmodel, stampsz=19):
        self.psfgridmodel = psfgridmodel
        self.default_stampsz = stampsz

    def __call__(self, col, row, stampsz=None, deriv=False):

        if stampsz is None:
            stampsz = self.default_stampsz

        parshape = numpy.broadcast(col, row).shape
        tparshape = parshape if len(parshape) > 0 else (1,)

        # numpy uses row, column notation
        rows, cols = np.indices((stampsz, stampsz)) - (np.array([stampsz, stampsz])-1)[:, None, None] / 2.

        # explicitly broadcast
        col = np.atleast_1d(col)
        row = np.atleast_1d(row)
        #rows = rows[:, :, None] + row[None, None, :]
        #cols = cols[:, :, None] + col[None, None, :]

        # photutils seems to use column, row notation
        # only works with photutils <= 1.6.0 - but is wrong there
        #stamps = self.psfgridmodel.evaluate(cols, rows, 1, col, row)
        # it returns something in (nstamps, row, col) shape
        # pretty sure that ought to be (col, row, nstamps) for crowdsource

        # andrew saydjari's version here:
        # it returns something in (nstamps, row, col) shape
        stamps = []
        for i in range(len(col)):
            # the +0.5 is required to actually center the PSF (empirically)
            #stamps.append(self.psfgridmodel.evaluate(cols+col[i]+0.5, rows+row[i]+0.5, 1, col[i], row[i]))
            # the above may have been true when we were using (incorrectly) offset PSFs
            stamps.append(self.psfgridmodel.evaluate(cols+col[i], rows+row[i], 1, col[i], row[i]))

        stamps = np.array(stamps)

        # for oversampled stamps, they may not be normalized
        stamps /= stamps.sum(axis=(1,2))[:,None,None]
        # this is evidently an incorrect transpose
        #stamps = np.transpose(stamps, axes=(0,2,1))

        if deriv:
            dpsfdrow, dpsfdcol = np.gradient(stamps, axis=(1, 2))

        ret = stamps
        if parshape != tparshape:
            ret = ret.reshape(stampsz, stampsz)
            if deriv:
                dpsfdrow = dpsfdrow.reshape(stampsz, stampsz)
                dpsfdcol = dpsfdcol.reshape(stampsz, stampsz)
        if deriv:
            ret = (ret, dpsfdcol, dpsfdrow)

        return ret

    def render_model(self, col, row, stampsz=None):
        """
        this function likely does nothing?
        """
        if stampsz is not None:
            self.stampsz = stampsz

        rows, cols = np.indices(self.stampsz, dtype=float) - (np.array(self.stampsz)-1)[:, None, None] / 2.

        return self.psfgridmodel.evaluate(cols, rows, 1, col, row).T.squeeze()




def inject_synthetic_stars(img):
    # use psf_grid in stpsf package to inject synthetic stars
    pass


def get_psf_model(filtername, proposal_id, field,
                  module,
                  obsdate=None,
                  target='w51',
                  stampsz=19,
                  oversample=1,
                  basepath='/orange/adamginsburg/jwst/'):
    """
    Return two types of PSF model, the first for DAOPhot and the second for Crowdsource
    """
    if filtername.upper() in ['F140M', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                              'F200W', 'F210M', 'F212N', 'F250M', 'F300M', 'F322W2',
                              'F335M', 'F356W', 'F360M', 'F410M', 'F430M', 'F444W',
                              'F460M', 'F466N', 'F480M']:
        instrument = 'NIRCam'
    else:
        instrument = 'MIRI'

    basepath = f'{basepath}/{target}'
    psf_dir = "."

    psf_fn = (
        f"{psf_dir}/"
        f"{instrument.lower()}_{module}_{filtername.lower()}_fovp101_samp4_npsf16.fits"
    )
    print('psf_fn', psf_fn, flush=True)
    if os.path.exists(psf_fn):
        print(f"Loading cached PSF grid: {psf_fn}", flush=True)
        grid = to_griddedpsfmodel(psf_fn)
        psf_model = WrappedPSFModel(grid, stampsz=stampsz)
        
    else:
        
        has_downloaded = False
        ntries = 0
        
        while not has_downloaded:
            ntries += 1
            try:
                print("Attempting to download WebbPSF data", flush=True)
                if filtername.upper() in ['F140M', 'F150W', 'F162M', 'F164N', 'F182M', 'F187N',
                                  'F200W', 'F210M', 'F212N', 'F250M', 'F300M', 'F322W2',
                                  'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W',
                                  'F460M', 'F466N', 'F480M']:
                    nrc = webbpsf.NIRCam()
                    #print('nrc=webbpsf.NIRCam()', flush=True)
                else:
                    nrc = webbpsf.MIRI()
                   # print('nrc=webbpsf.MIRI()', flush=True)
                nrc.load_wss_opd_by_date(f'{obsdate}T00:00:00')
                nrc.filter = filtername
                if module in ('nrca', 'nrcb'):
                    if 'F4' in filtername.upper() or 'F3' in filtername.upper():
                        nrc.detector = f'{module.upper()}5' # I think NRCA5 must be the "long" detector?
                    else:
                        nrc.detector = f'{module.upper()}1' #TODO: figure out a way to use all 4?
                    # default oversampling is 4
                    grid = nrc.psf_grid(num_psfs=16, all_detectors=False, verbose=True, save=True)
                elif 'mirimage' in module:
                    print('module', module, flush=True)
                    print(nrc.detector)
                    nrc.detector = 'MIRIM'
                    grid = nrc.psf_grid(num_psfs=16, all_detectors=False, verbose=True, save=True)
                else:
                    grid = nrc.psf_grid(num_psfs=16, all_detectors=True, verbose=True, save=True)
                has_downloaded = True
            except (urllib3.exceptions.ReadTimeoutError, requests.exceptions.ReadTimeout, requests.HTTPError) as ex:
                print(f"Failed to build PSF: {ex}", flush=True)
            except Exception as ex:
                print(ex, flush=True)
                if ntries > 10:
                    # avoid infinite loops
                    raise ValueError("Failed to download PSF, probably because of an error listed above")
                else:
                    continue
        if isinstance(grid, list):
            grid = grid[0]

        #yy, xx = np.indices([31,31], dtype=float)
        #grid.x_0 = grid.y_0 = 15.5
        #psf_model = crowdsource.psf.SimplePSF(stamp=grid(xx,yy))

        # bigger PSF probably needed
        yy, xx = np.indices([61, 61], dtype=float)
        grid.x_0 = grid.y_0 = 30
        psf_model = crowdsource.psf.SimplePSF(stamp=grid(xx, yy))

            
        
    return grid, psf_model
        
def load_data(filename):
    fh = fits.open(filename)
    im1 = fh
    data = im1['SCI'].data
    try:
        wht = im1['WHT'].data
    except KeyError:
        wht = None
    err = im1['ERR'].data
    instrument = im1[0].header['INSTRUME']
    telescope = im1[0].header['TELESCOP']
    obsdate = im1[0].header['DATE-OBS']
    return fh, im1, data, wht, err, instrument, telescope, obsdate    

def get_filenames(basepath, filtername, proposal_id, field, each_suffix, module, pupil='clear', visitid='001'):

    # jw01182004002_02101_00012_nrcalong_destreak_o004_crf.fits
    # jw02221001001_07101_00012_nrcalong_destreak_o001_crf.fits
    # jw02221001001_05101_00022_nrcb3_destreak_o001_crf.fits
        #jw06151002001_02101_00001_mirimage_i2d.fits

    glstr = f'{basepath}/{filtername}/pipeline/jw0{proposal_id}{field}*{module}*_{each_suffix}.fits'
    
  
    fglob = glob.glob(glstr)
    for st in fglob:
        print(st)
        if 'align' in st or 'uncal' in st:
            print(f"Removing {st} from glob string because it is an alignment file")
            fglob.remove(st)
    if len(fglob) == 0:
        raise ValueError(f"No matches found to {glstr}")
    else:
        return fglob


def inject_synthetic_stars(img, dqarr, psf_model, flux_range, ww,  filtername,fwhm_pix, grid_spacing=50, err=None, stampsz=19, nsigma=4):
    """
    Inject synthetic stars in a grid pattern across the image
    
    Parameters
    ----------
    img : ndarray
        The image data to inject stars into
    psf_model : WrappedPSFModel
        The PSF model to use for injection
    flux_range : array-like
        Array of flux values to test
    grid_spacing : int
        Spacing between injected stars in pixels
    err : ndarray, optional
        Error array for the image
    stampsz : int
        Size of PSF stamp
        
    Returns
    -------
    results : dict
        Dictionary containing injection and recovery information for each flux level
    """
    from photutils.detection import DAOStarFinder
    from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve_fft
    
    pixel_scale = ww.proj_plane_pixel_area()**0.5
    ny, nx = img.shape
    
    # Create grid of injection positions (avoiding edges)
    margin = stampsz
    x_positions = np.arange(margin, nx - margin, grid_spacing)
    y_positions = np.arange(margin, ny - margin, grid_spacing)
    xx, yy = np.meshgrid(x_positions, y_positions)
    
    # Flatten to get list of positions
    inject_x = xx.flatten()
    inject_y = yy.flatten()

    # create a 3d array of shape (len(x_positions), len(y_positions), len(flux_range)) to check which flux level the injected stars are detected at
    completeness_arr = np.full((len(y_positions), len(x_positions)), np.nan)

    print(f"Injecting {len(inject_x)} stars per flux level")
    edge_count =0
    
    completeness_values = []
    for iiflux, flux in enumerate(flux_range):
        print(f"\nTesting flux level: {flux:.2e}")
        
        # Create image with injected stars
        img_with_stars = img.copy()

        
        # Inject stars at grid positions (vectorized)
        psf_stamps = np.asarray(psf_model(inject_x, inject_y, stampsz=stampsz))

        # Force shape to (nstar, stamp_h, stamp_w)
        if psf_stamps.ndim == 2:
            psf_stamps = psf_stamps[None, :, :]
        elif psf_stamps.ndim == 4 and psf_stamps.shape[1] == 1:
            psf_stamps = psf_stamps[:, 0, :, :]
        elif psf_stamps.ndim != 3:
            raise ValueError(f"Unexpected PSF stamp shape: {psf_stamps.shape}")

        # Scale all stamps once per flux level
        scaled_flux = (flux / (pixel_scale)**2).to(u.MJy / u.sr)
        scaled_stamps = psf_stamps * scaled_flux.value

        # Use actual stamp shape (not assumed stampsz) so indices and values match
        nstar, stamp_h, stamp_w = scaled_stamps.shape
        y_lo = stamp_h // 2
        y_hi = stamp_h - y_lo
        x_lo = stamp_w // 2
        x_hi = stamp_w - x_lo
        # add saturated or nan pixels into the mask
        is_nan = np.isnan(img)
        is_saturated = (dqarr & dqflags.pixel['SATURATED']) != 0

        # Index into 2D masks at injection positions
        is_nan_at_pos = is_nan[inject_y.astype(int), inject_x.astype(int)]
        is_saturated_at_pos = is_saturated[inject_y.astype(int), inject_x.astype(int)]

        valid_mask = (
            (inject_y >= y_lo) & (inject_y < ny - (y_hi - 1)) &
            (inject_x >= x_lo) & (inject_x < nx - (x_hi - 1)) &
            (~is_saturated_at_pos) & (~is_nan_at_pos)
)

        edge_count = int(np.sum(~valid_mask))
        valid_x = inject_x[valid_mask].astype(int)
        valid_y = inject_y[valid_mask].astype(int)
        valid_stamps = scaled_stamps[valid_mask]

        
        if valid_stamps.size > 0:
            y_offsets = np.arange(-y_lo, y_hi, dtype=int)
            x_offsets = np.arange(-x_lo, x_hi, dtype=int)

            yy = valid_y[:, None, None] + y_offsets[None, :, None]
            xx = valid_x[:, None, None] + x_offsets[None, None, :]

            bshape = np.broadcast_shapes(yy.shape, xx.shape)
            n_index = int(np.prod(bshape))
            n_value = int(valid_stamps.size)

            if n_index != n_value:
                raise ValueError(
                    f"Index/value mismatch after broadcast: idx={n_index}, vals={n_value}, "
                    f"yy_shape={yy.shape}, xx_shape={xx.shape}, stamp_shape={valid_stamps.shape}"
                )

            np.add.at(img_with_stars, (yy, xx), valid_stamps)

            filtered_errest = np.nanmedian(err)
            daofind_tuned = DAOStarFinder(
                threshold=nsigma * filtered_errest,
                fwhm=fwhm_pix, roundhi=1.0, roundlo=-1.0,
                sharplo=0.30, sharphi=1.40
            )

            mask = np.isnan(img_with_stars)
            is_saturated = (dqarr & dqflags.pixel['SATURATED']) != 0
            is_saturated = ndimage.binary_dilation(is_saturated, iterations=1)
            img_ = img_with_stars.copy()
            img_[is_saturated] = np.nan
            mask |= is_saturated

            found_stars = daofind_tuned(img_, mask=mask)

        # Match injected to detected stars using ONLY valid injected positions
        n_injected = len(valid_x)

        if n_injected == 0:
            n_recovered = 0
            completeness = np.nan
            print("No valid injection positions at this flux level")
        elif (found_stars is not None) and (len(found_stars) > 0):
            tolerance = pixel_scale  # same angular unit family as sep from SkyCoord

            skycoord_found_stars = ww.pixel_to_world(found_stars['xcentroid'], found_stars['ycentroid'])
            skycoord_injected_valid = ww.pixel_to_world(valid_x, valid_y)

            nearest_inj_idx, sep, _ = skycoord_found_stars.match_to_catalog_sky(
                skycoord_injected_valid, nthneighbor=1
            )

            matched_det = sep < tolerance

            # Count unique injected stars recovered (prevents multi-detection overcount)
            recovered_valid_idx = np.unique(nearest_inj_idx[matched_det])
            n_recovered = recovered_valid_idx.size
            completeness = n_recovered / n_injected

            # Fill completeness map by injected-grid cell, not detected centroid cell
            valid_full_idx = np.flatnonzero(valid_mask)
            recovered_full_idx = valid_full_idx[recovered_valid_idx]
            grid_y = recovered_full_idx // len(x_positions)
            grid_x = recovered_full_idx % len(x_positions)
            completeness_arr[grid_y, grid_x] = flux

            print(f"Recovered {n_recovered} out of {n_injected} injected stars")
        else:
            n_recovered = 0
            completeness = 0.0
            print(f"Recovered {n_recovered} out of {n_injected} injected stars")
                
                # ...existing code...
        completeness_values.append(completeness)
    # convert check array to 2d array of shape (len(y_positions), len(x_positions)) where each element is the highest flux level at which the star was undetected


    # make an array of the same shape as completeness_arr where each element is the average pixel values within aperture at the position of the injected stars
    avg_flux_app_array = np.full((len(y_positions), len(x_positions)), np.nan)
    aperture_radius = 2 * fwhm_pix  # use FWHM as aperture radius
    y_indices, x_indices = np.indices(img.shape)
    for i in range(len(y_positions)):
        for j in range(len(x_positions)):
            x_center = x_positions[j]
            y_center = y_positions[i]
            r = np.sqrt((x_indices - x_center)**2 + (y_indices - y_center)**2)
            aperture_mask = r <= aperture_radius
            avg_flux_app_array[i, j] = np.nanmean(img[aperture_mask])

    return completeness_arr, completeness_values, avg_flux_app_array

def plot_global_completeness(completeness_values, flux_range, filename, savedir='/orange/adamginsburg/w51/jwst/plots/completeness_test/', th=0.8):
    import matplotlib.pyplot as plt
    from scipy.interpolate import PchipInterpolator

    # detection fraction is the fraction of the finite values in completeness_arr that are less than or equal to each flux level in flux_range
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.scatter(flux_range, completeness_values, color='blue', label='Completeness')
    pchip_func = PchipInterpolator(flux_range, completeness_values)
    x_smooth = np.logspace(np.log10(flux_range[0].value), np.log10(flux_range[-1].value), 100)
    ax.plot(x_smooth, pchip_func(x_smooth), color='red', ls='dotted')

    # plot horiontal line where dtection fraction is 0.9
    ax.axhline(th, color='gray', ls='dashed')
    # get the flux level where the detection fraction is 0.9
    flux_th = x_smooth[np.argmin(np.abs(pchip_func(x_smooth) - th))]
    print(f"{th*100}% completeness at flux level: {flux_th:.2e} Jy")
    ax.set_xscale('log')
    ax.set_xlabel('Flux (Jy)')
    ax.set_ylabel('Detection Fraction')
    ax.set_title(f'{filename.split("/")[-1]}')

    plt.savefig(f'{savedir}/{filename.split("/")[-1].replace(".fits", "_completeness")}.png')

    return flux_th


def main( stampsz=19, ):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--filtername", dest="filtername",
                    default='F140M',
                    help="filter name", metavar="filtername")
    parser.add_option("-m", "--modules", dest="modules",
                    #default='nrca,nrcb,merged,merged-reproject',
                    default='nrca,nrcb,merged',
                    help="module list", metavar="modules") 
    parser.add_option("--proposal_id", dest="proposal_id",
                    default='6151',
                    help="proposal_id", metavar="proposal_id")
    parser.add_option("--target", dest="target",
                    default='w51',
                    help="target", metavar="target")
    (options, args) = parser.parse_args()
    modules = options.modules.split(',')
    filtername = options.filtername
    proposal_id = options.proposal_id
    target = options.target
    nvisits = {'2221': {'brick': 1, 'cloudc': 2},
               '1182': {'brick': 2},
               '6151': {'w51': 1, 'w51_miri': 2}
               }
    field_to_reg_mapping = {'2221': {'001': 'brick', '002': 'cloudc'},
                            '1182': {'004': 'brick'},
                            '6151': {'001': 'w51', '002':'w51_miri'}}[proposal_id]
    reg_to_field_mapping = {v:k for k,v in field_to_reg_mapping.items()}
    field = reg_to_field_mapping[target]
    target_dir = target
    if proposal_id == '6151' and target == 'w51_miri':
        target_dir = 'w51'
    basepath = f'/orange/adamginsburg/jwst/{target_dir}/'



    # need to have incrementing _before_ test
    index = -1
    completeness_arr = []
    for module in modules:
        detector = module # no sub-detectors for long-NIRCAM
       
          
        for ii, visitid in enumerate(range(1, nvisits[proposal_id][target] + 1)):
               

            visitid = f'{visitid}'

            filenames = get_filenames(basepath, filtername, proposal_id,
                                                    field, visitid=visitid,
                                                    each_suffix='cal',
                                                    module=module, pupil='clear') 
            print('filenames', filenames, flush=True)
            for filename in filenames:
                index += 1
                # enable array jobs
                if os.getenv('SLURM_ARRAY_TASK_ID') is not None and int(os.getenv('SLURM_ARRAY_TASK_ID')) != index:
                    print(f'Task={os.getenv("SLURM_ARRAY_TASK_ID")} does not match index {index}')
                    continue
                fh, im1, data, wht, err, instrument, telescope, obsdate = load_data(filename)
                img = im1['SCI'].data
                # set up coordinate system
                ww = wcs.WCS(im1[1].header)
                pixscale = ww.proj_plane_pixel_area()**0.5
                cen = ww.pixel_to_world(im1[1].shape[1]/2, im1[1].shape[0]/2)
                
                grid, psf_model = get_psf_model(filtername, proposal_id, field,
                                                module=module,
                                                # if we're doing each exposure, we want the full grid
                                                target=target,
                                                obsdate=obsdate,
                                                basepath='/orange/adamginsburg/jwst/')
                
                # add grid on the bkg_img
            

                flux_range = np.logspace(-9, -2, 40) * u.Jy

                
                fwhm_tbl = Table.read(f'{basepath}/reduction/fwhm_table.ecsv')
                row = fwhm_tbl[fwhm_tbl['Filter'] == filtername]
                fwhm = fwhm_arcsec = float(row['PSF FWHM (arcsec)'][0])
                fwhm_pix = float(row['PSF FWHM (pixel)'][0])
                kernel = Gaussian2DKernel(x_stddev=fwhm_pix/2.355)
                mask = np.isnan(img)
                if 'DQ' in im1:
                    dqarr = im1['DQ'].data
                    is_saturated = (dqarr & dqflags.pixel['SATURATED']) != 0
                    # we want original image_with_synthetic_stars to be untouched for imshowing diagnostics etc.
                    img_ = img.copy()
                    img_[is_saturated] = np.nan
                    mask |= is_saturated
                else:
                    img_ = img

                nan_replaced_data = interpolate_replace_nans(img_, kernel, convolve=convolve_fft)

                completeness_map, completeness_values, avg_flux_app_array = inject_synthetic_stars(img_, dqarr, psf_model, flux_range, ww, filtername, fwhm_pix, grid_spacing=50, err=err, stampsz=stampsz)
                savedir = f'/orange/adamginsburg/w51/jwst/plots/completeness_test/{filtername}'
                if not os.path.exists(savedir):
                    os.makedirs(savedir)
                flux_80 = plot_global_completeness(completeness_values, flux_range, filename, savedir=savedir, th=0.8)

                # save flux_80 to a text file
                with open(f'{savedir}/{filename.split("/")[-1].replace(".fits", "_completeness.txt")}', 'w') as f:
                    f.write(f"{flux_80:.2e} Jy\n")    

if __name__ == "__main__":
    main()