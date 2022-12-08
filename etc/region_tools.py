from .io_helper import ofn_support

@ofn_support
def source_region(source_name, observation=None, radius=15):
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    import regions
    import astropy.units as u

    customSimbad = Simbad()
    customSimbad.get_votable_fields()
    result_table = customSimbad.query_object(source_name)
    #print(len(result_table))
    if len(result_table) != 1:
        raise Exception("More than one entry for "+str(source_name)+" in Simbad")
    result_table = result_table[0]
    center_sky = SkyCoord(result_table["RA"], result_table["DEC"], unit=(u.hourangle, u.deg), frame='fk5')

    region_sky = regions.CircleSkyRegion(center=center_sky, radius=radius * u.arcsec)
    return region_sky.serialize(format='ds9')


    
if __name__ == "__main__":
  #x = source_region("Arcturus",ofn="test.reg")    
  #print(x)
  pass
