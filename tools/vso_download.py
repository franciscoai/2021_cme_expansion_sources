import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a

attrs_time = a.Time('2005/01/01 00:10', '2005/01/01 00:15')
result = Fido.search(attrs_time, a.Instrument.eit)

print(result)

downloaded_files = Fido.fetch(result, path='/gehme/scratch/cremades')
print(downloaded_files)