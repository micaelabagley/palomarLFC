pro imalign, wispfield
; align images in IDL via hastrom.pro

; rename old i band image to keep for later
spawn, 'mv ' +wispfield+ '/' +wispfield+ '_i.fits ' +wispfield+ '/' $
        +wispfield+'_i_old.fits'
; read in images
g = mrdfits(wispfield+'/'+wispfield+'_g.fits', 0, hg)
oldi = mrdfits(wispfield+'/'+wispfield+'_i_old.fits', 0, hi)

; align i band to g band
hastrom, oldi, hi, newi, newhd, hg, interp=2, missing=0
; write image out
writefits, wispfield+'_i.fits', newi, newhd

return
end
