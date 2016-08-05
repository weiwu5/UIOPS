function timehhmmss=insec2hhmmss(insec)

timehhmmss = floor(insec/3600)*10000+floor(mod(insec,3600)/60)*100+floor(mod(insec,60));

end