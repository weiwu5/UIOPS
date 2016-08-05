function insec=hhmmss2insec(timehhmmss)

insec = floor(timehhmmss/10000)*3600+floor(mod(timehhmmss,10000)/100)*60+floor(mod(timehhmmss,100));

end