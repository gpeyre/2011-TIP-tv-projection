FILE = FadiliPeyreTVProj

postscript:
	make -f Makefile tex
	make -f Makefile dvips

pdf:
	make -f Makefile tex
	make -f Makefile dvips
	make -f Makefile ps2pdf
	rm ${FILE}.ps
tex:
	latex ${FILE}.tex
dvips:
	dvips -o ${FILE}.ps -Ppdf -G0 ${FILE}.dvi
ps2pdf:
	ps2pdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress ${FILE}.ps
clean:
	rm -f ${FILE}.{aux,dvi,log,los,out,tpm}
