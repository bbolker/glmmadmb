install: glmmadmb
	mv glmmadmb ../bin/linux32

exe: glmmadmb.tpl
	admb -r glmmadmb

clean:
	rm -f glmmadmb.htp glmmadmb.o glmmadmb.cpp b1 b2 s1 s2 eigv.rpt variance fmin.log
