OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(1.7280581) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.04714) q[0];
sx q[0];
rz(-1.5507076) q[0];
sx q[0];
rz(1.8909341) q[0];
rz(-pi) q[1];
rz(1.9102664) q[2];
sx q[2];
rz(-1.4541153) q[2];
sx q[2];
rz(1.1543903) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0800228) q[1];
sx q[1];
rz(-2.1407619) q[1];
sx q[1];
rz(-1.391173) q[1];
rz(-pi) q[2];
rz(-2.324027) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(-0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(-3.1175933) q[2];
rz(0.35605797) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(0.63496494) q[0];
rz(0.7629281) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-2.2244804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13697727) q[0];
sx q[0];
rz(-3.1023623) q[0];
sx q[0];
rz(2.0580388) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88998743) q[2];
sx q[2];
rz(-1.1900848) q[2];
sx q[2];
rz(-3.0842881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8134384) q[1];
sx q[1];
rz(-2.0106689) q[1];
sx q[1];
rz(0.085170345) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8648326) q[3];
sx q[3];
rz(-0.64084478) q[3];
sx q[3];
rz(-1.3230192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1941173) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(2.4083162) q[2];
rz(1.0670079) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5697524) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(-2.6133614) q[0];
rz(-0.86604467) q[1];
sx q[1];
rz(-1.8128315) q[1];
sx q[1];
rz(2.6285062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24241867) q[0];
sx q[0];
rz(-2.0089825) q[0];
sx q[0];
rz(-3.0946381) q[0];
rz(-pi) q[1];
rz(2.8305574) q[2];
sx q[2];
rz(-2.1967109) q[2];
sx q[2];
rz(-0.011510465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.60050234) q[1];
sx q[1];
rz(-1.2230495) q[1];
sx q[1];
rz(-3.0060124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5993027) q[3];
sx q[3];
rz(-1.0712726) q[3];
sx q[3];
rz(-1.0699501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(0.038724381) q[2];
rz(-0.45768467) q[3];
sx q[3];
rz(-2.3364412) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(-1.1210972) q[0];
rz(-0.7577678) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(-2.4574492) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093903784) q[0];
sx q[0];
rz(-0.8465811) q[0];
sx q[0];
rz(0.27064497) q[0];
rz(-pi) q[1];
rz(-3.1314881) q[2];
sx q[2];
rz(-1.5636383) q[2];
sx q[2];
rz(2.7814076) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0348548) q[1];
sx q[1];
rz(-1.8632006) q[1];
sx q[1];
rz(0.3180983) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97400333) q[3];
sx q[3];
rz(-1.1349548) q[3];
sx q[3];
rz(-1.7929994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.1476436) q[3];
sx q[3];
rz(0.014613541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358661) q[0];
sx q[0];
rz(-1.7638693) q[0];
sx q[0];
rz(-2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(1.6872663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73390612) q[0];
sx q[0];
rz(-0.21557325) q[0];
sx q[0];
rz(0.13224299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5227484) q[2];
sx q[2];
rz(-1.4312051) q[2];
sx q[2];
rz(-1.0447811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1087423) q[1];
sx q[1];
rz(-1.8287484) q[1];
sx q[1];
rz(2.5775419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84564836) q[3];
sx q[3];
rz(-2.925784) q[3];
sx q[3];
rz(-2.8874021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5877567) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(-0.25911123) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(-2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3542341) q[0];
sx q[0];
rz(-0.43208313) q[0];
sx q[0];
rz(2.4969192) q[0];
rz(-1.8395754) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(2.7929746) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520599) q[0];
sx q[0];
rz(-2.4675998) q[0];
sx q[0];
rz(-2.3537082) q[0];
x q[1];
rz(0.63747823) q[2];
sx q[2];
rz(-2.5036252) q[2];
sx q[2];
rz(-1.5359985) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8585414) q[1];
sx q[1];
rz(-1.9894454) q[1];
sx q[1];
rz(-2.8178697) q[1];
x q[2];
rz(1.5362605) q[3];
sx q[3];
rz(-1.4538962) q[3];
sx q[3];
rz(-1.7021321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-2.1861173) q[2];
sx q[2];
rz(-3.1032584) q[2];
rz(0.96758715) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2751145) q[0];
sx q[0];
rz(-0.51946467) q[0];
sx q[0];
rz(-0.073609784) q[0];
rz(1.4069125) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(-2.4085192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778962) q[0];
sx q[0];
rz(-2.1025582) q[0];
sx q[0];
rz(0.97490208) q[0];
x q[1];
rz(2.5400852) q[2];
sx q[2];
rz(-2.2081293) q[2];
sx q[2];
rz(-1.21797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0907369) q[1];
sx q[1];
rz(-2.6466718) q[1];
sx q[1];
rz(-0.75116091) q[1];
x q[2];
rz(0.2352318) q[3];
sx q[3];
rz(-1.3266272) q[3];
sx q[3];
rz(-1.5489123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0343895) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(-3.1305195) q[2];
rz(-0.30558807) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(-1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082064) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(2.407684) q[0];
rz(2.5607196) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(3.0427921) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042511333) q[0];
sx q[0];
rz(-2.3020199) q[0];
sx q[0];
rz(0.73584105) q[0];
rz(-pi) q[1];
rz(0.96510624) q[2];
sx q[2];
rz(-0.98585502) q[2];
sx q[2];
rz(-0.044986812) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4186814) q[1];
sx q[1];
rz(-2.0127814) q[1];
sx q[1];
rz(1.8479455) q[1];
rz(-pi) q[2];
rz(-3.0116073) q[3];
sx q[3];
rz(-1.49125) q[3];
sx q[3];
rz(1.6787488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29844478) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(-1.0996381) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(-2.5449424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(1.9644568) q[0];
rz(-1.6674532) q[1];
sx q[1];
rz(-2.2087704) q[1];
sx q[1];
rz(-1.4449545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6048353) q[0];
sx q[0];
rz(-2.6492181) q[0];
sx q[0];
rz(-2.5144666) q[0];
rz(0.58622245) q[2];
sx q[2];
rz(-0.86044035) q[2];
sx q[2];
rz(-1.9682056) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8813546) q[1];
sx q[1];
rz(-1.5020276) q[1];
sx q[1];
rz(-2.0860904) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9111432) q[3];
sx q[3];
rz(-0.7045463) q[3];
sx q[3];
rz(-1.4662389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9745447) q[2];
sx q[2];
rz(-1.5571152) q[2];
sx q[2];
rz(2.9696999) q[2];
rz(2.3520172) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(1.8062228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563141) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(-0.54157448) q[0];
rz(1.9821292) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(2.3084739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1492958) q[0];
sx q[0];
rz(-2.160666) q[0];
sx q[0];
rz(1.9043395) q[0];
rz(-pi) q[1];
rz(2.59899) q[2];
sx q[2];
rz(-1.2719947) q[2];
sx q[2];
rz(1.8712559) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6852452) q[1];
sx q[1];
rz(-2.4083071) q[1];
sx q[1];
rz(-1.9128146) q[1];
rz(-pi) q[2];
rz(3.1396486) q[3];
sx q[3];
rz(-0.84315791) q[3];
sx q[3];
rz(1.9950642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1200166) q[2];
sx q[2];
rz(-2.3638201) q[2];
sx q[2];
rz(0.94338256) q[2];
rz(1.0776445) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(2.8130892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(-2.9019451) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(0.69375615) q[2];
sx q[2];
rz(-2.0907331) q[2];
sx q[2];
rz(1.5246684) q[2];
rz(-2.4619674) q[3];
sx q[3];
rz(-2.7503895) q[3];
sx q[3];
rz(2.214307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
