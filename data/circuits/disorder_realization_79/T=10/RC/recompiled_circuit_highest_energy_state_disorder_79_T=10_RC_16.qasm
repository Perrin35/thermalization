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
rz(0.41843721) q[0];
sx q[0];
rz(-0.96324459) q[0];
sx q[0];
rz(0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(-1.5348943) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14102916) q[0];
sx q[0];
rz(-1.805843) q[0];
sx q[0];
rz(-1.1147333) q[0];
rz(1.9594749) q[2];
sx q[2];
rz(-2.3797894) q[2];
sx q[2];
rz(-2.3488059) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71434778) q[1];
sx q[1];
rz(-2.3934919) q[1];
sx q[1];
rz(2.1022878) q[1];
x q[2];
rz(1.4356218) q[3];
sx q[3];
rz(-1.3692642) q[3];
sx q[3];
rz(-0.38231787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(-2.3857462) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.2874228) q[0];
rz(-2.6584794) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(-1.2189254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0201976) q[0];
sx q[0];
rz(-1.4123865) q[0];
sx q[0];
rz(1.3729457) q[0];
rz(-pi) q[1];
rz(-2.4697815) q[2];
sx q[2];
rz(-1.7423358) q[2];
sx q[2];
rz(-1.951527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9049553) q[1];
sx q[1];
rz(-1.4262916) q[1];
sx q[1];
rz(1.852024) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1514436) q[3];
sx q[3];
rz(-2.4046728) q[3];
sx q[3];
rz(-1.1518971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(0.62180579) q[2];
rz(-0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(-0.21632347) q[0];
rz(2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(0.52282202) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7824088) q[0];
sx q[0];
rz(-0.45219996) q[0];
sx q[0];
rz(-2.3308999) q[0];
rz(-pi) q[1];
rz(0.39117809) q[2];
sx q[2];
rz(-1.7157946) q[2];
sx q[2];
rz(-2.7716293) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4292517) q[1];
sx q[1];
rz(-2.0852226) q[1];
sx q[1];
rz(-1.3144668) q[1];
x q[2];
rz(0.82562311) q[3];
sx q[3];
rz(-2.7047727) q[3];
sx q[3];
rz(0.88479155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(-1.6240906) q[2];
rz(-0.46487871) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.6200861) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9972123) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(1.602518) q[0];
rz(2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(-1.8446406) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813445) q[0];
sx q[0];
rz(-1.362365) q[0];
sx q[0];
rz(-0.84020241) q[0];
rz(-pi) q[1];
rz(-0.06242604) q[2];
sx q[2];
rz(-2.2656144) q[2];
sx q[2];
rz(-2.7165987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8140111) q[1];
sx q[1];
rz(-1.1594698) q[1];
sx q[1];
rz(-0.41366215) q[1];
rz(-1.8861214) q[3];
sx q[3];
rz(-1.0389757) q[3];
sx q[3];
rz(-1.9813615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79232717) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(-1.0456592) q[2];
rz(-2.9271434) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3359208) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(2.1995811) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(1.5325783) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62270861) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(2.4469482) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4759859) q[2];
sx q[2];
rz(-0.22826787) q[2];
sx q[2];
rz(0.98787243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10380393) q[1];
sx q[1];
rz(-2.6409864) q[1];
sx q[1];
rz(1.5816825) q[1];
rz(-1.665409) q[3];
sx q[3];
rz(-1.4615371) q[3];
sx q[3];
rz(-3.008568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4970826) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-2.5010342) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(-0.21387771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17992526) q[0];
sx q[0];
rz(-2.7583987) q[0];
sx q[0];
rz(-0.92143329) q[0];
x q[1];
rz(-2.7340545) q[2];
sx q[2];
rz(-2.354051) q[2];
sx q[2];
rz(-2.119273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7678309) q[1];
sx q[1];
rz(-1.6400178) q[1];
sx q[1];
rz(2.2009785) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0908666) q[3];
sx q[3];
rz(-1.6328535) q[3];
sx q[3];
rz(-2.7620535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(0.17769979) q[2];
rz(1.3972345) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(-2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419256) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.4615321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34237) q[0];
sx q[0];
rz(-0.91081753) q[0];
sx q[0];
rz(-0.8578542) q[0];
rz(2.2143557) q[2];
sx q[2];
rz(-2.3935648) q[2];
sx q[2];
rz(-3.1035556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.098432861) q[1];
sx q[1];
rz(-2.5665356) q[1];
sx q[1];
rz(-0.21456031) q[1];
x q[2];
rz(0.57638611) q[3];
sx q[3];
rz(-1.295305) q[3];
sx q[3];
rz(-0.34289962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(2.936787) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(-1.144145) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5789216) q[0];
sx q[0];
rz(-1.5454588) q[0];
sx q[0];
rz(-2.6549005) q[0];
rz(-pi) q[1];
rz(-0.099191908) q[2];
sx q[2];
rz(-1.2481523) q[2];
sx q[2];
rz(-1.1485554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28539666) q[1];
sx q[1];
rz(-0.90803972) q[1];
sx q[1];
rz(1.0876571) q[1];
x q[2];
rz(1.9162634) q[3];
sx q[3];
rz(-1.3059214) q[3];
sx q[3];
rz(0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3860151) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(1.8899274) q[2];
rz(-1.7898611) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(1.2563323) q[0];
rz(-0.095257692) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(-0.5307861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40079257) q[0];
sx q[0];
rz(-0.2966899) q[0];
sx q[0];
rz(0.99389561) q[0];
rz(-pi) q[1];
rz(2.8235648) q[2];
sx q[2];
rz(-1.6075168) q[2];
sx q[2];
rz(1.6317692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3717613) q[1];
sx q[1];
rz(-2.8690845) q[1];
sx q[1];
rz(-2.647577) q[1];
x q[2];
rz(0.23654273) q[3];
sx q[3];
rz(-0.45748392) q[3];
sx q[3];
rz(1.2381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.329616) q[2];
sx q[2];
rz(-1.5584385) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(-2.3234308) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5166017) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(2.7428108) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(3.0922281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731664) q[0];
sx q[0];
rz(-3.0169562) q[0];
sx q[0];
rz(-2.0813372) q[0];
rz(-pi) q[1];
rz(-1.8334421) q[2];
sx q[2];
rz(-1.4672973) q[2];
sx q[2];
rz(1.4348794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4409742) q[1];
sx q[1];
rz(-1.6153533) q[1];
sx q[1];
rz(2.150321) q[1];
x q[2];
rz(-0.63189854) q[3];
sx q[3];
rz(-1.4057341) q[3];
sx q[3];
rz(-0.89323211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(-1.6309942) q[2];
rz(-0.92783582) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(2.0768968) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(-2.0666368) q[2];
sx q[2];
rz(-1.1100162) q[2];
sx q[2];
rz(0.77592862) q[2];
rz(2.06729) q[3];
sx q[3];
rz(-2.5357694) q[3];
sx q[3];
rz(2.3800935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
