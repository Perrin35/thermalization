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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0005635) q[0];
sx q[0];
rz(-1.3357497) q[0];
sx q[0];
rz(1.1147333) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1821177) q[2];
sx q[2];
rz(-0.76180327) q[2];
sx q[2];
rz(-2.3488059) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71434778) q[1];
sx q[1];
rz(-2.3934919) q[1];
sx q[1];
rz(-2.1022878) q[1];
rz(-pi) q[2];
rz(0.58310572) q[3];
sx q[3];
rz(-2.8994377) q[3];
sx q[3];
rz(-0.21447578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(2.3857462) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(-2.1006179) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(1.2874228) q[0];
rz(2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(-1.2189254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201976) q[0];
sx q[0];
rz(-1.7292062) q[0];
sx q[0];
rz(1.3729457) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67181113) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(-1.1900657) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.375735) q[1];
sx q[1];
rz(-1.2925783) q[1];
sx q[1];
rz(2.9912659) q[1];
rz(-pi) q[2];
rz(0.35393012) q[3];
sx q[3];
rz(-0.90995379) q[3];
sx q[3];
rz(-0.61001813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-2.5197869) q[2];
rz(2.5032737) q[3];
sx q[3];
rz(-0.74898762) q[3];
sx q[3];
rz(2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(-2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-0.52282202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030634681) q[0];
sx q[0];
rz(-1.2485663) q[0];
sx q[0];
rz(-2.8186174) q[0];
x q[1];
rz(-0.39117809) q[2];
sx q[2];
rz(-1.4257981) q[2];
sx q[2];
rz(0.36996335) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.712341) q[1];
sx q[1];
rz(-2.0852226) q[1];
sx q[1];
rz(1.8271258) q[1];
rz(-2.3159695) q[3];
sx q[3];
rz(-0.43681991) q[3];
sx q[3];
rz(-0.88479155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(-1.5175021) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-2.2113694) q[3];
sx q[3];
rz(1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9972123) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(-1.602518) q[0];
rz(2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(1.296952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.813445) q[0];
sx q[0];
rz(-1.362365) q[0];
sx q[0];
rz(0.84020241) q[0];
rz(-pi) q[1];
rz(1.6455075) q[2];
sx q[2];
rz(-0.69715188) q[2];
sx q[2];
rz(-0.32767228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32758157) q[1];
sx q[1];
rz(-1.1594698) q[1];
sx q[1];
rz(0.41366215) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8861214) q[3];
sx q[3];
rz(-1.0389757) q[3];
sx q[3];
rz(-1.1602311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(-0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8056718) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(0.94201159) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(-1.6090144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.596622) q[0];
sx q[0];
rz(-0.76966296) q[0];
sx q[0];
rz(-0.53588698) q[0];
x q[1];
rz(-2.9608594) q[2];
sx q[2];
rz(-1.4305947) q[2];
sx q[2];
rz(3.0716346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11621257) q[1];
sx q[1];
rz(-1.0702225) q[1];
sx q[1];
rz(-0.0059555014) q[1];
x q[2];
rz(2.4306253) q[3];
sx q[3];
rz(-0.14440726) q[3];
sx q[3];
rz(2.292423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4970826) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(0.55822462) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(-0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(2.5010342) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(2.9277149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86588496) q[0];
sx q[0];
rz(-1.8731706) q[0];
sx q[0];
rz(0.23909607) q[0];
x q[1];
rz(-0.74486129) q[2];
sx q[2];
rz(-1.2861041) q[2];
sx q[2];
rz(-0.25279754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7678309) q[1];
sx q[1];
rz(-1.5015748) q[1];
sx q[1];
rz(0.94061416) q[1];
rz(1.6951896) q[3];
sx q[3];
rz(-0.52342192) q[3];
sx q[3];
rz(-1.0833797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(-1.7443582) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(-0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419256) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(1.4615321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7992226) q[0];
sx q[0];
rz(-0.91081753) q[0];
sx q[0];
rz(-2.2837385) q[0];
rz(-pi) q[1];
rz(-2.2093532) q[2];
sx q[2];
rz(-1.9912212) q[2];
sx q[2];
rz(-2.1115542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4883621) q[1];
sx q[1];
rz(-1.4547336) q[1];
sx q[1];
rz(0.56452063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6631132) q[3];
sx q[3];
rz(-2.5095486) q[3];
sx q[3];
rz(-1.5173591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(0.20480569) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(-1.9974476) q[1];
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
rz(-1.8587684) q[2];
sx q[2];
rz(-0.33703732) q[2];
sx q[2];
rz(0.84442838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42061372) q[1];
sx q[1];
rz(-0.79809626) q[1];
sx q[1];
rz(0.53687232) q[1];
rz(-0.89544483) q[3];
sx q[3];
rz(-0.43206462) q[3];
sx q[3];
rz(-0.91401446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-2.2224865) q[3];
sx q[3];
rz(2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018472483) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(1.8852604) q[0];
rz(-3.046335) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(0.5307861) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7408001) q[0];
sx q[0];
rz(-0.2966899) q[0];
sx q[0];
rz(-2.147697) q[0];
rz(-pi) q[1];
rz(1.5321391) q[2];
sx q[2];
rz(-1.2529904) q[2];
sx q[2];
rz(0.048887756) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3717613) q[1];
sx q[1];
rz(-0.27250817) q[1];
sx q[1];
rz(0.49401562) q[1];
x q[2];
rz(1.68566) q[3];
sx q[3];
rz(-1.126976) q[3];
sx q[3];
rz(2.165909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8119767) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(-0.81816188) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5166017) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(-0.3987819) q[0];
rz(1.7575691) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(-0.049364518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731664) q[0];
sx q[0];
rz(-3.0169562) q[0];
sx q[0];
rz(1.0602555) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3081506) q[2];
sx q[2];
rz(-1.6742953) q[2];
sx q[2];
rz(-1.7067133) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2422707) q[1];
sx q[1];
rz(-0.99192109) q[1];
sx q[1];
rz(-0.053236628) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27487288) q[3];
sx q[3];
rz(-0.65023732) q[3];
sx q[3];
rz(2.6848328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265726) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(1.0749558) q[2];
sx q[2];
rz(-1.1100162) q[2];
sx q[2];
rz(0.77592862) q[2];
rz(-0.31872411) q[3];
sx q[3];
rz(-2.095185) q[3];
sx q[3];
rz(-0.17879055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
