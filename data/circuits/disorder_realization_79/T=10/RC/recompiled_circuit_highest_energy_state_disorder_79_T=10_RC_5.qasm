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
rz(1.6066983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0005635) q[0];
sx q[0];
rz(-1.3357497) q[0];
sx q[0];
rz(-2.0268593) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7947172) q[2];
sx q[2];
rz(-2.2636608) q[2];
sx q[2];
rz(-1.833806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1032627) q[1];
sx q[1];
rz(-0.94417773) q[1];
sx q[1];
rz(-0.4396529) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20333692) q[3];
sx q[3];
rz(-1.7032188) q[3];
sx q[3];
rz(1.2156957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-0.55880004) q[2];
sx q[2];
rz(0.99754769) q[2];
rz(-2.3857462) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-1.8541699) q[0];
rz(2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.9226673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12139509) q[0];
sx q[0];
rz(-1.7292062) q[0];
sx q[0];
rz(-1.7686469) q[0];
x q[1];
rz(-0.27147175) q[2];
sx q[2];
rz(-0.69005943) q[2];
sx q[2];
rz(-0.59218237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23663737) q[1];
sx q[1];
rz(-1.4262916) q[1];
sx q[1];
rz(-1.852024) q[1];
rz(1.1514436) q[3];
sx q[3];
rz(-0.73691982) q[3];
sx q[3];
rz(1.9896955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-0.62180579) q[2];
rz(2.5032737) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8353552) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(0.52282202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4958333) q[0];
sx q[0];
rz(-1.8765939) q[0];
sx q[0];
rz(-1.9093139) q[0];
rz(-0.36575138) q[2];
sx q[2];
rz(-0.41588441) q[2];
sx q[2];
rz(-1.5378086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.26979687) q[1];
sx q[1];
rz(-1.3482454) q[1];
sx q[1];
rz(0.52877626) q[1];
rz(-pi) q[2];
rz(-2.834972) q[3];
sx q[3];
rz(-1.2546179) q[3];
sx q[3];
rz(3.0512323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9598976) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(-1.5175021) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9972123) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(1.5390747) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(-1.296952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015681277) q[0];
sx q[0];
rz(-0.75443479) q[0];
sx q[0];
rz(-1.2638646) q[0];
rz(-pi) q[1];
rz(1.4960852) q[2];
sx q[2];
rz(-0.69715188) q[2];
sx q[2];
rz(-2.8139204) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9819543) q[1];
sx q[1];
rz(-0.57483638) q[1];
sx q[1];
rz(0.82623085) q[1];
x q[2];
rz(0.4850895) q[3];
sx q[3];
rz(-0.61044932) q[3];
sx q[3];
rz(1.7318673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(-2.0959334) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(-1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8056718) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(-0.94201159) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(-1.5325783) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9056775) q[0];
sx q[0];
rz(-0.92936838) q[0];
sx q[0];
rz(1.1113313) q[0];
x q[1];
rz(-2.9608594) q[2];
sx q[2];
rz(-1.7109979) q[2];
sx q[2];
rz(-3.0716346) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6841507) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(1.0702151) q[1];
rz(-pi) q[2];
rz(0.71096731) q[3];
sx q[3];
rz(-0.14440726) q[3];
sx q[3];
rz(0.84916964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64451009) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(-2.8572148) q[2];
rz(-0.55822462) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.072642) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(0.64055842) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(2.9277149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9616674) q[0];
sx q[0];
rz(-0.38319394) q[0];
sx q[0];
rz(0.92143329) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40753813) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(-2.119273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37376172) q[1];
sx q[1];
rz(-1.5015748) q[1];
sx q[1];
rz(0.94061416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6951896) q[3];
sx q[3];
rz(-0.52342192) q[3];
sx q[3];
rz(2.058213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(1.3972345) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419256) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(-1.1055111) q[0];
rz(-2.8583177) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.4615321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7532336) q[0];
sx q[0];
rz(-2.2112911) q[0];
sx q[0];
rz(2.4413013) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2093532) q[2];
sx q[2];
rz(-1.9912212) q[2];
sx q[2];
rz(1.0300385) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6532306) q[1];
sx q[1];
rz(-1.4547336) q[1];
sx q[1];
rz(0.56452063) q[1];
rz(-pi) q[2];
rz(-2.5652065) q[3];
sx q[3];
rz(-1.295305) q[3];
sx q[3];
rz(2.798693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0561515) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(-2.936787) q[2];
rz(1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(1.9974476) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56267101) q[0];
sx q[0];
rz(-1.5454588) q[0];
sx q[0];
rz(-2.6549005) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8587684) q[2];
sx q[2];
rz(-2.8045553) q[2];
sx q[2];
rz(-2.2971643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5976286) q[1];
sx q[1];
rz(-1.1958599) q[1];
sx q[1];
rz(0.72245325) q[1];
rz(-pi) q[2];
rz(1.9162634) q[3];
sx q[3];
rz(-1.3059214) q[3];
sx q[3];
rz(0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(1.8899274) q[2];
rz(-1.7898611) q[3];
sx q[3];
rz(-2.2224865) q[3];
sx q[3];
rz(2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(3.046335) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(0.5307861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149218) q[0];
sx q[0];
rz(-1.7309395) q[0];
sx q[0];
rz(1.3199575) q[0];
rz(0.31802788) q[2];
sx q[2];
rz(-1.5340759) q[2];
sx q[2];
rz(1.6317692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3717613) q[1];
sx q[1];
rz(-0.27250817) q[1];
sx q[1];
rz(0.49401562) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44638951) q[3];
sx q[3];
rz(-1.6744895) q[3];
sx q[3];
rz(-2.5959792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8119767) q[2];
sx q[2];
rz(-1.5584385) q[2];
sx q[2];
rz(0.54083332) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(-2.5087859) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5166017) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(-2.7428108) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-0.049364518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20959768) q[0];
sx q[0];
rz(-1.6315797) q[0];
sx q[0];
rz(-1.4619191) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0344459) q[2];
sx q[2];
rz(-1.309589) q[2];
sx q[2];
rz(-3.0334453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80222622) q[1];
sx q[1];
rz(-0.58103937) q[1];
sx q[1];
rz(1.4895579) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(-2.2137568) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(1.0646959) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(-2.3774316) q[2];
sx q[2];
rz(-0.66351009) q[2];
sx q[2];
rz(-1.4828975) q[2];
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
