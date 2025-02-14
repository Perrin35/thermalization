OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.339191) q[0];
sx q[0];
rz(-1.1561166) q[0];
sx q[0];
rz(-2.6240786) q[0];
rz(1.327688) q[1];
sx q[1];
rz(-2.2496536) q[1];
sx q[1];
rz(-2.5995624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5919246) q[0];
sx q[0];
rz(-2.3506563) q[0];
sx q[0];
rz(2.8998081) q[0];
x q[1];
rz(1.8250685) q[2];
sx q[2];
rz(-0.82077998) q[2];
sx q[2];
rz(-0.17707846) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71859502) q[1];
sx q[1];
rz(-1.7684775) q[1];
sx q[1];
rz(-0.79841465) q[1];
x q[2];
rz(0.63707085) q[3];
sx q[3];
rz(-0.28709295) q[3];
sx q[3];
rz(-1.4305103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0077670495) q[2];
sx q[2];
rz(-1.0453036) q[2];
sx q[2];
rz(-0.93057752) q[2];
rz(0.13115701) q[3];
sx q[3];
rz(-0.37808642) q[3];
sx q[3];
rz(-1.4906496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59830484) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(0.92567956) q[1];
sx q[1];
rz(-1.8857748) q[1];
sx q[1];
rz(-1.6474887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0605436) q[0];
sx q[0];
rz(-0.88969066) q[0];
sx q[0];
rz(-0.42411719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41183128) q[2];
sx q[2];
rz(-2.4253325) q[2];
sx q[2];
rz(-1.0541944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32116649) q[1];
sx q[1];
rz(-1.120472) q[1];
sx q[1];
rz(1.1857595) q[1];
rz(-pi) q[2];
rz(2.8374313) q[3];
sx q[3];
rz(-1.5790325) q[3];
sx q[3];
rz(-1.7913763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1458448) q[2];
sx q[2];
rz(-0.44439849) q[2];
sx q[2];
rz(0.90163976) q[2];
rz(1.3580648) q[3];
sx q[3];
rz(-1.8243022) q[3];
sx q[3];
rz(-0.54734126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3667592) q[0];
sx q[0];
rz(-1.9623373) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(-2.215812) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-2.944223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.853123) q[0];
sx q[0];
rz(-0.98308167) q[0];
sx q[0];
rz(-2.9766686) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80855753) q[2];
sx q[2];
rz(-2.7312091) q[2];
sx q[2];
rz(1.9661599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0377522) q[1];
sx q[1];
rz(-1.7985797) q[1];
sx q[1];
rz(0.040018602) q[1];
x q[2];
rz(1.7719384) q[3];
sx q[3];
rz(-1.0363058) q[3];
sx q[3];
rz(2.8902658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7388514) q[2];
sx q[2];
rz(-2.1612942) q[2];
sx q[2];
rz(0.19035467) q[2];
rz(0.061577408) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(-2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9848118) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(-2.5132827) q[0];
rz(1.5208987) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(-0.77888387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951866) q[0];
sx q[0];
rz(-2.4281326) q[0];
sx q[0];
rz(-0.97795217) q[0];
x q[1];
rz(-1.2912108) q[2];
sx q[2];
rz(-1.2233073) q[2];
sx q[2];
rz(2.8742032) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4877704) q[1];
sx q[1];
rz(-1.5842591) q[1];
sx q[1];
rz(3.0938992) q[1];
x q[2];
rz(0.811399) q[3];
sx q[3];
rz(-0.86199844) q[3];
sx q[3];
rz(-0.78042316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4817619) q[2];
sx q[2];
rz(-1.0783106) q[2];
sx q[2];
rz(-4.3241186e-05) q[2];
rz(-2.0569885) q[3];
sx q[3];
rz(-1.9856039) q[3];
sx q[3];
rz(2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4701009) q[0];
sx q[0];
rz(-1.6974314) q[0];
sx q[0];
rz(1.2985562) q[0];
rz(2.5647054) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(2.6947122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53696886) q[0];
sx q[0];
rz(-2.4659221) q[0];
sx q[0];
rz(-2.8506822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2715029) q[2];
sx q[2];
rz(-1.6422049) q[2];
sx q[2];
rz(1.4324607) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.943534) q[1];
sx q[1];
rz(-1.7732278) q[1];
sx q[1];
rz(1.1620896) q[1];
rz(-1.5112259) q[3];
sx q[3];
rz(-1.7559663) q[3];
sx q[3];
rz(-2.6310754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7070861) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(2.9217829) q[2];
rz(-1.6541727) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(0.69948227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952482) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(-2.2440946) q[0];
rz(-1.7706361) q[1];
sx q[1];
rz(-2.1015344) q[1];
sx q[1];
rz(0.8955566) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33226704) q[0];
sx q[0];
rz(-0.73609645) q[0];
sx q[0];
rz(-2.9311137) q[0];
x q[1];
rz(-2.3562548) q[2];
sx q[2];
rz(-0.78162748) q[2];
sx q[2];
rz(-1.8434032) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0296391) q[1];
sx q[1];
rz(-1.9432507) q[1];
sx q[1];
rz(1.2122494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79781161) q[3];
sx q[3];
rz(-2.5113002) q[3];
sx q[3];
rz(2.6137461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.36876496) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(0.15711288) q[2];
rz(2.469192) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(3.0290643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86033487) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(-1.6819287) q[0];
rz(2.7809987) q[1];
sx q[1];
rz(-1.6723885) q[1];
sx q[1];
rz(1.8416539) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111709) q[0];
sx q[0];
rz(-1.9864559) q[0];
sx q[0];
rz(1.4372197) q[0];
rz(-pi) q[1];
rz(1.1088013) q[2];
sx q[2];
rz(-1.0014152) q[2];
sx q[2];
rz(1.5608567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2092525) q[1];
sx q[1];
rz(-2.8289218) q[1];
sx q[1];
rz(0.31950848) q[1];
rz(-1.9770558) q[3];
sx q[3];
rz(-0.99525577) q[3];
sx q[3];
rz(1.6353324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2099057) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(-2.9295861) q[2];
rz(2.0906406) q[3];
sx q[3];
rz(-1.3764328) q[3];
sx q[3];
rz(3.0534548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.5365005) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(-0.17298803) q[0];
rz(1.1146924) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(-3.0317422) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6914586) q[0];
sx q[0];
rz(-2.5485253) q[0];
sx q[0];
rz(2.791138) q[0];
rz(-2.3187056) q[2];
sx q[2];
rz(-1.7639065) q[2];
sx q[2];
rz(-2.2597974) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6348656) q[1];
sx q[1];
rz(-2.8322729) q[1];
sx q[1];
rz(-0.97815467) q[1];
x q[2];
rz(-0.55446378) q[3];
sx q[3];
rz(-1.4687456) q[3];
sx q[3];
rz(-0.79342519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7383808) q[2];
sx q[2];
rz(-1.2028368) q[2];
sx q[2];
rz(0.66486764) q[2];
rz(0.46857771) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(0.81336462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-2.5316694) q[0];
sx q[0];
rz(0.74158057) q[0];
rz(1.954151) q[1];
sx q[1];
rz(-1.6148184) q[1];
sx q[1];
rz(2.5724519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305107) q[0];
sx q[0];
rz(-1.2449322) q[0];
sx q[0];
rz(1.8799856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9958565) q[2];
sx q[2];
rz(-2.6837641) q[2];
sx q[2];
rz(-0.65953461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4148767) q[1];
sx q[1];
rz(-2.0316213) q[1];
sx q[1];
rz(-2.0445092) q[1];
x q[2];
rz(-2.3145202) q[3];
sx q[3];
rz(-0.88445348) q[3];
sx q[3];
rz(-0.65580578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6471214) q[2];
sx q[2];
rz(-0.98860604) q[2];
sx q[2];
rz(-2.3753601) q[2];
rz(1.1726441) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(2.9197781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93429339) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(-2.9497414) q[0];
rz(2.3244997) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(0.70768913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1199214) q[0];
sx q[0];
rz(-1.3967529) q[0];
sx q[0];
rz(2.8928323) q[0];
x q[1];
rz(0.3336256) q[2];
sx q[2];
rz(-2.1364436) q[2];
sx q[2];
rz(-0.27496342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.094511112) q[1];
sx q[1];
rz(-2.1489177) q[1];
sx q[1];
rz(2.6972527) q[1];
x q[2];
rz(-0.97174208) q[3];
sx q[3];
rz(-1.069456) q[3];
sx q[3];
rz(2.3399692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2231458) q[2];
sx q[2];
rz(-2.196849) q[2];
sx q[2];
rz(2.4165912) q[2];
rz(-2.9265192) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(-0.71896368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216777) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(1.1014145) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(-3.0240622) q[2];
sx q[2];
rz(-2.1030269) q[2];
sx q[2];
rz(3.0602958) q[2];
rz(1.5982895) q[3];
sx q[3];
rz(-2.0381654) q[3];
sx q[3];
rz(0.5034133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
