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
rz(2.763971) q[0];
sx q[0];
rz(-0.42855898) q[0];
sx q[0];
rz(0.23634401) q[0];
rz(1.4689245) q[1];
sx q[1];
rz(4.7279141) q[1];
sx q[1];
rz(9.5864819) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16888389) q[0];
sx q[0];
rz(-2.5492269) q[0];
sx q[0];
rz(-1.7667207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9457372) q[2];
sx q[2];
rz(-1.6278978) q[2];
sx q[2];
rz(1.3916935) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0031944) q[1];
sx q[1];
rz(-1.5439543) q[1];
sx q[1];
rz(3.1346727) q[1];
x q[2];
rz(1.3507648) q[3];
sx q[3];
rz(-1.2041948) q[3];
sx q[3];
rz(-2.0442968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24254313) q[2];
sx q[2];
rz(-2.264302) q[2];
sx q[2];
rz(-2.7522932) q[2];
rz(-0.23046514) q[3];
sx q[3];
rz(-0.018229818) q[3];
sx q[3];
rz(0.61483312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56735754) q[0];
sx q[0];
rz(-0.95139545) q[0];
sx q[0];
rz(-1.6472598) q[0];
rz(1.5556473) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(-1.1344604) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0527549) q[0];
sx q[0];
rz(-0.87273926) q[0];
sx q[0];
rz(-1.237117) q[0];
rz(-pi) q[1];
rz(-2.6642642) q[2];
sx q[2];
rz(-2.3466718) q[2];
sx q[2];
rz(-1.8177462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74090577) q[1];
sx q[1];
rz(-2.4016651) q[1];
sx q[1];
rz(0.94018971) q[1];
x q[2];
rz(2.1389009) q[3];
sx q[3];
rz(-2.3062996) q[3];
sx q[3];
rz(-2.7138674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.050934164) q[2];
sx q[2];
rz(-0.91973534) q[2];
sx q[2];
rz(-1.301379) q[2];
rz(2.0672412) q[3];
sx q[3];
rz(-0.31000546) q[3];
sx q[3];
rz(1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.018547) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(-2.5240335) q[0];
rz(-1.1005719) q[1];
sx q[1];
rz(-0.01958422) q[1];
sx q[1];
rz(0.40357959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8698533) q[0];
sx q[0];
rz(-1.4560486) q[0];
sx q[0];
rz(3.1310351) q[0];
rz(1.4694655) q[2];
sx q[2];
rz(-1.0472385) q[2];
sx q[2];
rz(0.061274139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8051648) q[1];
sx q[1];
rz(-1.4138828) q[1];
sx q[1];
rz(0.13586951) q[1];
x q[2];
rz(-0.37637122) q[3];
sx q[3];
rz(-2.5507002) q[3];
sx q[3];
rz(1.2749689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5201716) q[2];
sx q[2];
rz(-1.6941864) q[2];
sx q[2];
rz(-0.78110313) q[2];
rz(-2.805294) q[3];
sx q[3];
rz(-1.4399485) q[3];
sx q[3];
rz(-2.4380016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4562562) q[0];
sx q[0];
rz(-0.51613134) q[0];
sx q[0];
rz(-1.9235032) q[0];
rz(-1.3457899) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(-1.9075314) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505554) q[0];
sx q[0];
rz(-0.19007401) q[0];
sx q[0];
rz(0.051817373) q[0];
rz(-pi) q[1];
rz(-1.2786364) q[2];
sx q[2];
rz(-2.1641762) q[2];
sx q[2];
rz(2.0403543) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30013903) q[1];
sx q[1];
rz(-2.6965504) q[1];
sx q[1];
rz(2.1883394) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59438057) q[3];
sx q[3];
rz(-3.1295589) q[3];
sx q[3];
rz(1.8569225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1951083) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(1.9931603) q[2];
rz(-0.75442433) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(0.26328009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36963439) q[0];
sx q[0];
rz(-3.0535871) q[0];
sx q[0];
rz(-0.60246402) q[0];
rz(-0.98214904) q[1];
sx q[1];
rz(-3.1352477) q[1];
sx q[1];
rz(-0.51012653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5969491) q[0];
sx q[0];
rz(-1.7644797) q[0];
sx q[0];
rz(-0.18904833) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0544364) q[2];
sx q[2];
rz(-1.4498324) q[2];
sx q[2];
rz(2.6270514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1044961) q[1];
sx q[1];
rz(-1.949661) q[1];
sx q[1];
rz(-2.563743) q[1];
rz(1.5701775) q[3];
sx q[3];
rz(-2.3851756) q[3];
sx q[3];
rz(-1.7572614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.043857418) q[2];
sx q[2];
rz(-1.1172224) q[2];
sx q[2];
rz(-1.0350234) q[2];
rz(1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(2.5586832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3007091) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-2.050052) q[0];
rz(-0.40065271) q[1];
sx q[1];
rz(-3.1392097) q[1];
sx q[1];
rz(-0.56652743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3728722) q[0];
sx q[0];
rz(-2.3864288) q[0];
sx q[0];
rz(1.374097) q[0];
x q[1];
rz(-2.3344757) q[2];
sx q[2];
rz(-2.2645778) q[2];
sx q[2];
rz(-1.6859695) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2926082) q[1];
sx q[1];
rz(-1.0142393) q[1];
sx q[1];
rz(1.3054331) q[1];
rz(-pi) q[2];
rz(-2.5689396) q[3];
sx q[3];
rz(-1.6653456) q[3];
sx q[3];
rz(-1.270164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94622508) q[2];
sx q[2];
rz(-0.75104284) q[2];
sx q[2];
rz(-1.5283778) q[2];
rz(1.001312) q[3];
sx q[3];
rz(-2.2712205) q[3];
sx q[3];
rz(-2.9010469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85750759) q[0];
sx q[0];
rz(-2.3473098) q[0];
sx q[0];
rz(1.9965782) q[0];
rz(-1.5223632) q[1];
sx q[1];
rz(-0.018773627) q[1];
sx q[1];
rz(-1.1633263) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0534819) q[0];
sx q[0];
rz(-1.6170752) q[0];
sx q[0];
rz(1.2812216) q[0];
rz(0.14919632) q[2];
sx q[2];
rz(-0.90221221) q[2];
sx q[2];
rz(-0.72982349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12194828) q[1];
sx q[1];
rz(-2.5202529) q[1];
sx q[1];
rz(-1.7985733) q[1];
x q[2];
rz(0.16086732) q[3];
sx q[3];
rz(-0.64023521) q[3];
sx q[3];
rz(1.7750027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5160949) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(-2.7682448) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(-1.1628304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9173376) q[0];
sx q[0];
rz(-2.6118216) q[0];
sx q[0];
rz(0.24714558) q[0];
rz(0.22357926) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(-1.5162969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151536) q[0];
sx q[0];
rz(-1.5064539) q[0];
sx q[0];
rz(1.6872726) q[0];
rz(1.9886279) q[2];
sx q[2];
rz(-2.2925809) q[2];
sx q[2];
rz(-1.4773475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5943971) q[1];
sx q[1];
rz(-1.5548348) q[1];
sx q[1];
rz(-0.39401024) q[1];
rz(-3.0706625) q[3];
sx q[3];
rz(-0.59535691) q[3];
sx q[3];
rz(-1.8728674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.197) q[2];
sx q[2];
rz(-0.71114117) q[2];
sx q[2];
rz(1.5468583) q[2];
rz(-2.366015) q[3];
sx q[3];
rz(-1.2838793) q[3];
sx q[3];
rz(1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.813886) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(1.1073329) q[0];
rz(0.92512098) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(0.75606871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297345) q[0];
sx q[0];
rz(-1.2639771) q[0];
sx q[0];
rz(2.1351249) q[0];
rz(-pi) q[1];
rz(2.101641) q[2];
sx q[2];
rz(-1.565298) q[2];
sx q[2];
rz(-1.0571684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1518095) q[1];
sx q[1];
rz(-1.0149455) q[1];
sx q[1];
rz(1.226306) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1875137) q[3];
sx q[3];
rz(-2.5616259) q[3];
sx q[3];
rz(-0.43646508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5209311) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(-1.0010285) q[2];
rz(1.3641317) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(-2.6076243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022024632) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(-2.6764828) q[0];
rz(1.3348835) q[1];
sx q[1];
rz(-0.3723793) q[1];
sx q[1];
rz(-1.575527) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0571363) q[0];
sx q[0];
rz(-1.6142862) q[0];
sx q[0];
rz(1.3410596) q[0];
rz(-pi) q[1];
rz(-2.8303873) q[2];
sx q[2];
rz(-2.0796144) q[2];
sx q[2];
rz(2.6151163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3153957) q[1];
sx q[1];
rz(-0.0028571833) q[1];
sx q[1];
rz(-0.14551659) q[1];
rz(-pi) q[2];
rz(-0.09327831) q[3];
sx q[3];
rz(-0.3808379) q[3];
sx q[3];
rz(2.8990428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2472725) q[2];
sx q[2];
rz(-3.0970116) q[2];
sx q[2];
rz(-2.0559922) q[2];
rz(1.2224489) q[3];
sx q[3];
rz(-2.6294851) q[3];
sx q[3];
rz(-1.2781757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6923675) q[0];
sx q[0];
rz(-1.4341555) q[0];
sx q[0];
rz(1.7644802) q[0];
rz(1.5583246) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(-3.0408411) q[2];
sx q[2];
rz(-1.5660847) q[2];
sx q[2];
rz(-1.2790578) q[2];
rz(1.1261945) q[3];
sx q[3];
rz(-1.3570519) q[3];
sx q[3];
rz(-2.9828664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
