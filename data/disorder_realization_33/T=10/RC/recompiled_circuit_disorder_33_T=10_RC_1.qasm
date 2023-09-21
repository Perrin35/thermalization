OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(1.8338058) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966272) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(-2.3636723) q[0];
x q[1];
rz(2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(3.090976) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1882602) q[1];
sx q[1];
rz(-2.7213875) q[1];
sx q[1];
rz(2.5145867) q[1];
x q[2];
rz(2.9458463) q[3];
sx q[3];
rz(-1.9276852) q[3];
sx q[3];
rz(-1.9407879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(2.1291389) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5179948) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(2.2395796) q[0];
rz(-pi) q[1];
rz(0.88044135) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(1.134269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20570457) q[1];
sx q[1];
rz(-1.7824031) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(-pi) q[2];
rz(2.2651477) q[3];
sx q[3];
rz(-1.0507686) q[3];
sx q[3];
rz(0.84393535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(-1.8151981) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1754477) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.8885814) q[0];
rz(0.18205299) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(-1.6857266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4438666) q[1];
sx q[1];
rz(-1.4523456) q[1];
sx q[1];
rz(0.5603793) q[1];
rz(2.695735) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(-0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-0.50022593) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.4979699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7786176) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(2.1182563) q[0];
x q[1];
rz(-1.8565606) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(2.6464268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3264309) q[1];
sx q[1];
rz(-1.3109428) q[1];
sx q[1];
rz(1.3304779) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6763776) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(-2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(1.6197846) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(1.048208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74131504) q[0];
sx q[0];
rz(-1.8198697) q[0];
sx q[0];
rz(-1.2988017) q[0];
rz(-0.016096073) q[2];
sx q[2];
rz(-2.490009) q[2];
sx q[2];
rz(-2.1799257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9427467) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(1.320977) q[1];
rz(-2.1861595) q[3];
sx q[3];
rz(-1.9121998) q[3];
sx q[3];
rz(-1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(-2.5559588) q[0];
x q[1];
rz(-2.1075222) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(1.5922286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9784769) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-2.5549868) q[1];
rz(0.31452175) q[3];
sx q[3];
rz(-2.5701227) q[3];
sx q[3];
rz(-2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(-2.9061785) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670358) q[0];
sx q[0];
rz(-2.0445163) q[0];
sx q[0];
rz(1.0780328) q[0];
x q[1];
rz(-2.3115736) q[2];
sx q[2];
rz(-1.1168715) q[2];
sx q[2];
rz(1.9759076) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9069179) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-3.0659552) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7191914) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(-0.15448031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(-1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.4461393) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(-1.0780225) q[0];
x q[1];
rz(1.1649706) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(0.42025987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20128076) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(-2.8708354) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73236671) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(-0.68294169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.2040899) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4708913) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(-0.41505138) q[0];
rz(-1.2665389) q[2];
sx q[2];
rz(-0.90196246) q[2];
sx q[2];
rz(-2.5077016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6955399) q[1];
sx q[1];
rz(-2.3791168) q[1];
sx q[1];
rz(-1.6474849) q[1];
rz(0.44378186) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(-1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.2822255) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(1.4984109) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949085) q[0];
sx q[0];
rz(-2.2241728) q[0];
sx q[0];
rz(-0.16434591) q[0];
rz(-2.1543703) q[2];
sx q[2];
rz(-2.3505031) q[2];
sx q[2];
rz(0.9466048) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0018113) q[1];
sx q[1];
rz(-1.0011295) q[1];
sx q[1];
rz(-3.0210178) q[1];
rz(-pi) q[2];
rz(-1.3528353) q[3];
sx q[3];
rz(-2.2721604) q[3];
sx q[3];
rz(-0.71818128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-0.031899115) q[2];
sx q[2];
rz(-2.1733641) q[2];
sx q[2];
rz(2.7410438) q[2];
rz(3.071143) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
