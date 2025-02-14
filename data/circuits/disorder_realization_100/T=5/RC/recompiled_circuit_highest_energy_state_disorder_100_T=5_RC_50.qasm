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
rz(1.963653) q[0];
sx q[0];
rz(-0.093129245) q[0];
sx q[0];
rz(-2.4094474) q[0];
rz(-1.5692476) q[1];
sx q[1];
rz(-0.88580004) q[1];
sx q[1];
rz(-2.7028309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2193091) q[0];
sx q[0];
rz(-1.6262591) q[0];
sx q[0];
rz(1.903141) q[0];
rz(-pi) q[1];
rz(2.8813527) q[2];
sx q[2];
rz(-1.2965045) q[2];
sx q[2];
rz(0.14014527) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6631515) q[1];
sx q[1];
rz(-1.3125129) q[1];
sx q[1];
rz(1.6186991) q[1];
x q[2];
rz(0.4583764) q[3];
sx q[3];
rz(-1.5273558) q[3];
sx q[3];
rz(0.14940748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0977122) q[2];
sx q[2];
rz(-3.1011797) q[2];
sx q[2];
rz(2.2229693) q[2];
rz(1.0527323) q[3];
sx q[3];
rz(-0.91679263) q[3];
sx q[3];
rz(0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87616462) q[0];
sx q[0];
rz(-0.81887236) q[0];
sx q[0];
rz(-2.2692666) q[0];
rz(-2.1288952) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(2.8025119) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10386183) q[0];
sx q[0];
rz(-1.4418238) q[0];
sx q[0];
rz(-0.43855389) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0929957) q[2];
sx q[2];
rz(-1.9777918) q[2];
sx q[2];
rz(-2.5841449) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87552947) q[1];
sx q[1];
rz(-1.5632946) q[1];
sx q[1];
rz(0.041635978) q[1];
rz(0.12437625) q[3];
sx q[3];
rz(-1.2487097) q[3];
sx q[3];
rz(0.45201221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33112153) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(2.5362711) q[2];
rz(0.33453861) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(-0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89638585) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(-1.6687923) q[0];
rz(2.8145166) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(-0.27007857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.208858) q[0];
sx q[0];
rz(-1.7678102) q[0];
sx q[0];
rz(-2.6126562) q[0];
rz(-pi) q[1];
rz(-1.7350586) q[2];
sx q[2];
rz(-0.43385071) q[2];
sx q[2];
rz(2.1151154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1835598) q[1];
sx q[1];
rz(-2.2366675) q[1];
sx q[1];
rz(-2.141365) q[1];
x q[2];
rz(1.5746501) q[3];
sx q[3];
rz(-1.4906297) q[3];
sx q[3];
rz(-1.7146669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.104287) q[2];
sx q[2];
rz(-1.9992) q[2];
sx q[2];
rz(1.8541065) q[2];
rz(-1.9262975) q[3];
sx q[3];
rz(-1.3853962) q[3];
sx q[3];
rz(-0.17770411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.52008587) q[0];
sx q[0];
rz(-2.6984213) q[0];
sx q[0];
rz(-2.8218414) q[0];
rz(2.7301835) q[1];
sx q[1];
rz(-1.5965867) q[1];
sx q[1];
rz(-3.06126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92341766) q[0];
sx q[0];
rz(-1.8785155) q[0];
sx q[0];
rz(-0.20826343) q[0];
rz(0.14239399) q[2];
sx q[2];
rz(-2.3300397) q[2];
sx q[2];
rz(2.591557) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7388725) q[1];
sx q[1];
rz(-0.79258535) q[1];
sx q[1];
rz(1.9566105) q[1];
x q[2];
rz(-2.5307054) q[3];
sx q[3];
rz(-1.3502933) q[3];
sx q[3];
rz(1.0569416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8722998) q[2];
sx q[2];
rz(-2.9041957) q[2];
sx q[2];
rz(-0.044895127) q[2];
rz(-0.53523713) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1503247) q[0];
sx q[0];
rz(-0.70882216) q[0];
sx q[0];
rz(2.2659361) q[0];
rz(-0.21370299) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(-1.5692086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0914577) q[0];
sx q[0];
rz(-2.7612855) q[0];
sx q[0];
rz(-2.1652392) q[0];
rz(-pi) q[1];
rz(-2.6781668) q[2];
sx q[2];
rz(-2.6604386) q[2];
sx q[2];
rz(-1.1347186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6900269) q[1];
sx q[1];
rz(-0.70575778) q[1];
sx q[1];
rz(-0.16132055) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8205958) q[3];
sx q[3];
rz(-1.4813652) q[3];
sx q[3];
rz(0.018370779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18465061) q[2];
sx q[2];
rz(-0.1447548) q[2];
sx q[2];
rz(-2.0733898) q[2];
rz(-0.60643658) q[3];
sx q[3];
rz(-1.2111827) q[3];
sx q[3];
rz(0.010312168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73685) q[0];
sx q[0];
rz(-2.8917942) q[0];
sx q[0];
rz(1.7506208) q[0];
rz(2.6178316) q[1];
sx q[1];
rz(-2.5660089) q[1];
sx q[1];
rz(-1.9578804) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7203248) q[0];
sx q[0];
rz(-2.8479332) q[0];
sx q[0];
rz(1.378242) q[0];
x q[1];
rz(2.1864088) q[2];
sx q[2];
rz(-3.0490626) q[2];
sx q[2];
rz(-1.9513771) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5038915) q[1];
sx q[1];
rz(-0.63234416) q[1];
sx q[1];
rz(-1.31193) q[1];
rz(-2.2565485) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(0.96641738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70048731) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(-1.0318476) q[2];
rz(-1.6796238) q[3];
sx q[3];
rz(-0.69457355) q[3];
sx q[3];
rz(1.3612548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4782891) q[0];
sx q[0];
rz(-0.42600584) q[0];
sx q[0];
rz(-1.7167094) q[0];
rz(2.5583963) q[1];
sx q[1];
rz(-1.9164663) q[1];
sx q[1];
rz(0.057417631) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17744495) q[0];
sx q[0];
rz(-0.76745196) q[0];
sx q[0];
rz(1.5657809) q[0];
x q[1];
rz(-1.0591828) q[2];
sx q[2];
rz(-1.0213189) q[2];
sx q[2];
rz(-1.8715931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7023041) q[1];
sx q[1];
rz(-1.3953475) q[1];
sx q[1];
rz(-2.5006341) q[1];
rz(-pi) q[2];
rz(0.52773169) q[3];
sx q[3];
rz(-2.4411628) q[3];
sx q[3];
rz(-2.79984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8251557) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(-1.3081029) q[2];
rz(1.8209275) q[3];
sx q[3];
rz(-0.84669176) q[3];
sx q[3];
rz(0.86110419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2945781) q[0];
sx q[0];
rz(-0.68462831) q[0];
sx q[0];
rz(1.8735877) q[0];
rz(1.1216724) q[1];
sx q[1];
rz(-1.2753762) q[1];
sx q[1];
rz(-0.65006382) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0065816) q[0];
sx q[0];
rz(-2.0913634) q[0];
sx q[0];
rz(-0.4659981) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5154326) q[2];
sx q[2];
rz(-0.78842064) q[2];
sx q[2];
rz(-1.5453218) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44733455) q[1];
sx q[1];
rz(-1.9617394) q[1];
sx q[1];
rz(-2.5503134) q[1];
rz(2.1459747) q[3];
sx q[3];
rz(-1.1523968) q[3];
sx q[3];
rz(-2.0783238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0012297) q[2];
sx q[2];
rz(-1.4823806) q[2];
sx q[2];
rz(-2.7492827) q[2];
rz(-2.9366734) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(0.71995455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1392764) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(-1.683715) q[0];
rz(-0.32683364) q[1];
sx q[1];
rz(-1.146233) q[1];
sx q[1];
rz(1.0439509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8736105) q[0];
sx q[0];
rz(-1.3820433) q[0];
sx q[0];
rz(2.0006936) q[0];
rz(-2.6828917) q[2];
sx q[2];
rz(-1.6965116) q[2];
sx q[2];
rz(-1.9551376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3957246) q[1];
sx q[1];
rz(-1.9650794) q[1];
sx q[1];
rz(3.0529725) q[1];
x q[2];
rz(-0.92071988) q[3];
sx q[3];
rz(-1.3082005) q[3];
sx q[3];
rz(2.289734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8331208) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(0.92588818) q[2];
rz(-2.0738257) q[3];
sx q[3];
rz(-1.1238778) q[3];
sx q[3];
rz(-1.7759751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85497102) q[0];
sx q[0];
rz(-1.4737031) q[0];
sx q[0];
rz(0.58706748) q[0];
rz(0.58285561) q[1];
sx q[1];
rz(-2.8522377) q[1];
sx q[1];
rz(-2.6606182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078048006) q[0];
sx q[0];
rz(-2.6780811) q[0];
sx q[0];
rz(2.6592451) q[0];
x q[1];
rz(1.9312385) q[2];
sx q[2];
rz(-1.3925465) q[2];
sx q[2];
rz(-2.0958063) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96718699) q[1];
sx q[1];
rz(-0.71651006) q[1];
sx q[1];
rz(-0.67542507) q[1];
x q[2];
rz(1.770788) q[3];
sx q[3];
rz(-1.1781577) q[3];
sx q[3];
rz(0.61905608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19758548) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(1.8120922) q[2];
rz(-0.78209376) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(-1.9935002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44163497) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(-1.6089454) q[1];
sx q[1];
rz(-1.9736704) q[1];
sx q[1];
rz(-0.71221487) q[1];
rz(-0.69923007) q[2];
sx q[2];
rz(-0.86730994) q[2];
sx q[2];
rz(2.2326474) q[2];
rz(-2.6431177) q[3];
sx q[3];
rz(-1.7117906) q[3];
sx q[3];
rz(-1.6556647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
