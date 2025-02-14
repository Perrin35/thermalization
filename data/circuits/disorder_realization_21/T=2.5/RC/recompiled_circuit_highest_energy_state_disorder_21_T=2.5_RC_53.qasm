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
rz(1.709197) q[0];
sx q[0];
rz(3.4442918) q[0];
sx q[0];
rz(11.533307) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(-2.9472247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56563604) q[0];
sx q[0];
rz(-0.64467421) q[0];
sx q[0];
rz(2.8863532) q[0];
rz(-pi) q[1];
rz(1.871045) q[2];
sx q[2];
rz(-2.2451043) q[2];
sx q[2];
rz(2.6943494) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7558507) q[1];
sx q[1];
rz(-0.057690851) q[1];
sx q[1];
rz(0.92924802) q[1];
rz(-pi) q[2];
rz(-0.62567775) q[3];
sx q[3];
rz(-1.2114269) q[3];
sx q[3];
rz(-2.7070482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4091461) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(0.5644325) q[2];
rz(-1.8730646) q[3];
sx q[3];
rz(-3.1334435) q[3];
sx q[3];
rz(0.90659365) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0821575) q[0];
sx q[0];
rz(-0.0075639021) q[0];
sx q[0];
rz(2.2285158) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(0.93904644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093151) q[0];
sx q[0];
rz(-2.1255204) q[0];
sx q[0];
rz(-1.1251262) q[0];
rz(-3.1153283) q[2];
sx q[2];
rz(-2.677478) q[2];
sx q[2];
rz(2.8796706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8814358) q[1];
sx q[1];
rz(-1.9828567) q[1];
sx q[1];
rz(-2.7753633) q[1];
rz(-1.1640074) q[3];
sx q[3];
rz(-1.6322246) q[3];
sx q[3];
rz(0.2964501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82274503) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(1.5626296) q[3];
sx q[3];
rz(-3.1171418) q[3];
sx q[3];
rz(-1.8306556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398294) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(-2.7318562) q[0];
rz(0.01092625) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(-1.797537) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15191244) q[0];
sx q[0];
rz(-1.5912959) q[0];
sx q[0];
rz(2.0818039) q[0];
x q[1];
rz(3.0777099) q[2];
sx q[2];
rz(-1.6860355) q[2];
sx q[2];
rz(-1.2196397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2933767) q[1];
sx q[1];
rz(-2.8753706) q[1];
sx q[1];
rz(1.8195527) q[1];
rz(-3.1317409) q[3];
sx q[3];
rz(-1.6860699) q[3];
sx q[3];
rz(1.7719763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4695796) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(0.043188485) q[2];
rz(-2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(1.9388916) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7286872) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(-1.8034978) q[0];
rz(0.69522229) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-0.41740886) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090583853) q[0];
sx q[0];
rz(-0.90253297) q[0];
sx q[0];
rz(-2.9221852) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9193255) q[2];
sx q[2];
rz(-0.67643316) q[2];
sx q[2];
rz(-1.2471022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35170947) q[1];
sx q[1];
rz(-1.369006) q[1];
sx q[1];
rz(-0.92243299) q[1];
rz(1.1379477) q[3];
sx q[3];
rz(-1.4420849) q[3];
sx q[3];
rz(0.52279982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7030316) q[2];
sx q[2];
rz(-0.87623864) q[2];
sx q[2];
rz(1.38928) q[2];
rz(2.321068) q[3];
sx q[3];
rz(-1.5912143) q[3];
sx q[3];
rz(0.70782053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97838068) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(-2.1782844) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(0.16545573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2509569) q[0];
sx q[0];
rz(-1.5444618) q[0];
sx q[0];
rz(-0.043553003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9257912) q[2];
sx q[2];
rz(-0.7264464) q[2];
sx q[2];
rz(-2.1997228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7041118) q[1];
sx q[1];
rz(-0.81161122) q[1];
sx q[1];
rz(2.981217) q[1];
x q[2];
rz(0.31807301) q[3];
sx q[3];
rz(-1.9847365) q[3];
sx q[3];
rz(1.3596168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.307622) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(2.0818254) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(-1.8701657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6178013) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(1.6049438) q[0];
rz(2.6887584) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.7013928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0217845) q[0];
sx q[0];
rz(-2.9077466) q[0];
sx q[0];
rz(-0.77063693) q[0];
x q[1];
rz(-1.3511436) q[2];
sx q[2];
rz(-1.7872602) q[2];
sx q[2];
rz(-0.6714657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3523941) q[1];
sx q[1];
rz(-1.4161613) q[1];
sx q[1];
rz(-1.5633538) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4000638) q[3];
sx q[3];
rz(-1.472982) q[3];
sx q[3];
rz(2.3119237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77233934) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(0.079027979) q[2];
rz(-2.2760271) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(-2.1178093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(-0.22501568) q[0];
rz(-2.8354538) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(0.87047815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891312) q[0];
sx q[0];
rz(-1.5232183) q[0];
sx q[0];
rz(3.070288) q[0];
rz(-pi) q[1];
rz(2.5092431) q[2];
sx q[2];
rz(-0.82201695) q[2];
sx q[2];
rz(-2.1123304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5600109) q[1];
sx q[1];
rz(-2.3696179) q[1];
sx q[1];
rz(-0.08666579) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0848177) q[3];
sx q[3];
rz(-1.2013913) q[3];
sx q[3];
rz(1.891275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(-1.0352943) q[2];
rz(2.3546442) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(2.8662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(-3.1159478) q[0];
rz(-1.2643087) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(-0.69902507) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2866488) q[0];
sx q[0];
rz(-0.35084769) q[0];
sx q[0];
rz(0.69706608) q[0];
rz(1.4003631) q[2];
sx q[2];
rz(-2.0336069) q[2];
sx q[2];
rz(-1.544726) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5311218) q[1];
sx q[1];
rz(-1.4969331) q[1];
sx q[1];
rz(-0.38458891) q[1];
rz(-pi) q[2];
rz(2.2421942) q[3];
sx q[3];
rz(-2.3672383) q[3];
sx q[3];
rz(0.13860656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25253025) q[2];
sx q[2];
rz(-0.74873304) q[2];
sx q[2];
rz(0.32386455) q[2];
rz(2.9845386) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(0.15373716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57482982) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(0.59004849) q[0];
rz(-0.7198965) q[1];
sx q[1];
rz(-0.059559278) q[1];
sx q[1];
rz(-0.86729008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7035422) q[0];
sx q[0];
rz(-0.64213404) q[0];
sx q[0];
rz(-2.9657528) q[0];
rz(2.8810416) q[2];
sx q[2];
rz(-0.98439596) q[2];
sx q[2];
rz(1.9908277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1220951) q[1];
sx q[1];
rz(-1.6830793) q[1];
sx q[1];
rz(3.0730559) q[1];
rz(-pi) q[2];
rz(0.85764472) q[3];
sx q[3];
rz(-1.4271171) q[3];
sx q[3];
rz(1.9471256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94816339) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(-0.11290045) q[2];
rz(-2.0924977) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.0650487) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(1.5034058) q[0];
rz(-2.5050971) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(-1.5719302) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2777247) q[0];
sx q[0];
rz(-1.4879328) q[0];
sx q[0];
rz(-0.079282916) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20508619) q[2];
sx q[2];
rz(-3.1383927) q[2];
sx q[2];
rz(0.56006685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.43409882) q[1];
sx q[1];
rz(-1.5710305) q[1];
sx q[1];
rz(-3.1386761) q[1];
rz(-1.4079553) q[3];
sx q[3];
rz(-1.0615088) q[3];
sx q[3];
rz(-2.7247938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87127176) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(1.2284944) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-2.9401949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0495618) q[0];
sx q[0];
rz(-1.9632388) q[0];
sx q[0];
rz(-3.1185246) q[0];
rz(1.6231712) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(-3.1186947) q[2];
sx q[2];
rz(-1.2780634) q[2];
sx q[2];
rz(-2.7414049) q[2];
rz(-1.1313698) q[3];
sx q[3];
rz(-1.4356339) q[3];
sx q[3];
rz(2.1826134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
