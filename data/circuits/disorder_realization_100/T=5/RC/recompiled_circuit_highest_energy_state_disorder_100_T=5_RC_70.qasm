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
rz(1.572345) q[1];
sx q[1];
rz(0.88580004) q[1];
sx q[1];
rz(9.8635397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8092361) q[0];
sx q[0];
rz(-1.9026103) q[0];
sx q[0];
rz(-0.058666243) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26023999) q[2];
sx q[2];
rz(-1.2965045) q[2];
sx q[2];
rz(-3.0014474) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6631515) q[1];
sx q[1];
rz(-1.3125129) q[1];
sx q[1];
rz(1.5228935) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6832163) q[3];
sx q[3];
rz(-1.6142369) q[3];
sx q[3];
rz(0.14940748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0977122) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(0.91862339) q[2];
rz(-2.0888603) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(2.2241433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265428) q[0];
sx q[0];
rz(-2.3227203) q[0];
sx q[0];
rz(-2.2692666) q[0];
rz(-2.1288952) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-0.33908078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.942303) q[0];
sx q[0];
rz(-2.6856519) q[0];
sx q[0];
rz(-2.845167) q[0];
x q[1];
rz(-2.6800687) q[2];
sx q[2];
rz(-1.095003) q[2];
sx q[2];
rz(2.3522289) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51710684) q[1];
sx q[1];
rz(-0.042306012) q[1];
sx q[1];
rz(-2.9632764) q[1];
rz(-pi) q[2];
rz(3.0172164) q[3];
sx q[3];
rz(-1.2487097) q[3];
sx q[3];
rz(-0.45201221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8104711) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(0.60532153) q[2];
rz(-0.33453861) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(-3.065006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2452068) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(1.6687923) q[0];
rz(0.32707602) q[1];
sx q[1];
rz(-1.8388803) q[1];
sx q[1];
rz(2.8715141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1028087) q[0];
sx q[0];
rz(-2.5804418) q[0];
sx q[0];
rz(2.7649242) q[0];
x q[1];
rz(0.075614838) q[2];
sx q[2];
rz(-1.9984198) q[2];
sx q[2];
rz(0.84578925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0097937219) q[1];
sx q[1];
rz(-1.1321308) q[1];
sx q[1];
rz(2.3906204) q[1];
x q[2];
rz(1.5669426) q[3];
sx q[3];
rz(-1.4906297) q[3];
sx q[3];
rz(1.7146669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0373056) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(-1.8541065) q[2];
rz(-1.2152952) q[3];
sx q[3];
rz(-1.7561965) q[3];
sx q[3];
rz(-0.17770411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6215068) q[0];
sx q[0];
rz(-2.6984213) q[0];
sx q[0];
rz(2.8218414) q[0];
rz(0.41140914) q[1];
sx q[1];
rz(-1.5965867) q[1];
sx q[1];
rz(-0.080332669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7112996) q[0];
sx q[0];
rz(-1.3724494) q[0];
sx q[0];
rz(-1.8848778) q[0];
x q[1];
rz(-1.4223584) q[2];
sx q[2];
rz(-0.76984902) q[2];
sx q[2];
rz(2.7968873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87819609) q[1];
sx q[1];
rz(-2.2913763) q[1];
sx q[1];
rz(0.36468585) q[1];
x q[2];
rz(-0.37253054) q[3];
sx q[3];
rz(-2.4969359) q[3];
sx q[3];
rz(-0.21110591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2692928) q[2];
sx q[2];
rz(-2.9041957) q[2];
sx q[2];
rz(-0.044895127) q[2];
rz(-0.53523713) q[3];
sx q[3];
rz(-1.6343296) q[3];
sx q[3];
rz(-3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.991268) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(-2.2659361) q[0];
rz(2.9278897) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(1.572384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.101709) q[0];
sx q[0];
rz(-1.7802165) q[0];
sx q[0];
rz(-1.2509734) q[0];
x q[1];
rz(-2.7046811) q[2];
sx q[2];
rz(-1.3624117) q[2];
sx q[2];
rz(2.2885099) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.24092785) q[1];
sx q[1];
rz(-0.87605175) q[1];
sx q[1];
rz(1.7068295) q[1];
x q[2];
rz(2.8205958) q[3];
sx q[3];
rz(-1.6602274) q[3];
sx q[3];
rz(-0.018370779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.956942) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(-1.0682028) q[2];
rz(-2.5351561) q[3];
sx q[3];
rz(-1.2111827) q[3];
sx q[3];
rz(3.1312805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73685) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(-1.7506208) q[0];
rz(2.6178316) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(-1.1837122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340296) q[0];
sx q[0];
rz(-1.5153756) q[0];
sx q[0];
rz(-1.2822654) q[0];
x q[1];
rz(0.053534075) q[2];
sx q[2];
rz(-1.495289) q[2];
sx q[2];
rz(1.8078505) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9549) q[1];
sx q[1];
rz(-0.96265332) q[1];
sx q[1];
rz(-2.9561839) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88504412) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(2.1751753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70048731) q[2];
sx q[2];
rz(-0.59868559) q[2];
sx q[2];
rz(-2.109745) q[2];
rz(-1.4619689) q[3];
sx q[3];
rz(-2.4470191) q[3];
sx q[3];
rz(-1.7803378) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66330355) q[0];
sx q[0];
rz(-0.42600584) q[0];
sx q[0];
rz(-1.4248832) q[0];
rz(-0.58319631) q[1];
sx q[1];
rz(-1.2251264) q[1];
sx q[1];
rz(3.084175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9641477) q[0];
sx q[0];
rz(-2.3741407) q[0];
sx q[0];
rz(-1.5657809) q[0];
rz(0.67441794) q[2];
sx q[2];
rz(-0.73235529) q[2];
sx q[2];
rz(-0.44841097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8805931) q[1];
sx q[1];
rz(-0.94124244) q[1];
sx q[1];
rz(1.7884607) q[1];
x q[2];
rz(1.1693277) q[3];
sx q[3];
rz(-0.9801995) q[3];
sx q[3];
rz(2.1484753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8251557) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(1.3081029) q[2];
rz(-1.8209275) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(-2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84701455) q[0];
sx q[0];
rz(-0.68462831) q[0];
sx q[0];
rz(-1.8735877) q[0];
rz(2.0199203) q[1];
sx q[1];
rz(-1.2753762) q[1];
sx q[1];
rz(0.65006382) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0065816) q[0];
sx q[0];
rz(-1.0502292) q[0];
sx q[0];
rz(-2.6755946) q[0];
rz(-2.103527) q[2];
sx q[2];
rz(-2.1830171) q[2];
sx q[2];
rz(2.39447) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44733455) q[1];
sx q[1];
rz(-1.1798533) q[1];
sx q[1];
rz(-0.59127929) q[1];
x q[2];
rz(-2.6542957) q[3];
sx q[3];
rz(-1.0505884) q[3];
sx q[3];
rz(0.76508026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.140363) q[2];
sx q[2];
rz(-1.659212) q[2];
sx q[2];
rz(2.7492827) q[2];
rz(2.9366734) q[3];
sx q[3];
rz(-0.50060087) q[3];
sx q[3];
rz(-2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1392764) q[0];
sx q[0];
rz(-1.5895695) q[0];
sx q[0];
rz(-1.683715) q[0];
rz(-0.32683364) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(2.0976417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69125873) q[0];
sx q[0];
rz(-2.674461) q[0];
sx q[0];
rz(-1.1410261) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7108261) q[2];
sx q[2];
rz(-1.1159889) q[2];
sx q[2];
rz(0.32250139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20919101) q[1];
sx q[1];
rz(-1.652601) q[1];
sx q[1];
rz(-1.9664759) q[1];
rz(-pi) q[2];
rz(-2.8159385) q[3];
sx q[3];
rz(-2.1950588) q[3];
sx q[3];
rz(-2.6175218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8331208) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(-2.2157045) q[2];
rz(-2.0738257) q[3];
sx q[3];
rz(-1.1238778) q[3];
sx q[3];
rz(1.3656176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2866216) q[0];
sx q[0];
rz(-1.6678896) q[0];
sx q[0];
rz(2.5545252) q[0];
rz(-0.58285561) q[1];
sx q[1];
rz(-2.8522377) q[1];
sx q[1];
rz(-0.48097441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60763393) q[0];
sx q[0];
rz(-1.9780428) q[0];
sx q[0];
rz(1.7986222) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19020657) q[2];
sx q[2];
rz(-1.9252732) q[2];
sx q[2];
rz(0.5917393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9945337) q[1];
sx q[1];
rz(-1.9939341) q[1];
sx q[1];
rz(0.59696859) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7417861) q[3];
sx q[3];
rz(-1.7553864) q[3];
sx q[3];
rz(2.1124482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9440072) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(2.3594989) q[3];
sx q[3];
rz(-2.8160281) q[3];
sx q[3];
rz(-1.1480924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
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
rz(0.73405042) q[2];
sx q[2];
rz(-2.0838336) q[2];
sx q[2];
rz(1.1600829) q[2];
rz(2.6431177) q[3];
sx q[3];
rz(-1.4298021) q[3];
sx q[3];
rz(1.485928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
