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
rz(-1.1779397) q[0];
sx q[0];
rz(-3.0484634) q[0];
sx q[0];
rz(2.4094474) q[0];
rz(1.572345) q[1];
sx q[1];
rz(0.88580004) q[1];
sx q[1];
rz(9.8635397) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51073815) q[0];
sx q[0];
rz(-0.33677142) q[0];
sx q[0];
rz(-1.4022409) q[0];
rz(-1.8541502) q[2];
sx q[2];
rz(-1.3204976) q[2];
sx q[2];
rz(1.3586501) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10459929) q[1];
sx q[1];
rz(-1.5244836) q[1];
sx q[1];
rz(-2.8830257) q[1];
rz(-pi) q[2];
rz(-3.0436727) q[3];
sx q[3];
rz(-0.46028462) q[3];
sx q[3];
rz(1.8079881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0977122) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(0.91862339) q[2];
rz(1.0527323) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(-0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.265428) q[0];
sx q[0];
rz(-2.3227203) q[0];
sx q[0];
rz(2.2692666) q[0];
rz(2.1288952) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-2.8025119) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10386183) q[0];
sx q[0];
rz(-1.4418238) q[0];
sx q[0];
rz(-0.43855389) q[0];
x q[1];
rz(-2.2834899) q[2];
sx q[2];
rz(-0.65026186) q[2];
sx q[2];
rz(1.5257143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6244858) q[1];
sx q[1];
rz(-3.0992866) q[1];
sx q[1];
rz(-2.9632764) q[1];
rz(-pi) q[2];
rz(-1.2148772) q[3];
sx q[3];
rz(-2.7971075) q[3];
sx q[3];
rz(-3.0657299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33112153) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(2.5362711) q[2];
rz(-2.807054) q[3];
sx q[3];
rz(-2.7693222) q[3];
sx q[3];
rz(3.065006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89638585) q[0];
sx q[0];
rz(-2.4672282) q[0];
sx q[0];
rz(-1.6687923) q[0];
rz(0.32707602) q[1];
sx q[1];
rz(-1.8388803) q[1];
sx q[1];
rz(-0.27007857) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03878394) q[0];
sx q[0];
rz(-2.5804418) q[0];
sx q[0];
rz(0.3766685) q[0];
x q[1];
rz(-0.075614838) q[2];
sx q[2];
rz(-1.9984198) q[2];
sx q[2];
rz(2.2958034) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1317989) q[1];
sx q[1];
rz(-2.0094618) q[1];
sx q[1];
rz(-0.75097221) q[1];
rz(1.5669426) q[3];
sx q[3];
rz(-1.4906297) q[3];
sx q[3];
rz(-1.4269258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.104287) q[2];
sx q[2];
rz(-1.1423926) q[2];
sx q[2];
rz(-1.8541065) q[2];
rz(-1.2152952) q[3];
sx q[3];
rz(-1.7561965) q[3];
sx q[3];
rz(2.9638885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52008587) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(0.31975123) q[0];
rz(0.41140914) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(-3.06126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31422752) q[0];
sx q[0];
rz(-0.36970678) q[0];
sx q[0];
rz(0.99401919) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14239399) q[2];
sx q[2];
rz(-2.3300397) q[2];
sx q[2];
rz(0.55003563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4027202) q[1];
sx q[1];
rz(-2.3490073) q[1];
sx q[1];
rz(-1.9566105) q[1];
rz(0.37253054) q[3];
sx q[3];
rz(-0.64465678) q[3];
sx q[3];
rz(2.9304867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8722998) q[2];
sx q[2];
rz(-2.9041957) q[2];
sx q[2];
rz(-0.044895127) q[2];
rz(2.6063555) q[3];
sx q[3];
rz(-1.6343296) q[3];
sx q[3];
rz(-3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1503247) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(-0.87565652) q[0];
rz(-0.21370299) q[1];
sx q[1];
rz(-0.49476606) q[1];
sx q[1];
rz(1.5692086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039883651) q[0];
sx q[0];
rz(-1.3613762) q[0];
sx q[0];
rz(-1.8906192) q[0];
rz(-pi) q[1];
rz(1.3415229) q[2];
sx q[2];
rz(-1.9976282) q[2];
sx q[2];
rz(-0.62139702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8991291) q[1];
sx q[1];
rz(-1.6751667) q[1];
sx q[1];
rz(-0.69931286) q[1];
rz(0.32099681) q[3];
sx q[3];
rz(-1.4813652) q[3];
sx q[3];
rz(-0.018370779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.956942) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(-1.0682028) q[2];
rz(-2.5351561) q[3];
sx q[3];
rz(-1.9304099) q[3];
sx q[3];
rz(-3.1312805) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(-1.3909719) q[0];
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
rz(-1.8075631) q[0];
sx q[0];
rz(-1.5153756) q[0];
sx q[0];
rz(1.2822654) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.053534075) q[2];
sx q[2];
rz(-1.495289) q[2];
sx q[2];
rz(-1.8078505) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5038915) q[1];
sx q[1];
rz(-0.63234416) q[1];
sx q[1];
rz(-1.31193) q[1];
rz(-0.88504412) q[3];
sx q[3];
rz(-0.17286271) q[3];
sx q[3];
rz(0.96641738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4411053) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(-2.109745) q[2];
rz(1.4619689) q[3];
sx q[3];
rz(-0.69457355) q[3];
sx q[3];
rz(1.3612548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66330355) q[0];
sx q[0];
rz(-0.42600584) q[0];
sx q[0];
rz(-1.7167094) q[0];
rz(-2.5583963) q[1];
sx q[1];
rz(-1.2251264) q[1];
sx q[1];
rz(-3.084175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9641477) q[0];
sx q[0];
rz(-2.3741407) q[0];
sx q[0];
rz(1.5758118) q[0];
rz(-pi) q[1];
rz(-1.0591828) q[2];
sx q[2];
rz(-2.1202737) q[2];
sx q[2];
rz(1.8715931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0432644) q[1];
sx q[1];
rz(-0.66126138) q[1];
sx q[1];
rz(-2.8533881) q[1];
rz(-pi) q[2];
x q[2];
rz(2.613861) q[3];
sx q[3];
rz(-2.4411628) q[3];
sx q[3];
rz(-0.34175261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8251557) q[2];
sx q[2];
rz(-1.8437443) q[2];
sx q[2];
rz(1.3081029) q[2];
rz(1.3206652) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(0.86110419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945781) q[0];
sx q[0];
rz(-0.68462831) q[0];
sx q[0];
rz(-1.8735877) q[0];
rz(2.0199203) q[1];
sx q[1];
rz(-1.8662165) q[1];
sx q[1];
rz(2.4915288) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92597678) q[0];
sx q[0];
rz(-2.4575666) q[0];
sx q[0];
rz(0.90606545) q[0];
rz(-0.68393884) q[2];
sx q[2];
rz(-1.142148) q[2];
sx q[2];
rz(-2.6446083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6942581) q[1];
sx q[1];
rz(-1.9617394) q[1];
sx q[1];
rz(2.5503134) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88553377) q[3];
sx q[3];
rz(-0.69708744) q[3];
sx q[3];
rz(-0.05230418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.140363) q[2];
sx q[2];
rz(-1.659212) q[2];
sx q[2];
rz(-0.39230997) q[2];
rz(-0.20491925) q[3];
sx q[3];
rz(-0.50060087) q[3];
sx q[3];
rz(0.71995455) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0023163) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(1.4578777) q[0];
rz(-2.814759) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(-2.0976417) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69125873) q[0];
sx q[0];
rz(-0.46713167) q[0];
sx q[0];
rz(2.0005666) q[0];
rz(-pi) q[1];
rz(1.4307666) q[2];
sx q[2];
rz(-1.1159889) q[2];
sx q[2];
rz(0.32250139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9324016) q[1];
sx q[1];
rz(-1.652601) q[1];
sx q[1];
rz(1.1751168) q[1];
rz(2.8159385) q[3];
sx q[3];
rz(-0.94653385) q[3];
sx q[3];
rz(0.52407085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30847183) q[2];
sx q[2];
rz(-2.1127508) q[2];
sx q[2];
rz(2.2157045) q[2];
rz(-2.0738257) q[3];
sx q[3];
rz(-1.1238778) q[3];
sx q[3];
rz(-1.7759751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866216) q[0];
sx q[0];
rz(-1.4737031) q[0];
sx q[0];
rz(2.5545252) q[0];
rz(0.58285561) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(-0.48097441) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0635446) q[0];
sx q[0];
rz(-0.46351156) q[0];
sx q[0];
rz(-2.6592451) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9312385) q[2];
sx q[2];
rz(-1.3925465) q[2];
sx q[2];
rz(-1.0457863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1744057) q[1];
sx q[1];
rz(-2.4250826) q[1];
sx q[1];
rz(0.67542507) q[1];
x q[2];
rz(-1.3708047) q[3];
sx q[3];
rz(-1.963435) q[3];
sx q[3];
rz(2.5225366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9440072) q[2];
sx q[2];
rz(-0.53637594) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(-2.3594989) q[3];
sx q[3];
rz(-2.8160281) q[3];
sx q[3];
rz(1.1480924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999577) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(1.5326473) q[1];
sx q[1];
rz(-1.9736704) q[1];
sx q[1];
rz(-0.71221487) q[1];
rz(-0.73405042) q[2];
sx q[2];
rz(-1.0577591) q[2];
sx q[2];
rz(-1.9815097) q[2];
rz(-1.731012) q[3];
sx q[3];
rz(-2.0638777) q[3];
sx q[3];
rz(-0.0085245098) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
