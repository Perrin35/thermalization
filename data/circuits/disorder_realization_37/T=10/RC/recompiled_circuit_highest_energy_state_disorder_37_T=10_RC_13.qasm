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
rz(1.8062502) q[0];
sx q[0];
rz(-0.36646068) q[0];
sx q[0];
rz(0.45912826) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2278175) q[0];
sx q[0];
rz(-1.5019866) q[0];
sx q[0];
rz(-2.7976996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54643537) q[2];
sx q[2];
rz(-1.6197512) q[2];
sx q[2];
rz(-3.0374683) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3430378) q[1];
sx q[1];
rz(-0.75115582) q[1];
sx q[1];
rz(2.5704774) q[1];
rz(1.6123338) q[3];
sx q[3];
rz(-1.9088863) q[3];
sx q[3];
rz(-1.5134144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.090652466) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(0.92581785) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(0.085478641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31562284) q[0];
sx q[0];
rz(-1.96061) q[0];
sx q[0];
rz(-1.5743544) q[0];
x q[1];
rz(2.7932554) q[2];
sx q[2];
rz(-2.1536963) q[2];
sx q[2];
rz(-0.23886853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6798415) q[1];
sx q[1];
rz(-1.0437168) q[1];
sx q[1];
rz(-2.8746469) q[1];
rz(-pi) q[2];
x q[2];
rz(2.45473) q[3];
sx q[3];
rz(-1.9429038) q[3];
sx q[3];
rz(-2.5290979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(0.386664) q[2];
rz(-1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306408) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(0.49705848) q[0];
rz(-1.0423543) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(2.846948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7001273) q[0];
sx q[0];
rz(-1.5019121) q[0];
sx q[0];
rz(-2.8089351) q[0];
rz(0.23068409) q[2];
sx q[2];
rz(-2.4896224) q[2];
sx q[2];
rz(0.99562746) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55185071) q[1];
sx q[1];
rz(-2.183118) q[1];
sx q[1];
rz(0.53770868) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96986846) q[3];
sx q[3];
rz(-1.4517541) q[3];
sx q[3];
rz(-2.5210019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7872539) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(1.5798689) q[2];
rz(2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314986) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(0.22931799) q[0];
rz(1.5062821) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(-3.07952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4631043) q[0];
sx q[0];
rz(-2.0651428) q[0];
sx q[0];
rz(2.2114179) q[0];
rz(-1.1787291) q[2];
sx q[2];
rz(-1.3118366) q[2];
sx q[2];
rz(0.94253892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.58846274) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(-2.7314145) q[1];
rz(2.7113879) q[3];
sx q[3];
rz(-2.9211126) q[3];
sx q[3];
rz(-1.9752432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.092992358) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(1.830706) q[2];
rz(-3.1033031) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49486092) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(-0.89299655) q[0];
rz(1.0294186) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(3.0457048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77377787) q[0];
sx q[0];
rz(-0.68499631) q[0];
sx q[0];
rz(-2.2626586) q[0];
x q[1];
rz(-1.65972) q[2];
sx q[2];
rz(-1.5006353) q[2];
sx q[2];
rz(1.905575) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.694987) q[1];
sx q[1];
rz(-1.0195273) q[1];
sx q[1];
rz(2.0010082) q[1];
x q[2];
rz(-2.5378102) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15478495) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(2.7191539) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4130037) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(-1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(2.07043) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2334918) q[0];
sx q[0];
rz(-2.0530434) q[0];
sx q[0];
rz(1.890475) q[0];
rz(-1.8495249) q[2];
sx q[2];
rz(-2.0956925) q[2];
sx q[2];
rz(0.43290813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90729672) q[1];
sx q[1];
rz(-1.3456555) q[1];
sx q[1];
rz(2.2683558) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46573205) q[3];
sx q[3];
rz(-0.35251401) q[3];
sx q[3];
rz(-2.5172212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5693207) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(0.039483698) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(-2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(0.95034289) q[0];
rz(1.9901265) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.8340402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261953) q[0];
sx q[0];
rz(-1.2457677) q[0];
sx q[0];
rz(-0.23120489) q[0];
rz(-2.9588863) q[2];
sx q[2];
rz(-2.7596843) q[2];
sx q[2];
rz(0.71268247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79297355) q[1];
sx q[1];
rz(-1.7820616) q[1];
sx q[1];
rz(2.0729077) q[1];
rz(1.3780955) q[3];
sx q[3];
rz(-0.61470882) q[3];
sx q[3];
rz(-0.65700442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(0.55366984) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(-1.455201) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(2.8346862) q[0];
rz(-1.2762997) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(-1.9006405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6825323) q[0];
sx q[0];
rz(-1.1272361) q[0];
sx q[0];
rz(-2.9686767) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9996793) q[2];
sx q[2];
rz(-2.1253573) q[2];
sx q[2];
rz(-0.51046023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0047689) q[1];
sx q[1];
rz(-1.3826177) q[1];
sx q[1];
rz(2.8269563) q[1];
x q[2];
rz(1.9718902) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-2.9740262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3500195) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(-2.4604515) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(1.6397938) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660368) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(2.5545919) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(3.0618844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82725924) q[0];
sx q[0];
rz(-1.200935) q[0];
sx q[0];
rz(-0.24263675) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9201512) q[2];
sx q[2];
rz(-1.6802854) q[2];
sx q[2];
rz(1.2564645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4533055) q[1];
sx q[1];
rz(-2.2548179) q[1];
sx q[1];
rz(0.75835336) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2853197) q[3];
sx q[3];
rz(-0.82672182) q[3];
sx q[3];
rz(0.43275012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54032636) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(-0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.2905066) q[0];
rz(-3.105063) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.288702) q[0];
sx q[0];
rz(-1.1452617) q[0];
sx q[0];
rz(2.8136926) q[0];
x q[1];
rz(-1.9007334) q[2];
sx q[2];
rz(-1.82429) q[2];
sx q[2];
rz(1.1546749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.579291) q[1];
sx q[1];
rz(-0.74534432) q[1];
sx q[1];
rz(-1.3424804) q[1];
x q[2];
rz(-2.815991) q[3];
sx q[3];
rz(-0.34593098) q[3];
sx q[3];
rz(-2.9671362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2269939) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.8168137) q[2];
rz(2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.666438) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(-1.7186164) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-1.995261) q[2];
sx q[2];
rz(-2.2590315) q[2];
sx q[2];
rz(-3.0207241) q[2];
rz(1.1005836) q[3];
sx q[3];
rz(-0.96405021) q[3];
sx q[3];
rz(-1.8586803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
