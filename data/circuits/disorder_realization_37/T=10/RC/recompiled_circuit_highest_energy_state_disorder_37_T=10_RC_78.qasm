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
rz(-1.3353424) q[0];
sx q[0];
rz(-2.775132) q[0];
sx q[0];
rz(-0.45912826) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739959) q[0];
sx q[0];
rz(-1.2277506) q[0];
sx q[0];
rz(-1.497722) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0475869) q[2];
sx q[2];
rz(-0.54840198) q[2];
sx q[2];
rz(1.7552055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3430378) q[1];
sx q[1];
rz(-2.3904368) q[1];
sx q[1];
rz(0.57111528) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11756331) q[3];
sx q[3];
rz(-2.8010578) q[3];
sx q[3];
rz(-1.3887608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(2.8029627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(-0.92581785) q[0];
rz(-2.0766808) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(0.085478641) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31562284) q[0];
sx q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(-1.5743544) q[0];
rz(-pi) q[1];
rz(2.1824942) q[2];
sx q[2];
rz(-1.2818205) q[2];
sx q[2];
rz(2.0069569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.036219941) q[1];
sx q[1];
rz(-0.58507996) q[1];
sx q[1];
rz(-1.9963422) q[1];
rz(2.45473) q[3];
sx q[3];
rz(-1.9429038) q[3];
sx q[3];
rz(-2.5290979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.815879) q[2];
sx q[2];
rz(-1.6625983) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(-1.2169085) q[3];
sx q[3];
rz(-2.7972126) q[3];
sx q[3];
rz(0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3306408) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-0.49705848) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(-2.846948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.153107) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(1.4979304) q[0];
rz(-pi) q[1];
rz(1.7435837) q[2];
sx q[2];
rz(-0.9388939) q[2];
sx q[2];
rz(-1.8586707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5897419) q[1];
sx q[1];
rz(-0.95847469) q[1];
sx q[1];
rz(0.53770868) q[1];
rz(-pi) q[2];
rz(-2.1717242) q[3];
sx q[3];
rz(-1.4517541) q[3];
sx q[3];
rz(0.62059072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(-1.076237) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266066) q[0];
sx q[0];
rz(-1.6153233) q[0];
sx q[0];
rz(-2.9122747) q[0];
rz(-1.5062821) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(-3.07952) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-2.1770085) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(-0.073748253) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.049438795) q[1];
sx q[1];
rz(-2.2518932) q[1];
sx q[1];
rz(1.9309631) q[1];
rz(-pi) q[2];
rz(-1.4775949) q[3];
sx q[3];
rz(-1.7708994) q[3];
sx q[3];
rz(-2.4148108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0486003) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.3108866) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(-2.2485961) q[0];
rz(1.0294186) q[1];
sx q[1];
rz(-0.72975492) q[1];
sx q[1];
rz(3.0457048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22658081) q[0];
sx q[0];
rz(-1.9862729) q[0];
sx q[0];
rz(-2.1323432) q[0];
rz(-pi) q[1];
rz(-3.0711543) q[2];
sx q[2];
rz(-1.6595006) q[2];
sx q[2];
rz(-0.34102893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3600625) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(-2.5468154) q[1];
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
rz(2.9868077) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(2.7191539) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130037) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-2.2393176) q[1];
sx q[1];
rz(2.07043) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18497047) q[0];
sx q[0];
rz(-1.8529467) q[0];
sx q[0];
rz(2.6376702) q[0];
x q[1];
rz(-1.2920678) q[2];
sx q[2];
rz(-2.0956925) q[2];
sx q[2];
rz(-0.43290813) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90729672) q[1];
sx q[1];
rz(-1.3456555) q[1];
sx q[1];
rz(-0.87323685) q[1];
rz(-pi) q[2];
rz(1.4070687) q[3];
sx q[3];
rz(-1.2571954) q[3];
sx q[3];
rz(2.0255476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(3.102109) q[2];
rz(0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(2.6881645) q[3];
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
rz(-pi) q[3];
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
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(1.1514661) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(-1.8340402) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3111585) q[0];
sx q[0];
rz(-1.3518999) q[0];
sx q[0];
rz(-1.9040742) q[0];
rz(-1.4979532) q[2];
sx q[2];
rz(-1.1955639) q[2];
sx q[2];
rz(-2.2323687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.4129909) q[1];
sx q[1];
rz(-0.54122335) q[1];
sx q[1];
rz(1.1515929) q[1];
rz(3.0072104) q[3];
sx q[3];
rz(-2.1724977) q[3];
sx q[3];
rz(-0.4225522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-2.5879228) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(0.30690646) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.9006405) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0690445) q[0];
sx q[0];
rz(-0.47397754) q[0];
sx q[0];
rz(-1.9182253) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5979171) q[2];
sx q[2];
rz(-1.9321402) q[2];
sx q[2];
rz(-1.2966228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13682374) q[1];
sx q[1];
rz(-1.3826177) q[1];
sx q[1];
rz(2.8269563) q[1];
rz(-1.9718902) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3500195) q[2];
sx q[2];
rz(-0.80524033) q[2];
sx q[2];
rz(1.2903068) q[2];
rz(-0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(-2.5545919) q[1];
sx q[1];
rz(-0.78158164) q[1];
sx q[1];
rz(-3.0618844) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3143334) q[0];
sx q[0];
rz(-1.9406576) q[0];
sx q[0];
rz(2.8989559) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46405053) q[2];
sx q[2];
rz(-0.24663217) q[2];
sx q[2];
rz(2.3753948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6708128) q[1];
sx q[1];
rz(-0.97320405) q[1];
sx q[1];
rz(-2.2714991) q[1];
rz(1.8675141) q[3];
sx q[3];
rz(-0.78700262) q[3];
sx q[3];
rz(0.84144652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54032636) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5909083) q[3];
sx q[3];
rz(-1.0100789) q[3];
sx q[3];
rz(2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(-3.105063) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1636476) q[0];
sx q[0];
rz(-0.53102101) q[0];
sx q[0];
rz(-2.1885896) q[0];
rz(2.8743083) q[2];
sx q[2];
rz(-1.8898095) q[2];
sx q[2];
rz(0.50179447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2730305) q[1];
sx q[1];
rz(-0.84914637) q[1];
sx q[1];
rz(2.9356586) q[1];
x q[2];
rz(-0.32907069) q[3];
sx q[3];
rz(-1.6794723) q[3];
sx q[3];
rz(1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9145987) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(2.3850208) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(-2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.666438) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-2.4074211) q[2];
sx q[2];
rz(-1.894507) q[2];
sx q[2];
rz(1.4121216) q[2];
rz(0.57831709) q[3];
sx q[3];
rz(-0.74902799) q[3];
sx q[3];
rz(2.0109162) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
