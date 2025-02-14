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
rz(-0.65161172) q[1];
sx q[1];
rz(-1.7348644) q[1];
sx q[1];
rz(0.23174098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.988294) q[0];
sx q[0];
rz(-0.35044119) q[0];
sx q[0];
rz(0.20163433) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54643537) q[2];
sx q[2];
rz(-1.5218415) q[2];
sx q[2];
rz(-3.0374683) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3331795) q[1];
sx q[1];
rz(-1.1929379) q[1];
sx q[1];
rz(0.66587944) q[1];
rz(-2.8032325) q[3];
sx q[3];
rz(-1.6099811) q[3];
sx q[3];
rz(3.0704263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(1.862662) q[2];
rz(0.88979641) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(2.8029627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(3.056114) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2538214) q[0];
sx q[0];
rz(-1.5740875) q[0];
sx q[0];
rz(-0.38981593) q[0];
rz(2.1824942) q[2];
sx q[2];
rz(-1.8597721) q[2];
sx q[2];
rz(-2.0069569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6798415) q[1];
sx q[1];
rz(-2.0978758) q[1];
sx q[1];
rz(-0.26694571) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68686266) q[3];
sx q[3];
rz(-1.1986889) q[3];
sx q[3];
rz(2.5290979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-0.386664) q[2];
rz(1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(-2.9215422) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(0.49705848) q[0];
rz(-2.0992384) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(2.846948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9884856) q[0];
sx q[0];
rz(-1.9026347) q[0];
sx q[0];
rz(1.6436623) q[0];
rz(-pi) q[1];
rz(-1.7435837) q[2];
sx q[2];
rz(-2.2026988) q[2];
sx q[2];
rz(1.2829219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68874796) q[1];
sx q[1];
rz(-2.0032681) q[1];
sx q[1];
rz(-2.2562863) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9975843) q[3];
sx q[3];
rz(-2.1668808) q[3];
sx q[3];
rz(-2.2726187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3543388) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(1.5617237) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(0.92897433) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(-2.9122747) q[0];
rz(-1.5062821) q[1];
sx q[1];
rz(-1.458026) q[1];
sx q[1];
rz(3.07952) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6939225) q[0];
sx q[0];
rz(-1.0167443) q[0];
sx q[0];
rz(2.5497132) q[0];
x q[1];
rz(-2.1770085) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(3.0678444) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0921539) q[1];
sx q[1];
rz(-2.2518932) q[1];
sx q[1];
rz(-1.9309631) q[1];
x q[2];
rz(2.9406406) q[3];
sx q[3];
rz(-1.6621328) q[3];
sx q[3];
rz(-0.82543711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0486003) q[2];
sx q[2];
rz(-0.91602641) q[2];
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
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.49486092) q[0];
sx q[0];
rz(-2.4202388) q[0];
sx q[0];
rz(-2.2485961) q[0];
rz(-2.112174) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(-3.0457048) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77377787) q[0];
sx q[0];
rz(-0.68499631) q[0];
sx q[0];
rz(-2.2626586) q[0];
rz(-pi) q[1];
rz(0.90135677) q[2];
sx q[2];
rz(-3.0283805) q[2];
sx q[2];
rz(0.33153807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27285114) q[1];
sx q[1];
rz(-0.68528803) q[1];
sx q[1];
rz(-2.5455695) q[1];
rz(-pi) q[2];
rz(2.5378102) q[3];
sx q[3];
rz(-2.5386435) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(-2.7161157) q[2];
rz(-2.7191539) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4130037) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(-0.11716209) q[0];
rz(1.4940184) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-1.0711627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534096) q[0];
sx q[0];
rz(-2.5700535) q[0];
sx q[0];
rz(-2.6009212) q[0];
x q[1];
rz(-1.8495249) q[2];
sx q[2];
rz(-1.0459002) q[2];
sx q[2];
rz(-0.43290813) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90729672) q[1];
sx q[1];
rz(-1.7959372) q[1];
sx q[1];
rz(0.87323685) q[1];
rz(-pi) q[2];
rz(-0.46573205) q[3];
sx q[3];
rz(-2.7890786) q[3];
sx q[3];
rz(2.5172212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57227197) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-3.102109) q[2];
rz(-3.0651921) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(-2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70306784) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(-2.1912498) q[0];
rz(1.1514661) q[1];
sx q[1];
rz(-1.6116424) q[1];
sx q[1];
rz(-1.8340402) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261953) q[0];
sx q[0];
rz(-1.2457677) q[0];
sx q[0];
rz(2.9103878) q[0];
x q[1];
rz(2.9588863) q[2];
sx q[2];
rz(-0.38190834) q[2];
sx q[2];
rz(0.71268247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79297355) q[1];
sx q[1];
rz(-1.7820616) q[1];
sx q[1];
rz(1.068685) q[1];
rz(0.1343822) q[3];
sx q[3];
rz(-0.9690949) q[3];
sx q[3];
rz(2.7190405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-0.68009192) q[2];
sx q[2];
rz(-0.55366984) q[2];
rz(2.4456444) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(-1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(-0.30690646) q[0];
rz(-1.2762997) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.9006405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6825323) q[0];
sx q[0];
rz(-2.0143565) q[0];
sx q[0];
rz(0.17291594) q[0];
rz(-pi) q[1];
rz(0.5979171) q[2];
sx q[2];
rz(-1.9321402) q[2];
sx q[2];
rz(1.8449699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7684264) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(-1.7684446) q[1];
rz(-1.9718902) q[3];
sx q[3];
rz(-0.99906427) q[3];
sx q[3];
rz(0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3500195) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(1.2903068) q[2];
rz(-2.4604515) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(-1.9679029) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(-3.0618844) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3143334) q[0];
sx q[0];
rz(-1.200935) q[0];
sx q[0];
rz(2.8989559) q[0];
x q[1];
rz(2.6775421) q[2];
sx q[2];
rz(-2.8949605) q[2];
sx q[2];
rz(0.76619785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4533055) q[1];
sx q[1];
rz(-2.2548179) q[1];
sx q[1];
rz(-2.3832393) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8675141) q[3];
sx q[3];
rz(-2.35459) q[3];
sx q[3];
rz(-2.3001461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.8076753) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(-2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659336) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(1.851086) q[0];
rz(-0.036529649) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.999015) q[0];
sx q[0];
rz(-1.8685088) q[0];
sx q[0];
rz(-1.1243058) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2453111) q[2];
sx q[2];
rz(-0.41322979) q[2];
sx q[2];
rz(-2.9256224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86856213) q[1];
sx q[1];
rz(-2.2924463) q[1];
sx q[1];
rz(2.9356586) q[1];
x q[2];
rz(1.6855816) q[3];
sx q[3];
rz(-1.8978531) q[3];
sx q[3];
rz(-0.17010526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(2.3850208) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(1.1463317) q[2];
sx q[2];
rz(-2.2590315) q[2];
sx q[2];
rz(-3.0207241) q[2];
rz(2.041009) q[3];
sx q[3];
rz(-2.1775424) q[3];
sx q[3];
rz(1.2829124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
