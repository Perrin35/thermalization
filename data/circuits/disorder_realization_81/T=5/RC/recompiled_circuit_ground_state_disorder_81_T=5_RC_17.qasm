OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(0.42940816) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(-1.5996999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34719742) q[0];
sx q[0];
rz(-1.7849677) q[0];
sx q[0];
rz(1.2783575) q[0];
x q[1];
rz(-0.23878204) q[2];
sx q[2];
rz(-2.9936643) q[2];
sx q[2];
rz(1.5419568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1356537) q[1];
sx q[1];
rz(-2.7126813) q[1];
sx q[1];
rz(-3.0687544) q[1];
rz(-pi) q[2];
rz(2.7807063) q[3];
sx q[3];
rz(-1.8461707) q[3];
sx q[3];
rz(0.1584681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8371007) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(-0.692918) q[2];
rz(-0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3226586) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-2.0767427) q[0];
rz(2.5191567) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(2.650824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1932632) q[0];
sx q[0];
rz(-1.5800467) q[0];
sx q[0];
rz(-0.023764334) q[0];
rz(-1.347001) q[2];
sx q[2];
rz(-1.9836805) q[2];
sx q[2];
rz(-3.0996291) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1176396) q[1];
sx q[1];
rz(-0.7600541) q[1];
sx q[1];
rz(0.13848409) q[1];
rz(-0.048100483) q[3];
sx q[3];
rz(-1.401618) q[3];
sx q[3];
rz(2.1524803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6530767) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(2.901279) q[2];
rz(-0.49992418) q[3];
sx q[3];
rz(-0.58568716) q[3];
sx q[3];
rz(-1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847454) q[0];
sx q[0];
rz(-0.95504967) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(-2.6495972) q[1];
sx q[1];
rz(-1.0849846) q[1];
sx q[1];
rz(-1.6384151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41427754) q[0];
sx q[0];
rz(-0.87782598) q[0];
sx q[0];
rz(2.5056865) q[0];
rz(-pi) q[1];
rz(-0.12937029) q[2];
sx q[2];
rz(-0.40310848) q[2];
sx q[2];
rz(-2.108824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3635849) q[1];
sx q[1];
rz(-2.6202046) q[1];
sx q[1];
rz(2.2849977) q[1];
rz(-pi) q[2];
rz(3.1206661) q[3];
sx q[3];
rz(-1.4731579) q[3];
sx q[3];
rz(-0.00033631246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4846399) q[2];
sx q[2];
rz(-0.52565614) q[2];
sx q[2];
rz(3.0204115) q[2];
rz(-2.1335404) q[3];
sx q[3];
rz(-2.0262148) q[3];
sx q[3];
rz(1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0693531) q[0];
sx q[0];
rz(-3.1407052) q[0];
sx q[0];
rz(0.51373154) q[0];
rz(-2.1836102) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(-1.4428008) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4289249) q[0];
sx q[0];
rz(-2.6972572) q[0];
sx q[0];
rz(-2.7508468) q[0];
x q[1];
rz(-1.3495693) q[2];
sx q[2];
rz(-2.0906679) q[2];
sx q[2];
rz(1.2017565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.084667) q[1];
sx q[1];
rz(-2.5206893) q[1];
sx q[1];
rz(1.4292745) q[1];
rz(1.4738655) q[3];
sx q[3];
rz(-2.2301794) q[3];
sx q[3];
rz(1.7771378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55502597) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(1.8761934) q[2];
rz(-1.966656) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(1.9780212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555399) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(-1.8170005) q[0];
rz(2.1728204) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(2.103215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0875577) q[0];
sx q[0];
rz(-1.4097347) q[0];
sx q[0];
rz(2.7557082) q[0];
rz(-pi) q[1];
rz(0.81401396) q[2];
sx q[2];
rz(-2.2377439) q[2];
sx q[2];
rz(0.89165724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6250302) q[1];
sx q[1];
rz(-1.6119526) q[1];
sx q[1];
rz(-2.7709318) q[1];
x q[2];
rz(0.65870993) q[3];
sx q[3];
rz(-1.988387) q[3];
sx q[3];
rz(-2.5233248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6394627) q[2];
sx q[2];
rz(-2.837193) q[2];
sx q[2];
rz(-2.2501865) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(2.4752899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0797743) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(-0.99639446) q[0];
rz(0.90006104) q[1];
sx q[1];
rz(-2.3521017) q[1];
sx q[1];
rz(2.9642504) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9422007) q[0];
sx q[0];
rz(-0.15155242) q[0];
sx q[0];
rz(-1.569239) q[0];
rz(-pi) q[1];
rz(1.5992237) q[2];
sx q[2];
rz(-0.6387944) q[2];
sx q[2];
rz(0.5109238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6533493) q[1];
sx q[1];
rz(-2.0470139) q[1];
sx q[1];
rz(0.32250065) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7436036) q[3];
sx q[3];
rz(-2.6263642) q[3];
sx q[3];
rz(2.990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8511054) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(0.30430749) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(2.2296026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07311634) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-2.2090744) q[0];
rz(3.0724691) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(2.1014012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0737586) q[0];
sx q[0];
rz(-0.76601765) q[0];
sx q[0];
rz(0.037184663) q[0];
x q[1];
rz(1.8273152) q[2];
sx q[2];
rz(-1.7396915) q[2];
sx q[2];
rz(-1.0948563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55807796) q[1];
sx q[1];
rz(-0.78418523) q[1];
sx q[1];
rz(0.23434831) q[1];
x q[2];
rz(0.64994855) q[3];
sx q[3];
rz(-0.81226617) q[3];
sx q[3];
rz(-0.60571972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0065464) q[2];
sx q[2];
rz(-2.4032205) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(1.0829571) q[3];
sx q[3];
rz(-1.0561918) q[3];
sx q[3];
rz(2.1347031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873544) q[0];
sx q[0];
rz(-0.08389689) q[0];
sx q[0];
rz(-0.49302897) q[0];
rz(-0.21656491) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(2.1176178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39505491) q[0];
sx q[0];
rz(-0.97486541) q[0];
sx q[0];
rz(1.5563957) q[0];
x q[1];
rz(1.4570974) q[2];
sx q[2];
rz(-1.9010882) q[2];
sx q[2];
rz(-1.8735261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1962016) q[1];
sx q[1];
rz(-1.7111254) q[1];
sx q[1];
rz(0.36044557) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8423457) q[3];
sx q[3];
rz(-0.95029921) q[3];
sx q[3];
rz(2.3503691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0852802) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(-0.67576605) q[2];
rz(-2.7285649) q[3];
sx q[3];
rz(-1.690381) q[3];
sx q[3];
rz(-0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-2.4830723) q[0];
sx q[0];
rz(-3.107048) q[0];
rz(3.0870364) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(0.82130718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217348) q[0];
sx q[0];
rz(-0.99911753) q[0];
sx q[0];
rz(-1.4522533) q[0];
x q[1];
rz(-3.0376833) q[2];
sx q[2];
rz(-3.0523411) q[2];
sx q[2];
rz(-0.98429843) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9266204) q[1];
sx q[1];
rz(-1.4550303) q[1];
sx q[1];
rz(-0.15067071) q[1];
rz(-pi) q[2];
rz(2.2230439) q[3];
sx q[3];
rz(-1.9652742) q[3];
sx q[3];
rz(-0.6834417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36046946) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(0.13295573) q[2];
rz(2.7305056) q[3];
sx q[3];
rz(-1.5186331) q[3];
sx q[3];
rz(-2.0558004) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0947615) q[0];
sx q[0];
rz(-2.5346041) q[0];
sx q[0];
rz(1.8827615) q[0];
rz(1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(2.3283995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55594873) q[0];
sx q[0];
rz(-3.0559982) q[0];
sx q[0];
rz(-0.42285796) q[0];
rz(-3.0243479) q[2];
sx q[2];
rz(-2.239253) q[2];
sx q[2];
rz(0.14456597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.836536) q[1];
sx q[1];
rz(-2.146135) q[1];
sx q[1];
rz(-2.1161377) q[1];
x q[2];
rz(0.2701668) q[3];
sx q[3];
rz(-2.674683) q[3];
sx q[3];
rz(1.8570569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64858156) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(-2.0900334) q[2];
rz(2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(0.65783182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(1.3195994) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(-1.6435087) q[2];
sx q[2];
rz(-1.3617392) q[2];
sx q[2];
rz(-2.9206252) q[2];
rz(-2.9444957) q[3];
sx q[3];
rz(-1.6296248) q[3];
sx q[3];
rz(-1.5546391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
