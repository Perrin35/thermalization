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
rz(-1.4425059) q[0];
sx q[0];
rz(-1.2966172) q[0];
sx q[0];
rz(1.698864) q[0];
rz(-0.35248414) q[1];
sx q[1];
rz(-0.52630693) q[1];
sx q[1];
rz(1.7736645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57940021) q[0];
sx q[0];
rz(-0.44792563) q[0];
sx q[0];
rz(0.62904398) q[0];
rz(-pi) q[1];
rz(0.10124293) q[2];
sx q[2];
rz(-0.35738073) q[2];
sx q[2];
rz(-1.9592154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96208159) q[1];
sx q[1];
rz(-1.700811) q[1];
sx q[1];
rz(0.52501734) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62217766) q[3];
sx q[3];
rz(-2.2832979) q[3];
sx q[3];
rz(-1.5714558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0679396) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(1.2037207) q[2];
rz(-2.1667571) q[3];
sx q[3];
rz(-2.1741512) q[3];
sx q[3];
rz(1.8956641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-0.75342733) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(-0.88428307) q[0];
rz(0.70941225) q[1];
sx q[1];
rz(-1.8953036) q[1];
sx q[1];
rz(2.2394004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35036885) q[0];
sx q[0];
rz(-1.1205427) q[0];
sx q[0];
rz(1.7993991) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96376597) q[2];
sx q[2];
rz(-1.4861408) q[2];
sx q[2];
rz(-2.7805258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98595218) q[1];
sx q[1];
rz(-1.8690171) q[1];
sx q[1];
rz(-1.8267426) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0445116) q[3];
sx q[3];
rz(-0.15972129) q[3];
sx q[3];
rz(2.2217939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.117729) q[2];
sx q[2];
rz(-2.8459097) q[2];
sx q[2];
rz(1.4884865) q[2];
rz(1.9272517) q[3];
sx q[3];
rz(-1.4597273) q[3];
sx q[3];
rz(-1.2795718) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073256) q[0];
sx q[0];
rz(-1.7420344) q[0];
sx q[0];
rz(-2.9752327) q[0];
rz(-2.4436489) q[1];
sx q[1];
rz(-0.78015399) q[1];
sx q[1];
rz(1.4770329) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1082527) q[0];
sx q[0];
rz(-2.0128002) q[0];
sx q[0];
rz(-3.1382794) q[0];
x q[1];
rz(1.2899621) q[2];
sx q[2];
rz(-1.3094433) q[2];
sx q[2];
rz(0.42119831) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.229847) q[1];
sx q[1];
rz(-1.032879) q[1];
sx q[1];
rz(2.6078546) q[1];
rz(-pi) q[2];
rz(0.03607492) q[3];
sx q[3];
rz(-0.49866184) q[3];
sx q[3];
rz(-0.3968249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4878238) q[2];
sx q[2];
rz(-1.9705801) q[2];
sx q[2];
rz(-2.5269395) q[2];
rz(1.6807618) q[3];
sx q[3];
rz(-1.5771644) q[3];
sx q[3];
rz(3.0998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(1.5186658) q[0];
rz(-0.80865771) q[1];
sx q[1];
rz(-1.2963908) q[1];
sx q[1];
rz(0.048695806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0322389) q[0];
sx q[0];
rz(-1.236602) q[0];
sx q[0];
rz(-0.70698694) q[0];
rz(-pi) q[1];
rz(-2.8559309) q[2];
sx q[2];
rz(-2.7451519) q[2];
sx q[2];
rz(0.91914058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.16597) q[1];
sx q[1];
rz(-2.4104154) q[1];
sx q[1];
rz(-2.380409) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1761769) q[3];
sx q[3];
rz(-1.2413238) q[3];
sx q[3];
rz(-2.5191865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0398756) q[2];
sx q[2];
rz(-2.3287435) q[2];
sx q[2];
rz(1.2561049) q[2];
rz(0.060955437) q[3];
sx q[3];
rz(-0.90164369) q[3];
sx q[3];
rz(-0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561279) q[0];
sx q[0];
rz(-2.0138854) q[0];
sx q[0];
rz(-0.10966478) q[0];
rz(-1.0288641) q[1];
sx q[1];
rz(-1.2867462) q[1];
sx q[1];
rz(1.5922155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76413918) q[0];
sx q[0];
rz(-1.2887338) q[0];
sx q[0];
rz(-2.6741323) q[0];
rz(2.0711871) q[2];
sx q[2];
rz(-1.2493881) q[2];
sx q[2];
rz(-1.6381581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60051892) q[1];
sx q[1];
rz(-1.1364577) q[1];
sx q[1];
rz(0.51323607) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.034329359) q[3];
sx q[3];
rz(-1.5331556) q[3];
sx q[3];
rz(-0.08480367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5515543) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(2.6596098) q[2];
rz(0.41989741) q[3];
sx q[3];
rz(-2.0764543) q[3];
sx q[3];
rz(-0.70986748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14966203) q[0];
sx q[0];
rz(-2.2580632) q[0];
sx q[0];
rz(-0.71257198) q[0];
rz(-0.86992162) q[1];
sx q[1];
rz(-2.7674119) q[1];
sx q[1];
rz(-3.1178927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4159235) q[0];
sx q[0];
rz(-2.5953889) q[0];
sx q[0];
rz(-2.8234981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4392088) q[2];
sx q[2];
rz(-1.1105892) q[2];
sx q[2];
rz(0.61379877) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5620251) q[1];
sx q[1];
rz(-1.7365841) q[1];
sx q[1];
rz(0.001494249) q[1];
rz(-0.48051254) q[3];
sx q[3];
rz(-1.3837722) q[3];
sx q[3];
rz(0.20166892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80453834) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(1.1402593) q[2];
rz(-1.1537457) q[3];
sx q[3];
rz(-1.2146344) q[3];
sx q[3];
rz(-1.3212737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7158647) q[0];
sx q[0];
rz(-1.7540997) q[0];
sx q[0];
rz(-2.7737889) q[0];
rz(1.6416719) q[1];
sx q[1];
rz(-2.2190614) q[1];
sx q[1];
rz(-2.0466764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85830583) q[0];
sx q[0];
rz(-2.0525161) q[0];
sx q[0];
rz(2.879309) q[0];
x q[1];
rz(-0.55042437) q[2];
sx q[2];
rz(-1.4604521) q[2];
sx q[2];
rz(-2.7743024) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8290295) q[1];
sx q[1];
rz(-0.2364279) q[1];
sx q[1];
rz(-0.076024012) q[1];
rz(-pi) q[2];
rz(1.9743552) q[3];
sx q[3];
rz(-1.5903683) q[3];
sx q[3];
rz(0.38364601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0423923) q[2];
sx q[2];
rz(-1.3513869) q[2];
sx q[2];
rz(3.0858827) q[2];
rz(-2.4646711) q[3];
sx q[3];
rz(-2.5261295) q[3];
sx q[3];
rz(-0.87853471) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95213503) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(0.96543717) q[0];
rz(1.8769439) q[1];
sx q[1];
rz(-0.74600428) q[1];
sx q[1];
rz(2.8146578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069939216) q[0];
sx q[0];
rz(-2.7732458) q[0];
sx q[0];
rz(-1.8828431) q[0];
x q[1];
rz(-3.1189402) q[2];
sx q[2];
rz(-2.6955397) q[2];
sx q[2];
rz(1.866801) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85046834) q[1];
sx q[1];
rz(-1.7479436) q[1];
sx q[1];
rz(1.7825104) q[1];
rz(-pi) q[2];
rz(1.502874) q[3];
sx q[3];
rz(-2.0188031) q[3];
sx q[3];
rz(-0.94326708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38291976) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(-0.79603377) q[2];
rz(-3.054079) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(0.93761888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336808) q[0];
sx q[0];
rz(-1.3734564) q[0];
sx q[0];
rz(-0.32980907) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-0.70711702) q[1];
sx q[1];
rz(1.0940394) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1739376) q[0];
sx q[0];
rz(-0.47220818) q[0];
sx q[0];
rz(-0.081993563) q[0];
rz(-pi) q[1];
rz(3.0743672) q[2];
sx q[2];
rz(-1.1738699) q[2];
sx q[2];
rz(2.0527184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3628914) q[1];
sx q[1];
rz(-0.56462949) q[1];
sx q[1];
rz(-2.8933011) q[1];
x q[2];
rz(-1.1393896) q[3];
sx q[3];
rz(-1.6069222) q[3];
sx q[3];
rz(2.2808035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0426992) q[2];
sx q[2];
rz(-0.77184474) q[2];
sx q[2];
rz(-1.137255) q[2];
rz(-0.62816652) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(-2.1342733) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4758509) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(-1.092859) q[0];
rz(1.2650371) q[1];
sx q[1];
rz(-1.3053514) q[1];
sx q[1];
rz(1.0337894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40202478) q[0];
sx q[0];
rz(-1.0453859) q[0];
sx q[0];
rz(0.55243405) q[0];
rz(-pi) q[1];
rz(-1.8142419) q[2];
sx q[2];
rz(-2.3156791) q[2];
sx q[2];
rz(-1.5559529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4235735) q[1];
sx q[1];
rz(-1.1597654) q[1];
sx q[1];
rz(0.6842821) q[1];
rz(-pi) q[2];
rz(1.4321253) q[3];
sx q[3];
rz(-1.6059173) q[3];
sx q[3];
rz(2.7447774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.25896245) q[2];
sx q[2];
rz(-0.7236824) q[2];
sx q[2];
rz(-0.21931973) q[2];
rz(-0.80260459) q[3];
sx q[3];
rz(-1.31253) q[3];
sx q[3];
rz(0.32138166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6312859) q[0];
sx q[0];
rz(-1.8230556) q[0];
sx q[0];
rz(-2.0429116) q[0];
rz(-0.17768271) q[1];
sx q[1];
rz(-1.0164574) q[1];
sx q[1];
rz(-0.16998092) q[1];
rz(-1.646252) q[2];
sx q[2];
rz(-2.4377078) q[2];
sx q[2];
rz(-1.0000142) q[2];
rz(-0.11832506) q[3];
sx q[3];
rz(-1.7601624) q[3];
sx q[3];
rz(2.7047529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
