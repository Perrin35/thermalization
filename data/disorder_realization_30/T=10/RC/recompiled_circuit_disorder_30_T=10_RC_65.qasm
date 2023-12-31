OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40987504) q[0];
sx q[0];
rz(-0.048225064) q[0];
sx q[0];
rz(-1.6884786) q[0];
rz(-pi) q[1];
rz(0.33485246) q[2];
sx q[2];
rz(-1.1455673) q[2];
sx q[2];
rz(1.2644757) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1026417) q[1];
sx q[1];
rz(-1.1167553) q[1];
sx q[1];
rz(-0.42788831) q[1];
rz(3.0311534) q[3];
sx q[3];
rz(-0.38023284) q[3];
sx q[3];
rz(1.5766174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(0.28960323) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.7864236) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8913169) q[0];
sx q[0];
rz(-0.96069562) q[0];
sx q[0];
rz(0.079205714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(-3.0485857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(-2.0204087) q[1];
rz(-pi) q[2];
x q[2];
rz(0.077640688) q[3];
sx q[3];
rz(-1.75378) q[3];
sx q[3];
rz(1.3444572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168034) q[0];
sx q[0];
rz(-1.5644801) q[0];
sx q[0];
rz(1.9301231) q[0];
rz(0.89945729) q[2];
sx q[2];
rz(-0.81842917) q[2];
sx q[2];
rz(-2.1607272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6667337) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(0.95400793) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8912656) q[3];
sx q[3];
rz(-2.2776789) q[3];
sx q[3];
rz(1.3644497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-0.24969077) q[0];
rz(1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(-3.1304741) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7751559) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(-2.0648271) q[0];
rz(-pi) q[1];
rz(-2.3035994) q[2];
sx q[2];
rz(-1.0655155) q[2];
sx q[2];
rz(2.879564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74777714) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(1.0541037) q[1];
rz(-0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.8146851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622612) q[0];
sx q[0];
rz(-0.766303) q[0];
sx q[0];
rz(2.9192231) q[0];
x q[1];
rz(-2.1660216) q[2];
sx q[2];
rz(-0.86704463) q[2];
sx q[2];
rz(3.1095568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8618968) q[1];
sx q[1];
rz(-1.8021291) q[1];
sx q[1];
rz(0.33613236) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1633478) q[3];
sx q[3];
rz(-1.3031928) q[3];
sx q[3];
rz(1.0405376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(-1.2549531) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-0.83165283) q[0];
rz(-0.90066465) q[2];
sx q[2];
rz(-1.2256983) q[2];
sx q[2];
rz(-0.64054856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0259077) q[1];
sx q[1];
rz(-0.58870643) q[1];
sx q[1];
rz(3.1097417) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7498155) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(-2.9359407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(2.4672467) q[0];
rz(-2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(2.5792714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4903957) q[0];
sx q[0];
rz(-0.90264747) q[0];
sx q[0];
rz(3.1227123) q[0];
x q[1];
rz(-0.28568761) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(3.0722741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9971804) q[1];
sx q[1];
rz(-1.715853) q[1];
sx q[1];
rz(-0.1952862) q[1];
rz(-2.167422) q[3];
sx q[3];
rz(-1.6430292) q[3];
sx q[3];
rz(-0.40856397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.6202392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521261) q[0];
sx q[0];
rz(-2.9331101) q[0];
sx q[0];
rz(0.33574386) q[0];
x q[1];
rz(-2.7395758) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.5356262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7594787) q[1];
sx q[1];
rz(-1.6303051) q[1];
sx q[1];
rz(2.3630523) q[1];
x q[2];
rz(0.11407125) q[3];
sx q[3];
rz(-2.2469006) q[3];
sx q[3];
rz(-1.3249719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7662738) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(-0.71592285) q[0];
rz(-pi) q[1];
rz(-0.46220772) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(-0.30009899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7120413) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(0.045750381) q[1];
x q[2];
rz(-0.98468303) q[3];
sx q[3];
rz(-2.1911748) q[3];
sx q[3];
rz(-1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-0.45544004) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586944) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(-1.0844896) q[0];
rz(-pi) q[1];
rz(3.1206563) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(2.7431969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1315688) q[1];
sx q[1];
rz(-1.459517) q[1];
sx q[1];
rz(-0.93025031) q[1];
x q[2];
rz(-1.6845735) q[3];
sx q[3];
rz(-0.69996951) q[3];
sx q[3];
rz(-1.9663119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0192169) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-1.238263) q[2];
sx q[2];
rz(-2.5730972) q[2];
sx q[2];
rz(2.2039883) q[2];
rz(0.76392382) q[3];
sx q[3];
rz(-2.7803956) q[3];
sx q[3];
rz(-1.9261388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
