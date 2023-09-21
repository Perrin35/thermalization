OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(2.7471623) q[0];
rz(-3.1354304) q[1];
sx q[1];
rz(-2.8013464) q[1];
sx q[1];
rz(-1.9415829) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7317176) q[0];
sx q[0];
rz(-3.0933676) q[0];
sx q[0];
rz(1.6884786) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0179022) q[2];
sx q[2];
rz(-1.2667709) q[2];
sx q[2];
rz(0.44888874) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1026417) q[1];
sx q[1];
rz(-2.0248374) q[1];
sx q[1];
rz(2.7137043) q[1];
rz(-1.5267738) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.7864236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8913169) q[0];
sx q[0];
rz(-0.96069562) q[0];
sx q[0];
rz(-0.079205714) q[0];
rz(-0.16263527) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.4974809) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.558074) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(-3.0706057) q[1];
rz(-pi) q[2];
rz(-1.9676898) q[3];
sx q[3];
rz(-0.19860425) q[3];
sx q[3];
rz(-1.7484776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(0.64087254) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-2.0764988) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6247892) q[0];
sx q[0];
rz(-1.5771126) q[0];
sx q[0];
rz(-1.9301231) q[0];
rz(0.58653411) q[2];
sx q[2];
rz(-2.1792992) q[2];
sx q[2];
rz(-1.8412794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6667337) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(2.1875847) q[1];
x q[2];
rz(-0.84824003) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(-0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(-3.1304741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5866833) q[0];
sx q[0];
rz(-2.0218098) q[0];
sx q[0];
rz(-2.6888072) q[0];
rz(2.501802) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(-1.7196136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1383458) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(-2.7702615) q[1];
rz(-2.8233077) q[3];
sx q[3];
rz(-1.7537698) q[3];
sx q[3];
rz(2.0103679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3834284) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(1.3616256) q[0];
rz(0.58381501) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-2.3655287) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2697619) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(0.62015066) q[1];
x q[2];
rz(-2.8513961) q[3];
sx q[3];
rz(-1.1786596) q[3];
sx q[3];
rz(-0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9003446) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(2.3099398) q[0];
rz(-pi) q[1];
rz(0.90066465) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(-0.64054856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7129732) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-2.5531205) q[1];
x q[2];
rz(-1.3917771) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091244) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599265) q[0];
sx q[0];
rz(-2.473218) q[0];
sx q[0];
rz(-1.5468803) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4698896) q[2];
sx q[2];
rz(-1.2404053) q[2];
sx q[2];
rz(0.23320564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14441227) q[1];
sx q[1];
rz(-1.4257396) q[1];
sx q[1];
rz(-2.9463065) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.167422) q[3];
sx q[3];
rz(-1.4985634) q[3];
sx q[3];
rz(0.40856397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.6202392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521261) q[0];
sx q[0];
rz(-2.9331101) q[0];
sx q[0];
rz(0.33574386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40201681) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.6059665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.382114) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(-0.77854034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(-1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3753189) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(-2.4256698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46220772) q[2];
sx q[2];
rz(-2.2773829) q[2];
sx q[2];
rz(-2.8414937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4295514) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(-3.0958423) q[1];
rz(-0.98468303) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(-1.8970722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(0.49490067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5737168) q[0];
sx q[0];
rz(-1.084869) q[0];
sx q[0];
rz(0.042851187) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5779183) q[2];
sx q[2];
rz(-1.2429503) q[2];
sx q[2];
rz(2.7210826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1315688) q[1];
sx q[1];
rz(-1.6820757) q[1];
sx q[1];
rz(0.93025031) q[1];
rz(-pi) q[2];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(1.3235843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-1.238263) q[2];
sx q[2];
rz(-2.5730972) q[2];
sx q[2];
rz(2.2039883) q[2];
rz(-0.76392382) q[3];
sx q[3];
rz(-0.36119701) q[3];
sx q[3];
rz(1.2154538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
