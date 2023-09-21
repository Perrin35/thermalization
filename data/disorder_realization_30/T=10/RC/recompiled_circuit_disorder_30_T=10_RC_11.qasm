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
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(-1.2000097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29205706) q[0];
sx q[0];
rz(-1.5229051) q[0];
sx q[0];
rz(-0.0056664771) q[0];
rz(-2.0179022) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(-2.6927039) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(-1.0784472) q[1];
x q[2];
rz(-1.6148189) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8346943) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(0.28960323) q[2];
rz(0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(-1.7864236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8913169) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-3.0623869) q[0];
rz(-pi) q[1];
rz(2.9789574) q[2];
sx q[2];
rz(-1.6911104) q[2];
sx q[2];
rz(1.4974809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.6347486) q[1];
sx q[1];
rz(-1.121184) q[1];
rz(-pi) q[2];
rz(-1.7543206) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(-2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(-2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(2.5007201) q[0];
rz(0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0707866) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(1.552836) q[0];
rz(-pi) q[1];
rz(-2.5550585) q[2];
sx q[2];
rz(-0.9622935) q[2];
sx q[2];
rz(1.8412794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9948431) q[1];
sx q[1];
rz(-2.5238876) q[1];
sx q[1];
rz(1.5084933) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8912656) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(1.3644497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5549094) q[0];
sx q[0];
rz(-1.1197829) q[0];
sx q[0];
rz(0.45278544) q[0];
x q[1];
rz(-2.501802) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(1.421979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1383458) q[1];
sx q[1];
rz(-2.0577288) q[1];
sx q[1];
rz(0.37133118) q[1];
x q[2];
rz(2.8233077) q[3];
sx q[3];
rz(-1.7537698) q[3];
sx q[3];
rz(1.1312248) q[3];
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
rz(-0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0622612) q[0];
sx q[0];
rz(-0.766303) q[0];
sx q[0];
rz(-2.9192231) q[0];
rz(-pi) q[1];
rz(0.97557108) q[2];
sx q[2];
rz(-0.86704463) q[2];
sx q[2];
rz(3.1095568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2697619) q[1];
sx q[1];
rz(-0.40553906) q[1];
sx q[1];
rz(-2.521442) q[1];
rz(-pi) q[2];
rz(0.96552403) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(1.0799288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.241248) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(-0.83165283) q[0];
x q[1];
rz(-0.90066465) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(0.64054856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9876154) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(1.5495367) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4721018) q[3];
sx q[3];
rz(-1.7117501) q[3];
sx q[3];
rz(1.6657176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092098504) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(2.2390319) q[0];
rz(-2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-3.0722741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6866236) q[1];
sx q[1];
rz(-1.7640055) q[1];
sx q[1];
rz(1.7186233) q[1];
rz(0.087248487) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(-2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8521261) q[0];
sx q[0];
rz(-0.20848256) q[0];
sx q[0];
rz(0.33574386) q[0];
rz(-pi) q[1];
rz(0.40201681) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.5356262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24890451) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(-3.0569539) q[1];
x q[2];
rz(-3.0275214) q[3];
sx q[3];
rz(-0.89469203) q[3];
sx q[3];
rz(1.3249719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(0.35933581) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(-2.8709581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0918002) q[0];
sx q[0];
rz(-1.5349749) q[0];
sx q[0];
rz(1.5396176) q[0];
rz(-pi) q[1];
rz(0.46220772) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(0.30009899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4295514) q[1];
sx q[1];
rz(-1.5263285) q[1];
sx q[1];
rz(3.0958423) q[1];
x q[2];
rz(-0.6587894) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(2.7487019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(-0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(0.49490067) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56787581) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(3.0987415) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5636744) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(-0.42051007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4329263) q[1];
sx q[1];
rz(-0.64879829) q[1];
sx q[1];
rz(-1.385958) q[1];
x q[2];
rz(3.0462618) q[3];
sx q[3];
rz(-0.87626002) q[3];
sx q[3];
rz(-1.3235843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-1.9033296) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(2.3776688) q[3];
sx q[3];
rz(-0.36119701) q[3];
sx q[3];
rz(1.2154538) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];