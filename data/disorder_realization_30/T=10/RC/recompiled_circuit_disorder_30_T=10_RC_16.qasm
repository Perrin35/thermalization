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
rz(2.8770652) q[0];
sx q[0];
rz(9.8192083) q[0];
rz(-3.1354304) q[1];
sx q[1];
rz(-2.8013464) q[1];
sx q[1];
rz(-1.9415829) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7317176) q[0];
sx q[0];
rz(-3.0933676) q[0];
sx q[0];
rz(1.4531141) q[0];
rz(2.0179022) q[2];
sx q[2];
rz(-1.2667709) q[2];
sx q[2];
rz(-2.6927039) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0389509) q[1];
sx q[1];
rz(-2.0248374) q[1];
sx q[1];
rz(-2.7137043) q[1];
rz(-pi) q[2];
rz(1.6148189) q[3];
sx q[3];
rz(-1.1929973) q[3];
sx q[3];
rz(-1.4461185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(2.7976024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(-1.7864236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2750759) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(-0.95922031) q[0];
x q[1];
rz(1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.720571) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(-1.717091) q[1];
rz(-pi) q[2];
rz(1.7543206) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.552836) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87419072) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(-3.0490321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0703735) q[1];
sx q[1];
rz(-2.1871236) q[1];
sx q[1];
rz(-3.0973869) q[1];
rz(-pi) q[2];
rz(1.2885116) q[3];
sx q[3];
rz(-0.74263393) q[3];
sx q[3];
rz(1.4020855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(-2.4242145) q[2];
rz(-2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(3.1304741) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7751559) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(-2.0648271) q[0];
rz(-0.83799329) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(2.879564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3938155) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(-2.087489) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6075881) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(-0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(0.66343534) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
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
sx q[2];
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
rz(-1.8146851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3834284) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(1.3616256) q[0];
rz(2.3438498) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(1.1898578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2111645) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(1.1633478) q[3];
sx q[3];
rz(-1.8383998) q[3];
sx q[3];
rz(2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(1.1157657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(-2.1372165) q[0];
rz(-pi) q[1];
rz(2.240928) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(0.64054856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0259077) q[1];
sx q[1];
rz(-2.5528862) q[1];
sx q[1];
rz(-3.1097417) q[1];
rz(-pi) q[2];
rz(-2.4721018) q[3];
sx q[3];
rz(-1.7117501) q[3];
sx q[3];
rz(-1.475875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(-0.3195233) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(2.4672467) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(0.56232125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6816662) q[0];
sx q[0];
rz(-0.66837464) q[0];
sx q[0];
rz(1.5468803) q[0];
rz(-1.6717031) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(-2.908387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9971804) q[1];
sx q[1];
rz(-1.4257396) q[1];
sx q[1];
rz(-2.9463065) q[1];
rz(-pi) q[2];
rz(0.9741707) q[3];
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
rz(2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.5213535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28946653) q[0];
sx q[0];
rz(-2.9331101) q[0];
sx q[0];
rz(0.33574386) q[0];
rz(2.2659726) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(-2.3919174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.382114) q[1];
sx q[1];
rz(-1.6303051) q[1];
sx q[1];
rz(-0.77854034) q[1];
rz(-pi) q[2];
rz(-1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-0.27063453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04979241) q[0];
sx q[0];
rz(-1.6066178) q[0];
sx q[0];
rz(-1.5396176) q[0];
x q[1];
rz(2.052202) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(-0.35441986) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86079019) q[1];
sx q[1];
rz(-1.6165015) q[1];
sx q[1];
rz(-1.6153107) q[1];
rz(-2.4828033) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6311326) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4821645) q[0];
sx q[0];
rz(-0.48766252) q[0];
sx q[0];
rz(1.4898666) q[0];
rz(3.1206563) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(2.7431969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1315688) q[1];
sx q[1];
rz(-1.459517) q[1];
sx q[1];
rz(-2.2113423) q[1];
rz(-pi) q[2];
rz(-0.095330843) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(1.3235843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(2.949529) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
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
rz(1.8264063) q[3];
sx q[3];
rz(-1.3127463) q[3];
sx q[3];
rz(-1.128872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
