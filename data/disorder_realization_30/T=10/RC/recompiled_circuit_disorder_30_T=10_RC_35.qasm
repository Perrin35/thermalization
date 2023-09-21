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
rz(-1.2000097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40987504) q[0];
sx q[0];
rz(-0.048225064) q[0];
sx q[0];
rz(1.4531141) q[0];
x q[1];
rz(-0.33485246) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(-1.877117) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4760121) q[1];
sx q[1];
rz(-1.1886547) q[1];
sx q[1];
rz(-2.0631454) q[1];
rz(-pi) q[2];
rz(-0.37813152) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(0.10842987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.6480185) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(-1.3551691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38793135) q[0];
sx q[0];
rz(-2.5270215) q[0];
sx q[0];
rz(1.6835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(3.0485857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1234839) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(-2.0204087) q[1];
x q[2];
rz(1.387272) q[3];
sx q[3];
rz(-1.4944544) q[3];
sx q[3];
rz(-2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(0.74869853) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(2.0764988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0707866) q[0];
sx q[0];
rz(-0.35937989) q[0];
sx q[0];
rz(1.552836) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2674019) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(-3.0490321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0703735) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
rz(-0.84824003) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(-0.041785985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(2.1266134) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(3.1304741) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7751559) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(2.0648271) q[0];
rz(-pi) q[1];
rz(2.3035994) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(2.879564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3938155) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(2.087489) q[1];
rz(1.7632145) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(-0.49945143) q[3];
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
rz(2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304612) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.3269075) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793314) q[0];
sx q[0];
rz(-0.766303) q[0];
sx q[0];
rz(0.22236951) q[0];
rz(-pi) q[1];
rz(0.7977428) q[2];
sx q[2];
rz(-1.1290871) q[2];
sx q[2];
rz(1.9517348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9304282) q[1];
sx q[1];
rz(-1.2439562) q[1];
sx q[1];
rz(-1.3262761) q[1];
rz(-pi) q[2];
rz(2.8513961) q[3];
sx q[3];
rz(-1.1786596) q[3];
sx q[3];
rz(0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(1.1157657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(2.1372165) q[0];
rz(-2.240928) q[2];
sx q[2];
rz(-1.9158944) q[2];
sx q[2];
rz(-0.64054856) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-0.58847217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22478215) q[3];
sx q[3];
rz(-2.4596679) q[3];
sx q[3];
rz(3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091244) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(0.56232125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651197) q[0];
sx q[0];
rz(-0.90264747) q[0];
sx q[0];
rz(-3.1227123) q[0];
rz(-1.4698896) q[2];
sx q[2];
rz(-1.9011874) q[2];
sx q[2];
rz(-0.23320564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79531419) q[1];
sx q[1];
rz(-2.8988796) q[1];
sx q[1];
rz(-2.4962891) q[1];
rz(0.087248487) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-2.8015461) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(-1.5213535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(1.6403857) q[0];
rz(-2.0422158) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(-2.9277756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7594787) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(-2.3630523) q[1];
rz(-pi) q[2];
rz(1.429854) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(1.9977026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(-2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(0.27063453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3753189) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(0.71592285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6793849) q[2];
sx q[2];
rz(-0.86420977) q[2];
sx q[2];
rz(-2.8414937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.088085727) q[1];
sx q[1];
rz(-3.0778031) q[1];
sx q[1];
rz(2.3699058) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1569096) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(2.646692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65942818) q[0];
sx q[0];
rz(-2.6539301) q[0];
sx q[0];
rz(1.4898666) q[0];
rz(-pi) q[1];
rz(1.5636744) q[2];
sx q[2];
rz(-1.8986423) q[2];
sx q[2];
rz(0.42051007) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.478185) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(0.138476) q[1];
rz(0.095330843) q[3];
sx q[3];
rz(-0.87626002) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-2.1140425) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(-2.8752747) q[3];
sx q[3];
rz(-1.8177633) q[3];
sx q[3];
rz(0.37533356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
