OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.42263) q[0];
sx q[0];
rz(-2.842272) q[0];
sx q[0];
rz(2.646995) q[0];
rz(-1.9994796) q[1];
sx q[1];
rz(-2.1358868) q[1];
sx q[1];
rz(-1.1297273) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725192) q[0];
sx q[0];
rz(-1.7290218) q[0];
sx q[0];
rz(-2.4618966) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3006265) q[2];
sx q[2];
rz(-0.6305002) q[2];
sx q[2];
rz(-1.7918685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3898824) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(-0.33716476) q[1];
x q[2];
rz(2.1100329) q[3];
sx q[3];
rz(-2.4162216) q[3];
sx q[3];
rz(1.4998719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31892458) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(2.2887716) q[2];
rz(-1.3301814) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(-0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80832076) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.6774696) q[0];
rz(1.4785712) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(-1.9333855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.941146) q[0];
sx q[0];
rz(-1.0387372) q[0];
sx q[0];
rz(-1.0131628) q[0];
rz(-pi) q[1];
rz(-1.1651498) q[2];
sx q[2];
rz(-2.6313836) q[2];
sx q[2];
rz(-1.134553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8025) q[1];
sx q[1];
rz(-2.8475757) q[1];
sx q[1];
rz(1.1400998) q[1];
rz(-pi) q[2];
rz(0.45470806) q[3];
sx q[3];
rz(-1.4147926) q[3];
sx q[3];
rz(1.7059513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53604424) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.14979664) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-0.27161828) q[0];
rz(-2.244921) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(-1.8439878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51894278) q[0];
sx q[0];
rz(-1.1077767) q[0];
sx q[0];
rz(0.011269022) q[0];
x q[1];
rz(2.5141352) q[2];
sx q[2];
rz(-0.59840032) q[2];
sx q[2];
rz(-2.5638169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5601968) q[1];
sx q[1];
rz(-0.27276688) q[1];
sx q[1];
rz(1.1352989) q[1];
rz(-pi) q[2];
rz(0.0053169189) q[3];
sx q[3];
rz(-2.379619) q[3];
sx q[3];
rz(-0.055082037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74035949) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(-1.8476123) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(-0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2006705) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(-1.6773552) q[0];
rz(-1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(-0.11437036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3791311) q[0];
sx q[0];
rz(-1.987559) q[0];
sx q[0];
rz(-2.5083191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.52621) q[2];
sx q[2];
rz(-0.75591959) q[2];
sx q[2];
rz(2.4494954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83437551) q[1];
sx q[1];
rz(-2.465472) q[1];
sx q[1];
rz(1.6130877) q[1];
rz(1.5642197) q[3];
sx q[3];
rz(-1.7688171) q[3];
sx q[3];
rz(-1.1634367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6681246) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(-0.075695666) q[2];
rz(-0.43478742) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7441523) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(2.721526) q[0];
rz(2.9282667) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(-1.2423645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38753375) q[0];
sx q[0];
rz(-2.2926169) q[0];
sx q[0];
rz(-0.57768627) q[0];
rz(-pi) q[1];
rz(1.2790643) q[2];
sx q[2];
rz(-1.1559249) q[2];
sx q[2];
rz(1.3049098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27659163) q[1];
sx q[1];
rz(-1.5861771) q[1];
sx q[1];
rz(-1.6713784) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20988864) q[3];
sx q[3];
rz(-0.22612962) q[3];
sx q[3];
rz(1.1864288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3199978) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(1.0125259) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(2.4969126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(2.7110355) q[0];
rz(-0.14398362) q[1];
sx q[1];
rz(-2.0591683) q[1];
sx q[1];
rz(0.63124257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743606) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(1.505018) q[0];
x q[1];
rz(-1.2530042) q[2];
sx q[2];
rz(-2.7016407) q[2];
sx q[2];
rz(-0.41340128) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65578016) q[1];
sx q[1];
rz(-2.1596599) q[1];
sx q[1];
rz(-2.1766125) q[1];
rz(-pi) q[2];
rz(-2.0231163) q[3];
sx q[3];
rz(-1.8743487) q[3];
sx q[3];
rz(-1.1293052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97079078) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(1.4405174) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-1.1037339) q[3];
sx q[3];
rz(0.48719278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454527) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(0.92426306) q[0];
rz(-0.29280064) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(2.3407095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594052) q[0];
sx q[0];
rz(-2.6818083) q[0];
sx q[0];
rz(-0.96084382) q[0];
rz(-pi) q[1];
rz(2.0037193) q[2];
sx q[2];
rz(-0.70919631) q[2];
sx q[2];
rz(1.8877754) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1276857) q[1];
sx q[1];
rz(-0.10914762) q[1];
sx q[1];
rz(2.8692895) q[1];
rz(-pi) q[2];
rz(3.1195779) q[3];
sx q[3];
rz(-1.5436633) q[3];
sx q[3];
rz(-0.61329816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31876365) q[2];
sx q[2];
rz(-1.9209361) q[2];
sx q[2];
rz(-1.8611543) q[2];
rz(-0.66796962) q[3];
sx q[3];
rz(-2.6536055) q[3];
sx q[3];
rz(-0.41332301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.2385383) q[0];
sx q[0];
rz(1.1827693) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(-0.059344083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6366258) q[0];
sx q[0];
rz(-1.4511746) q[0];
sx q[0];
rz(-1.8455532) q[0];
x q[1];
rz(-2.5586023) q[2];
sx q[2];
rz(-2.3263479) q[2];
sx q[2];
rz(1.5423403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4044734) q[1];
sx q[1];
rz(-1.307202) q[1];
sx q[1];
rz(-0.10837676) q[1];
rz(-pi) q[2];
rz(0.077499302) q[3];
sx q[3];
rz(-1.8935673) q[3];
sx q[3];
rz(-2.4609044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5817029) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(1.3524559) q[2];
rz(1.3672359) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(0.8555612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2818114) q[0];
sx q[0];
rz(-2.3563522) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(-3.1135318) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(-1.185816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738041) q[0];
sx q[0];
rz(-1.8361183) q[0];
sx q[0];
rz(1.6173043) q[0];
rz(-1.469428) q[2];
sx q[2];
rz(-1.6497872) q[2];
sx q[2];
rz(-1.7192507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8227905) q[1];
sx q[1];
rz(-1.0981961) q[1];
sx q[1];
rz(2.5986555) q[1];
x q[2];
rz(0.61798685) q[3];
sx q[3];
rz(-1.3387965) q[3];
sx q[3];
rz(-2.8232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(2.5458941) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(-0.025207635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6237685) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(0.18381707) q[0];
rz(0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-0.26430166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562132) q[0];
sx q[0];
rz(-2.1341679) q[0];
sx q[0];
rz(-1.8895288) q[0];
rz(-2.0980706) q[2];
sx q[2];
rz(-0.95607483) q[2];
sx q[2];
rz(1.5103024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0288913) q[1];
sx q[1];
rz(-1.750578) q[1];
sx q[1];
rz(1.4110908) q[1];
rz(-pi) q[2];
rz(-1.4749583) q[3];
sx q[3];
rz(-1.0776099) q[3];
sx q[3];
rz(3.1214023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(-1.4043572) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(-3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514773) q[0];
sx q[0];
rz(-2.293312) q[0];
sx q[0];
rz(1.5207416) q[0];
rz(-0.71113853) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(-2.0696832) q[2];
sx q[2];
rz(-1.6098534) q[2];
sx q[2];
rz(0.44865566) q[2];
rz(1.0300954) q[3];
sx q[3];
rz(-2.0602915) q[3];
sx q[3];
rz(0.66707053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
