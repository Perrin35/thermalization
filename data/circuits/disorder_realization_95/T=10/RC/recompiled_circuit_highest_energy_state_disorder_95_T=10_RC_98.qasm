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
rz(0.79701841) q[0];
sx q[0];
rz(-2.2198644) q[0];
sx q[0];
rz(-1.8812688) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(0.33198196) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2352358) q[0];
sx q[0];
rz(-2.7180494) q[0];
sx q[0];
rz(-1.5111501) q[0];
rz(0.44029616) q[2];
sx q[2];
rz(-0.74661359) q[2];
sx q[2];
rz(2.8842215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1412072) q[1];
sx q[1];
rz(-1.5391333) q[1];
sx q[1];
rz(2.382676) q[1];
rz(-pi) q[2];
rz(1.6360248) q[3];
sx q[3];
rz(-2.0121775) q[3];
sx q[3];
rz(-0.38204398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6648286) q[2];
sx q[2];
rz(-2.2085184) q[2];
sx q[2];
rz(-2.942371) q[2];
rz(0.18579379) q[3];
sx q[3];
rz(-1.7265065) q[3];
sx q[3];
rz(1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.324447) q[0];
sx q[0];
rz(-2.6528093) q[0];
sx q[0];
rz(2.4031438) q[0];
rz(-1.6959408) q[1];
sx q[1];
rz(-1.7414469) q[1];
sx q[1];
rz(0.48283985) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8890742) q[0];
sx q[0];
rz(-1.7216428) q[0];
sx q[0];
rz(2.9993527) q[0];
rz(-pi) q[1];
rz(-1.1170399) q[2];
sx q[2];
rz(-0.38402176) q[2];
sx q[2];
rz(-0.45986097) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0625667) q[1];
sx q[1];
rz(-1.1084659) q[1];
sx q[1];
rz(-0.61840017) q[1];
x q[2];
rz(2.9779997) q[3];
sx q[3];
rz(-1.0181659) q[3];
sx q[3];
rz(1.9528509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26756755) q[2];
sx q[2];
rz(-1.678062) q[2];
sx q[2];
rz(-1.2919424) q[2];
rz(-2.0180295) q[3];
sx q[3];
rz(-2.4396887) q[3];
sx q[3];
rz(-1.7263713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8181151) q[0];
sx q[0];
rz(-0.98387843) q[0];
sx q[0];
rz(1.6997319) q[0];
rz(-1.2280751) q[1];
sx q[1];
rz(-1.2308729) q[1];
sx q[1];
rz(2.0679811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5283877) q[0];
sx q[0];
rz(-2.3484485) q[0];
sx q[0];
rz(-0.49164495) q[0];
rz(-pi) q[1];
rz(1.8262489) q[2];
sx q[2];
rz(-0.67136229) q[2];
sx q[2];
rz(0.76847149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6317128) q[1];
sx q[1];
rz(-0.86871618) q[1];
sx q[1];
rz(0.04336642) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32322901) q[3];
sx q[3];
rz(-2.6628837) q[3];
sx q[3];
rz(1.9883101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8492665) q[2];
sx q[2];
rz(-2.8385415) q[2];
sx q[2];
rz(-2.1185875) q[2];
rz(-0.060221378) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(0.72505081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66134727) q[0];
sx q[0];
rz(-2.9954973) q[0];
sx q[0];
rz(-0.51937854) q[0];
rz(0.6005148) q[1];
sx q[1];
rz(-2.0527716) q[1];
sx q[1];
rz(2.3868938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31904991) q[0];
sx q[0];
rz(-0.62915914) q[0];
sx q[0];
rz(0.99790093) q[0];
rz(-pi) q[1];
rz(1.7394051) q[2];
sx q[2];
rz(-2.4314432) q[2];
sx q[2];
rz(0.53393364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2024196) q[1];
sx q[1];
rz(-1.8853097) q[1];
sx q[1];
rz(1.0553204) q[1];
x q[2];
rz(-1.9672606) q[3];
sx q[3];
rz(-2.0959586) q[3];
sx q[3];
rz(-1.7523867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34919136) q[2];
sx q[2];
rz(-2.3442522) q[2];
sx q[2];
rz(0.1758197) q[2];
rz(-2.0154121) q[3];
sx q[3];
rz(-1.7837985) q[3];
sx q[3];
rz(0.064182909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675451) q[0];
sx q[0];
rz(-1.6459246) q[0];
sx q[0];
rz(-2.5597036) q[0];
rz(1.9582845) q[1];
sx q[1];
rz(-0.71154037) q[1];
sx q[1];
rz(1.8875095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81800705) q[0];
sx q[0];
rz(-0.76420751) q[0];
sx q[0];
rz(-1.0952522) q[0];
rz(2.6627867) q[2];
sx q[2];
rz(-1.1421912) q[2];
sx q[2];
rz(-1.3474825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0662567) q[1];
sx q[1];
rz(-1.7793223) q[1];
sx q[1];
rz(-1.8715402) q[1];
rz(-1.759911) q[3];
sx q[3];
rz(-2.0678036) q[3];
sx q[3];
rz(-0.86732098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1648272) q[2];
sx q[2];
rz(-2.3259614) q[2];
sx q[2];
rz(-2.7867479) q[2];
rz(-1.1497078) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(-0.84158516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717764) q[0];
sx q[0];
rz(-0.27158296) q[0];
sx q[0];
rz(-0.040104453) q[0];
rz(-1.4647723) q[1];
sx q[1];
rz(-2.1720839) q[1];
sx q[1];
rz(0.47223314) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036333648) q[0];
sx q[0];
rz(-1.8324295) q[0];
sx q[0];
rz(-0.18966578) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5750999) q[2];
sx q[2];
rz(-2.488778) q[2];
sx q[2];
rz(3.0497361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4716954) q[1];
sx q[1];
rz(-1.9006839) q[1];
sx q[1];
rz(-2.8389588) q[1];
rz(-pi) q[2];
rz(-1.8566441) q[3];
sx q[3];
rz(-0.82995358) q[3];
sx q[3];
rz(0.85808104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8670696) q[2];
sx q[2];
rz(-2.3752866) q[2];
sx q[2];
rz(1.6609894) q[2];
rz(0.04145043) q[3];
sx q[3];
rz(-0.71235123) q[3];
sx q[3];
rz(-2.1672772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.239045) q[0];
sx q[0];
rz(-0.79896611) q[0];
sx q[0];
rz(0.7582742) q[0];
rz(-0.43765086) q[1];
sx q[1];
rz(-1.9270908) q[1];
sx q[1];
rz(0.95380107) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4364638) q[0];
sx q[0];
rz(-2.7565694) q[0];
sx q[0];
rz(1.692311) q[0];
x q[1];
rz(-0.96289159) q[2];
sx q[2];
rz(-2.1959119) q[2];
sx q[2];
rz(1.6776379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0702677) q[1];
sx q[1];
rz(-2.4355781) q[1];
sx q[1];
rz(1.6672177) q[1];
x q[2];
rz(0.28080583) q[3];
sx q[3];
rz(-1.1698517) q[3];
sx q[3];
rz(-1.795648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7268251) q[2];
sx q[2];
rz(-2.813952) q[2];
sx q[2];
rz(-0.3978351) q[2];
rz(1.2658524) q[3];
sx q[3];
rz(-1.419302) q[3];
sx q[3];
rz(2.6771767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9809113) q[0];
sx q[0];
rz(-2.2137764) q[0];
sx q[0];
rz(2.8371147) q[0];
rz(-2.0092087) q[1];
sx q[1];
rz(-0.67190036) q[1];
sx q[1];
rz(2.0679881) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17504263) q[0];
sx q[0];
rz(-1.5538408) q[0];
sx q[0];
rz(-1.5938363) q[0];
x q[1];
rz(1.2108285) q[2];
sx q[2];
rz(-1.9641293) q[2];
sx q[2];
rz(-0.29811146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5636812) q[1];
sx q[1];
rz(-1.8133687) q[1];
sx q[1];
rz(0.96176186) q[1];
x q[2];
rz(1.520548) q[3];
sx q[3];
rz(-2.2138688) q[3];
sx q[3];
rz(-2.1724043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4094746) q[2];
sx q[2];
rz(-1.858859) q[2];
sx q[2];
rz(-0.049962433) q[2];
rz(0.91160715) q[3];
sx q[3];
rz(-0.92602366) q[3];
sx q[3];
rz(2.2302506) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6066345) q[0];
sx q[0];
rz(-1.7462523) q[0];
sx q[0];
rz(0.362679) q[0];
rz(-0.72049385) q[1];
sx q[1];
rz(-0.81169218) q[1];
sx q[1];
rz(1.1788751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900951) q[0];
sx q[0];
rz(-1.0637366) q[0];
sx q[0];
rz(0.66452615) q[0];
rz(0.32068738) q[2];
sx q[2];
rz(-1.4960519) q[2];
sx q[2];
rz(0.2439258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2612103) q[1];
sx q[1];
rz(-0.19621135) q[1];
sx q[1];
rz(-2.0878778) q[1];
rz(-pi) q[2];
rz(-1.9835155) q[3];
sx q[3];
rz(-1.0636119) q[3];
sx q[3];
rz(0.87178236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1371586) q[2];
sx q[2];
rz(-0.68778554) q[2];
sx q[2];
rz(1.9384109) q[2];
rz(0.016599003) q[3];
sx q[3];
rz(-1.6417475) q[3];
sx q[3];
rz(-1.2703007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59852973) q[0];
sx q[0];
rz(-0.5235343) q[0];
sx q[0];
rz(-1.6981) q[0];
rz(2.6632925) q[1];
sx q[1];
rz(-0.6820448) q[1];
sx q[1];
rz(-1.2102478) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87916527) q[0];
sx q[0];
rz(-2.0602003) q[0];
sx q[0];
rz(0.64206815) q[0];
x q[1];
rz(-1.7487594) q[2];
sx q[2];
rz(-2.2134668) q[2];
sx q[2];
rz(2.8699574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4251488) q[1];
sx q[1];
rz(-1.8952491) q[1];
sx q[1];
rz(-1.6258045) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7948512) q[3];
sx q[3];
rz(-1.2478078) q[3];
sx q[3];
rz(0.59524465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18834194) q[2];
sx q[2];
rz(-0.61389273) q[2];
sx q[2];
rz(-2.4284412) q[2];
rz(0.6330511) q[3];
sx q[3];
rz(-0.96708599) q[3];
sx q[3];
rz(-0.87575325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0093832) q[0];
sx q[0];
rz(-1.3306946) q[0];
sx q[0];
rz(2.7015986) q[0];
rz(-0.72883365) q[1];
sx q[1];
rz(-2.3757917) q[1];
sx q[1];
rz(1.052312) q[1];
rz(1.960878) q[2];
sx q[2];
rz(-1.6772391) q[2];
sx q[2];
rz(-0.46628484) q[2];
rz(-0.36076693) q[3];
sx q[3];
rz(-0.84070342) q[3];
sx q[3];
rz(-2.5221205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
