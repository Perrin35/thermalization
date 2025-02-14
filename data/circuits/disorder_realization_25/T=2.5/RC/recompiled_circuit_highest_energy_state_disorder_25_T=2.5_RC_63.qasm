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
rz(2.0410886) q[0];
sx q[0];
rz(-0.70900822) q[0];
sx q[0];
rz(-2.3472743) q[0];
rz(0.89214832) q[1];
sx q[1];
rz(-1.661707) q[1];
sx q[1];
rz(2.7977112) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1411439) q[0];
sx q[0];
rz(-1.6578995) q[0];
sx q[0];
rz(-2.8892975) q[0];
rz(2.6394541) q[2];
sx q[2];
rz(-1.0912899) q[2];
sx q[2];
rz(1.5554847) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81864322) q[1];
sx q[1];
rz(-1.4122077) q[1];
sx q[1];
rz(-1.4399546) q[1];
x q[2];
rz(0.12012336) q[3];
sx q[3];
rz(-1.7645482) q[3];
sx q[3];
rz(1.4258949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5188562) q[2];
sx q[2];
rz(-1.3750261) q[2];
sx q[2];
rz(-2.1905017) q[2];
rz(-2.2226492) q[3];
sx q[3];
rz(-0.94075847) q[3];
sx q[3];
rz(-2.9520891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6491991) q[0];
sx q[0];
rz(-2.352582) q[0];
sx q[0];
rz(-1.5845818) q[0];
rz(-0.7400662) q[1];
sx q[1];
rz(-2.2941755) q[1];
sx q[1];
rz(-1.7795631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131378) q[0];
sx q[0];
rz(-1.0759996) q[0];
sx q[0];
rz(2.6005756) q[0];
rz(1.9744423) q[2];
sx q[2];
rz(-1.4440618) q[2];
sx q[2];
rz(2.6069146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6891514) q[1];
sx q[1];
rz(-1.930962) q[1];
sx q[1];
rz(-2.4634201) q[1];
rz(-pi) q[2];
x q[2];
rz(0.858694) q[3];
sx q[3];
rz(-1.5017548) q[3];
sx q[3];
rz(0.25180975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43807277) q[2];
sx q[2];
rz(-1.1601996) q[2];
sx q[2];
rz(0.6456708) q[2];
rz(0.034363184) q[3];
sx q[3];
rz(-2.2785701) q[3];
sx q[3];
rz(-0.063095108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768169) q[0];
sx q[0];
rz(-2.9314628) q[0];
sx q[0];
rz(-2.3576376) q[0];
rz(-0.3166554) q[1];
sx q[1];
rz(-1.4953934) q[1];
sx q[1];
rz(2.1276316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74815566) q[0];
sx q[0];
rz(-1.5542826) q[0];
sx q[0];
rz(-0.17513398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9892817) q[2];
sx q[2];
rz(-0.51562247) q[2];
sx q[2];
rz(-1.2116694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55401504) q[1];
sx q[1];
rz(-1.8178175) q[1];
sx q[1];
rz(0.26320158) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6211735) q[3];
sx q[3];
rz(-1.0380967) q[3];
sx q[3];
rz(2.1100489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.97518146) q[2];
sx q[2];
rz(-1.760249) q[2];
sx q[2];
rz(-0.4057917) q[2];
rz(2.8864012) q[3];
sx q[3];
rz(-1.2602592) q[3];
sx q[3];
rz(2.3658128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.880421) q[0];
sx q[0];
rz(-0.71735993) q[0];
sx q[0];
rz(-0.98947155) q[0];
rz(-2.5148892) q[1];
sx q[1];
rz(-0.53755212) q[1];
sx q[1];
rz(-0.42868838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1470448) q[0];
sx q[0];
rz(-0.8126077) q[0];
sx q[0];
rz(-1.7269686) q[0];
rz(-1.7103042) q[2];
sx q[2];
rz(-1.9006471) q[2];
sx q[2];
rz(-2.863186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3544412) q[1];
sx q[1];
rz(-2.0729182) q[1];
sx q[1];
rz(0.67014613) q[1];
x q[2];
rz(2.8270014) q[3];
sx q[3];
rz(-0.73986161) q[3];
sx q[3];
rz(3.0219697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2428212) q[2];
sx q[2];
rz(-0.12671825) q[2];
sx q[2];
rz(-0.80647331) q[2];
rz(1.1558695) q[3];
sx q[3];
rz(-2.3838796) q[3];
sx q[3];
rz(0.19494593) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72215801) q[0];
sx q[0];
rz(-2.2787155) q[0];
sx q[0];
rz(-0.3062329) q[0];
rz(-0.58079863) q[1];
sx q[1];
rz(-0.88169801) q[1];
sx q[1];
rz(0.64540783) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54242486) q[0];
sx q[0];
rz(-1.4087447) q[0];
sx q[0];
rz(-0.57181825) q[0];
rz(-0.077351582) q[2];
sx q[2];
rz(-1.0707375) q[2];
sx q[2];
rz(1.3429221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4731406) q[1];
sx q[1];
rz(-2.4048971) q[1];
sx q[1];
rz(-0.84676844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.018086334) q[3];
sx q[3];
rz(-1.340217) q[3];
sx q[3];
rz(0.54508506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8530276) q[2];
sx q[2];
rz(-1.7662851) q[2];
sx q[2];
rz(2.1993401) q[2];
rz(-0.45091584) q[3];
sx q[3];
rz(-0.98603606) q[3];
sx q[3];
rz(-0.26421419) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0838202) q[0];
sx q[0];
rz(-1.271893) q[0];
sx q[0];
rz(0.94495946) q[0];
rz(0.68929976) q[1];
sx q[1];
rz(-1.1079451) q[1];
sx q[1];
rz(-1.1030997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5430008) q[0];
sx q[0];
rz(-0.75816407) q[0];
sx q[0];
rz(-2.7757194) q[0];
rz(-1.1464153) q[2];
sx q[2];
rz(-1.7327899) q[2];
sx q[2];
rz(-1.4206518) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1018033) q[1];
sx q[1];
rz(-1.1032133) q[1];
sx q[1];
rz(2.1469248) q[1];
rz(-0.96007993) q[3];
sx q[3];
rz(-1.999049) q[3];
sx q[3];
rz(0.48510374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5276706) q[2];
sx q[2];
rz(-3.0752144) q[2];
sx q[2];
rz(1.4572246) q[2];
rz(-0.97113222) q[3];
sx q[3];
rz(-1.6077653) q[3];
sx q[3];
rz(-0.11626135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74294746) q[0];
sx q[0];
rz(-1.255641) q[0];
sx q[0];
rz(2.9435112) q[0];
rz(-2.6598341) q[1];
sx q[1];
rz(-1.878673) q[1];
sx q[1];
rz(-1.6162965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71512521) q[0];
sx q[0];
rz(-1.5670527) q[0];
sx q[0];
rz(3.1327898) q[0];
rz(-pi) q[1];
rz(0.003633066) q[2];
sx q[2];
rz(-1.1143408) q[2];
sx q[2];
rz(2.9287287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35514606) q[1];
sx q[1];
rz(-1.0736559) q[1];
sx q[1];
rz(-1.0721779) q[1];
rz(-pi) q[2];
rz(1.7500739) q[3];
sx q[3];
rz(-2.6411169) q[3];
sx q[3];
rz(-1.2795606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6306448) q[2];
sx q[2];
rz(-1.3628565) q[2];
sx q[2];
rz(-0.93713588) q[2];
rz(2.7142094) q[3];
sx q[3];
rz(-1.00777) q[3];
sx q[3];
rz(1.5905323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.19510022) q[0];
sx q[0];
rz(-1.6747403) q[0];
sx q[0];
rz(1.7505769) q[0];
rz(1.4637949) q[1];
sx q[1];
rz(-2.548806) q[1];
sx q[1];
rz(-2.6189199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928888) q[0];
sx q[0];
rz(-0.33451053) q[0];
sx q[0];
rz(-1.7728895) q[0];
rz(-pi) q[1];
rz(2.9982469) q[2];
sx q[2];
rz(-2.0269577) q[2];
sx q[2];
rz(0.83269115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9620888) q[1];
sx q[1];
rz(-0.25280646) q[1];
sx q[1];
rz(2.4180555) q[1];
rz(2.8925965) q[3];
sx q[3];
rz(-2.0820302) q[3];
sx q[3];
rz(1.5863904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78397822) q[2];
sx q[2];
rz(-1.7366333) q[2];
sx q[2];
rz(3.0493375) q[2];
rz(0.15010321) q[3];
sx q[3];
rz(-1.9893407) q[3];
sx q[3];
rz(0.86103719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1152231) q[0];
sx q[0];
rz(-2.8404591) q[0];
sx q[0];
rz(1.2979771) q[0];
rz(0.84833974) q[1];
sx q[1];
rz(-1.6611165) q[1];
sx q[1];
rz(-2.8678105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0299579) q[0];
sx q[0];
rz(-2.9024501) q[0];
sx q[0];
rz(-0.14683036) q[0];
rz(0.3329113) q[2];
sx q[2];
rz(-0.71117461) q[2];
sx q[2];
rz(-0.43363562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1787415) q[1];
sx q[1];
rz(-0.76707572) q[1];
sx q[1];
rz(0.13654582) q[1];
x q[2];
rz(-1.7785092) q[3];
sx q[3];
rz(-1.150032) q[3];
sx q[3];
rz(1.3265691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0724642) q[2];
sx q[2];
rz(-1.0651361) q[2];
sx q[2];
rz(-0.96159846) q[2];
rz(1.0570863) q[3];
sx q[3];
rz(-1.0870533) q[3];
sx q[3];
rz(-0.45007733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3707651) q[0];
sx q[0];
rz(-1.9285044) q[0];
sx q[0];
rz(-0.83592498) q[0];
rz(-1.4113374) q[1];
sx q[1];
rz(-0.58914271) q[1];
sx q[1];
rz(0.020847598) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6792471) q[0];
sx q[0];
rz(-1.6007375) q[0];
sx q[0];
rz(-1.7172555) q[0];
rz(2.8785275) q[2];
sx q[2];
rz(-2.6537958) q[2];
sx q[2];
rz(-2.8939025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67733586) q[1];
sx q[1];
rz(-1.3078823) q[1];
sx q[1];
rz(-2.8347375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18414149) q[3];
sx q[3];
rz(-1.8098347) q[3];
sx q[3];
rz(-2.4381541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6622322) q[2];
sx q[2];
rz(-1.4630829) q[2];
sx q[2];
rz(-1.4987017) q[2];
rz(-0.44954506) q[3];
sx q[3];
rz(-0.42732626) q[3];
sx q[3];
rz(0.27165616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0736921) q[0];
sx q[0];
rz(-1.3561159) q[0];
sx q[0];
rz(2.7664716) q[0];
rz(0.66315229) q[1];
sx q[1];
rz(-1.8795492) q[1];
sx q[1];
rz(2.1660027) q[1];
rz(1.3781056) q[2];
sx q[2];
rz(-1.4507254) q[2];
sx q[2];
rz(-2.7117827) q[2];
rz(-1.8835486) q[3];
sx q[3];
rz(-1.1328075) q[3];
sx q[3];
rz(-1.9771747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
