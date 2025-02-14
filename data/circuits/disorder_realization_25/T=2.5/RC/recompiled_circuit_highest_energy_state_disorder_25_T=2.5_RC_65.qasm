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
rz(0.79431835) q[0];
rz(0.89214832) q[1];
sx q[1];
rz(4.6214784) q[1];
sx q[1];
rz(5.9393039) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958069) q[0];
sx q[0];
rz(-2.8749879) q[0];
sx q[0];
rz(-2.8050799) q[0];
rz(1.0353824) q[2];
sx q[2];
rz(-2.0119925) q[2];
sx q[2];
rz(0.23278415) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3229494) q[1];
sx q[1];
rz(-1.7293849) q[1];
sx q[1];
rz(-1.701638) q[1];
rz(3.0214693) q[3];
sx q[3];
rz(-1.7645482) q[3];
sx q[3];
rz(1.7156977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6227365) q[2];
sx q[2];
rz(-1.7665665) q[2];
sx q[2];
rz(2.1905017) q[2];
rz(2.2226492) q[3];
sx q[3];
rz(-2.2008342) q[3];
sx q[3];
rz(-2.9520891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6491991) q[0];
sx q[0];
rz(-2.352582) q[0];
sx q[0];
rz(1.5570109) q[0];
rz(2.4015265) q[1];
sx q[1];
rz(-0.84741712) q[1];
sx q[1];
rz(1.7795631) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82845485) q[0];
sx q[0];
rz(-2.065593) q[0];
sx q[0];
rz(-0.54101701) q[0];
rz(-pi) q[1];
rz(-1.9744423) q[2];
sx q[2];
rz(-1.4440618) q[2];
sx q[2];
rz(0.53467804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3950065) q[1];
sx q[1];
rz(-0.94330245) q[1];
sx q[1];
rz(-1.1203517) q[1];
rz(0.091097761) q[3];
sx q[3];
rz(-0.86074867) q[3];
sx q[3];
rz(1.8820845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7035199) q[2];
sx q[2];
rz(-1.1601996) q[2];
sx q[2];
rz(-2.4959219) q[2];
rz(0.034363184) q[3];
sx q[3];
rz(-0.86302257) q[3];
sx q[3];
rz(0.063095108) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768169) q[0];
sx q[0];
rz(-0.21012981) q[0];
sx q[0];
rz(2.3576376) q[0];
rz(0.3166554) q[1];
sx q[1];
rz(-1.4953934) q[1];
sx q[1];
rz(1.0139611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72958845) q[0];
sx q[0];
rz(-2.9656898) q[0];
sx q[0];
rz(0.094502016) q[0];
rz(-pi) q[1];
rz(1.0929669) q[2];
sx q[2];
rz(-1.3690557) q[2];
sx q[2];
rz(0.7283177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0590208) q[1];
sx q[1];
rz(-1.3157651) q[1];
sx q[1];
rz(1.8262556) q[1];
rz(2.2714642) q[3];
sx q[3];
rz(-2.4149272) q[3];
sx q[3];
rz(0.18584968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97518146) q[2];
sx q[2];
rz(-1.3813436) q[2];
sx q[2];
rz(-2.735801) q[2];
rz(0.25519145) q[3];
sx q[3];
rz(-1.8813335) q[3];
sx q[3];
rz(-0.7757799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.880421) q[0];
sx q[0];
rz(-2.4242327) q[0];
sx q[0];
rz(2.1521211) q[0];
rz(0.62670341) q[1];
sx q[1];
rz(-2.6040405) q[1];
sx q[1];
rz(0.42868838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31590474) q[0];
sx q[0];
rz(-1.4576212) q[0];
sx q[0];
rz(2.3772888) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4312885) q[2];
sx q[2];
rz(-1.9006471) q[2];
sx q[2];
rz(-2.863186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7871514) q[1];
sx q[1];
rz(-2.0729182) q[1];
sx q[1];
rz(0.67014613) q[1];
rz(-0.71487244) q[3];
sx q[3];
rz(-1.7809516) q[3];
sx q[3];
rz(-1.6870267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2428212) q[2];
sx q[2];
rz(-3.0148744) q[2];
sx q[2];
rz(0.80647331) q[2];
rz(1.9857231) q[3];
sx q[3];
rz(-2.3838796) q[3];
sx q[3];
rz(2.9466467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72215801) q[0];
sx q[0];
rz(-2.2787155) q[0];
sx q[0];
rz(0.3062329) q[0];
rz(-2.560794) q[1];
sx q[1];
rz(-2.2598946) q[1];
sx q[1];
rz(0.64540783) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5991678) q[0];
sx q[0];
rz(-1.732848) q[0];
sx q[0];
rz(2.5697744) q[0];
x q[1];
rz(1.0694765) q[2];
sx q[2];
rz(-1.638661) q[2];
sx q[2];
rz(2.9508638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20504643) q[1];
sx q[1];
rz(-1.0433739) q[1];
sx q[1];
rz(-0.54103367) q[1];
rz(-pi) q[2];
rz(-1.801412) q[3];
sx q[3];
rz(-1.5531887) q[3];
sx q[3];
rz(1.0215774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8530276) q[2];
sx q[2];
rz(-1.3753076) q[2];
sx q[2];
rz(-2.1993401) q[2];
rz(0.45091584) q[3];
sx q[3];
rz(-0.98603606) q[3];
sx q[3];
rz(-2.8773785) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577724) q[0];
sx q[0];
rz(-1.8696996) q[0];
sx q[0];
rz(0.94495946) q[0];
rz(-0.68929976) q[1];
sx q[1];
rz(-1.1079451) q[1];
sx q[1];
rz(1.1030997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24352989) q[0];
sx q[0];
rz(-1.8193428) q[0];
sx q[0];
rz(0.72406475) q[0];
x q[1];
rz(-0.1774474) q[2];
sx q[2];
rz(-1.1523231) q[2];
sx q[2];
rz(-0.22289101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0660432) q[1];
sx q[1];
rz(-0.72490059) q[1];
sx q[1];
rz(-2.3182475) q[1];
rz(0.96007993) q[3];
sx q[3];
rz(-1.999049) q[3];
sx q[3];
rz(2.6564889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61392203) q[2];
sx q[2];
rz(-0.066378243) q[2];
sx q[2];
rz(1.684368) q[2];
rz(-2.1704604) q[3];
sx q[3];
rz(-1.6077653) q[3];
sx q[3];
rz(-3.0253313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3986452) q[0];
sx q[0];
rz(-1.255641) q[0];
sx q[0];
rz(2.9435112) q[0];
rz(-0.48175851) q[1];
sx q[1];
rz(-1.878673) q[1];
sx q[1];
rz(1.6162965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4264674) q[0];
sx q[0];
rz(-1.5670527) q[0];
sx q[0];
rz(0.0088028839) q[0];
rz(-0.003633066) q[2];
sx q[2];
rz(-2.0272519) q[2];
sx q[2];
rz(-0.21286392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7864466) q[1];
sx q[1];
rz(-1.0736559) q[1];
sx q[1];
rz(-2.0694147) q[1];
rz(-pi) q[2];
rz(2.0644998) q[3];
sx q[3];
rz(-1.6564661) q[3];
sx q[3];
rz(0.44891294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5109479) q[2];
sx q[2];
rz(-1.7787361) q[2];
sx q[2];
rz(0.93713588) q[2];
rz(-0.42738327) q[3];
sx q[3];
rz(-2.1338227) q[3];
sx q[3];
rz(-1.5905323) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9464924) q[0];
sx q[0];
rz(-1.4668523) q[0];
sx q[0];
rz(-1.7505769) q[0];
rz(1.4637949) q[1];
sx q[1];
rz(-2.548806) q[1];
sx q[1];
rz(-2.6189199) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36230817) q[0];
sx q[0];
rz(-1.8982419) q[0];
sx q[0];
rz(0.069652005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9982469) q[2];
sx q[2];
rz(-2.0269577) q[2];
sx q[2];
rz(0.83269115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9620888) q[1];
sx q[1];
rz(-2.8887862) q[1];
sx q[1];
rz(2.4180555) q[1];
rz(-pi) q[2];
rz(2.0955047) q[3];
sx q[3];
rz(-1.7874092) q[3];
sx q[3];
rz(-0.10817402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3576144) q[2];
sx q[2];
rz(-1.7366333) q[2];
sx q[2];
rz(-0.092255175) q[2];
rz(-2.9914894) q[3];
sx q[3];
rz(-1.9893407) q[3];
sx q[3];
rz(0.86103719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1152231) q[0];
sx q[0];
rz(-2.8404591) q[0];
sx q[0];
rz(-1.2979771) q[0];
rz(2.2932529) q[1];
sx q[1];
rz(-1.6611165) q[1];
sx q[1];
rz(2.8678105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1116347) q[0];
sx q[0];
rz(-2.9024501) q[0];
sx q[0];
rz(2.9947623) q[0];
rz(2.8086814) q[2];
sx q[2];
rz(-2.430418) q[2];
sx q[2];
rz(-0.43363562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9628512) q[1];
sx q[1];
rz(-0.76707572) q[1];
sx q[1];
rz(3.0050468) q[1];
x q[2];
rz(-1.3630834) q[3];
sx q[3];
rz(-1.150032) q[3];
sx q[3];
rz(1.8150235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0724642) q[2];
sx q[2];
rz(-1.0651361) q[2];
sx q[2];
rz(0.96159846) q[2];
rz(-1.0570863) q[3];
sx q[3];
rz(-1.0870533) q[3];
sx q[3];
rz(-2.6915153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3707651) q[0];
sx q[0];
rz(-1.2130883) q[0];
sx q[0];
rz(0.83592498) q[0];
rz(-1.4113374) q[1];
sx q[1];
rz(-2.5524499) q[1];
sx q[1];
rz(-0.020847598) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3086714) q[0];
sx q[0];
rz(-0.14946689) q[0];
sx q[0];
rz(-1.368379) q[0];
rz(1.4336939) q[2];
sx q[2];
rz(-2.0404271) q[2];
sx q[2];
rz(0.048222311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9350773) q[1];
sx q[1];
rz(-0.40136101) q[1];
sx q[1];
rz(0.727833) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18414149) q[3];
sx q[3];
rz(-1.8098347) q[3];
sx q[3];
rz(0.70343859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6622322) q[2];
sx q[2];
rz(-1.4630829) q[2];
sx q[2];
rz(-1.642891) q[2];
rz(2.6920476) q[3];
sx q[3];
rz(-2.7142664) q[3];
sx q[3];
rz(-0.27165616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0679006) q[0];
sx q[0];
rz(-1.3561159) q[0];
sx q[0];
rz(2.7664716) q[0];
rz(0.66315229) q[1];
sx q[1];
rz(-1.8795492) q[1];
sx q[1];
rz(2.1660027) q[1];
rz(-1.0085874) q[2];
sx q[2];
rz(-2.9149483) q[2];
sx q[2];
rz(-1.6915773) q[2];
rz(-1.258044) q[3];
sx q[3];
rz(-2.0087852) q[3];
sx q[3];
rz(1.164418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
