OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(1.2553296) q[0];
sx q[0];
rz(9.9160739) q[0];
rz(-3.4499912) q[1];
sx q[1];
rz(1.8412794) q[1];
sx q[1];
rz(12.507764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97293905) q[0];
sx q[0];
rz(-2.0117451) q[0];
sx q[0];
rz(-0.95885528) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6858133) q[2];
sx q[2];
rz(-2.3037065) q[2];
sx q[2];
rz(-1.3742336) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3112009) q[1];
sx q[1];
rz(-1.4102542) q[1];
sx q[1];
rz(2.83648) q[1];
x q[2];
rz(0.079777579) q[3];
sx q[3];
rz(-1.2523796) q[3];
sx q[3];
rz(1.7169881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0356174) q[2];
sx q[2];
rz(-2.3507037) q[2];
sx q[2];
rz(-2.7604575) q[2];
rz(3.1086339) q[3];
sx q[3];
rz(-1.4573174) q[3];
sx q[3];
rz(1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9829262) q[0];
sx q[0];
rz(-0.52678147) q[0];
sx q[0];
rz(2.7016933) q[0];
rz(0.061821763) q[1];
sx q[1];
rz(-1.5301751) q[1];
sx q[1];
rz(2.8672112) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8500811) q[0];
sx q[0];
rz(-0.78118229) q[0];
sx q[0];
rz(-2.2456332) q[0];
rz(-pi) q[1];
rz(-1.1438649) q[2];
sx q[2];
rz(-1.7795658) q[2];
sx q[2];
rz(-0.17886272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.741863) q[1];
sx q[1];
rz(-1.6464087) q[1];
sx q[1];
rz(2.0194598) q[1];
x q[2];
rz(2.5952789) q[3];
sx q[3];
rz(-1.7726745) q[3];
sx q[3];
rz(3.0936732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2939833) q[2];
sx q[2];
rz(-0.01394883) q[2];
sx q[2];
rz(0.070276109) q[2];
rz(-0.62486068) q[3];
sx q[3];
rz(-1.3288386) q[3];
sx q[3];
rz(3.0666053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4570419) q[0];
sx q[0];
rz(-1.6383189) q[0];
sx q[0];
rz(-0.3918089) q[0];
rz(0.99569544) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(0.044274274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366695) q[0];
sx q[0];
rz(-1.2630487) q[0];
sx q[0];
rz(-0.97881628) q[0];
rz(-0.70719126) q[2];
sx q[2];
rz(-0.91333616) q[2];
sx q[2];
rz(0.75188504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7024732) q[1];
sx q[1];
rz(-2.8608192) q[1];
sx q[1];
rz(0.72661119) q[1];
rz(-pi) q[2];
rz(1.6674394) q[3];
sx q[3];
rz(-2.1657888) q[3];
sx q[3];
rz(-0.71221029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5135045) q[2];
sx q[2];
rz(-0.91614437) q[2];
sx q[2];
rz(-1.8872895) q[2];
rz(3.0301376) q[3];
sx q[3];
rz(-2.1696551) q[3];
sx q[3];
rz(2.286262) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3067538) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(2.2656031) q[0];
rz(-0.39729473) q[1];
sx q[1];
rz(-1.5700424) q[1];
sx q[1];
rz(-1.3321336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080058424) q[0];
sx q[0];
rz(-1.3803015) q[0];
sx q[0];
rz(-0.47866042) q[0];
rz(-pi) q[1];
rz(0.45355637) q[2];
sx q[2];
rz(-2.1906914) q[2];
sx q[2];
rz(0.93728055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.024272763) q[1];
sx q[1];
rz(-1.7242715) q[1];
sx q[1];
rz(1.6502817) q[1];
rz(-pi) q[2];
rz(-0.2355742) q[3];
sx q[3];
rz(-1.0456418) q[3];
sx q[3];
rz(0.88334879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1490271) q[2];
sx q[2];
rz(-1.6529447) q[2];
sx q[2];
rz(-1.3943025) q[2];
rz(-2.1483138) q[3];
sx q[3];
rz(-1.1634049) q[3];
sx q[3];
rz(-1.8786028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(-2.2827523) q[0];
rz(2.7186588) q[1];
sx q[1];
rz(-2.6507288) q[1];
sx q[1];
rz(1.4609969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92991023) q[0];
sx q[0];
rz(-1.0722463) q[0];
sx q[0];
rz(-1.7105667) q[0];
rz(-pi) q[1];
rz(-2.016045) q[2];
sx q[2];
rz(-1.5684874) q[2];
sx q[2];
rz(0.18773676) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8940058) q[1];
sx q[1];
rz(-1.3787706) q[1];
sx q[1];
rz(-0.2094261) q[1];
rz(0.14270881) q[3];
sx q[3];
rz(-1.0694711) q[3];
sx q[3];
rz(-0.10979788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4451663) q[2];
sx q[2];
rz(-2.7327765) q[2];
sx q[2];
rz(1.6430631) q[2];
rz(2.0131352) q[3];
sx q[3];
rz(-1.107629) q[3];
sx q[3];
rz(-2.4988153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6978825) q[0];
sx q[0];
rz(-1.8118129) q[0];
sx q[0];
rz(-0.71769303) q[0];
rz(2.7806661) q[1];
sx q[1];
rz(-0.91054994) q[1];
sx q[1];
rz(0.31401971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234834) q[0];
sx q[0];
rz(-0.49719329) q[0];
sx q[0];
rz(0.086693356) q[0];
x q[1];
rz(1.6144361) q[2];
sx q[2];
rz(-0.72481643) q[2];
sx q[2];
rz(2.5284655) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2870658) q[1];
sx q[1];
rz(-0.57537938) q[1];
sx q[1];
rz(0.96436149) q[1];
rz(-pi) q[2];
rz(0.76992294) q[3];
sx q[3];
rz(-1.439538) q[3];
sx q[3];
rz(-0.24410393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6641984) q[2];
sx q[2];
rz(-1.3496642) q[2];
sx q[2];
rz(-2.9273709) q[2];
rz(-0.91841206) q[3];
sx q[3];
rz(-2.1954506) q[3];
sx q[3];
rz(-1.4001747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65900954) q[0];
sx q[0];
rz(-0.56203401) q[0];
sx q[0];
rz(0.62579489) q[0];
rz(-1.91045) q[1];
sx q[1];
rz(-1.8742671) q[1];
sx q[1];
rz(2.321718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5441262) q[0];
sx q[0];
rz(-1.0080999) q[0];
sx q[0];
rz(-1.94348) q[0];
rz(2.0076114) q[2];
sx q[2];
rz(-0.90825082) q[2];
sx q[2];
rz(2.8134741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4266127) q[1];
sx q[1];
rz(-2.8016114) q[1];
sx q[1];
rz(1.3424385) q[1];
rz(0.19015892) q[3];
sx q[3];
rz(-3.0453322) q[3];
sx q[3];
rz(-3.0031693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95255533) q[2];
sx q[2];
rz(-1.2574715) q[2];
sx q[2];
rz(2.7541568) q[2];
rz(-2.6881325) q[3];
sx q[3];
rz(-0.87415868) q[3];
sx q[3];
rz(2.9816154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74359918) q[0];
sx q[0];
rz(-2.0026119) q[0];
sx q[0];
rz(-2.0177662) q[0];
rz(0.75346142) q[1];
sx q[1];
rz(-1.9517784) q[1];
sx q[1];
rz(1.168728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92896398) q[0];
sx q[0];
rz(-0.88200404) q[0];
sx q[0];
rz(0.11879158) q[0];
x q[1];
rz(0.0059466023) q[2];
sx q[2];
rz(-1.1134673) q[2];
sx q[2];
rz(-1.3471239) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.035188474) q[1];
sx q[1];
rz(-0.6995844) q[1];
sx q[1];
rz(-2.0856218) q[1];
rz(2.6791689) q[3];
sx q[3];
rz(-2.1743757) q[3];
sx q[3];
rz(-1.0360595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1184825) q[2];
sx q[2];
rz(-1.0175984) q[2];
sx q[2];
rz(-2.6712096) q[2];
rz(-1.6810301) q[3];
sx q[3];
rz(-1.2602256) q[3];
sx q[3];
rz(-2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7505782) q[0];
sx q[0];
rz(-1.7711201) q[0];
sx q[0];
rz(1.4960666) q[0];
rz(-2.0059313) q[1];
sx q[1];
rz(-1.8541502) q[1];
sx q[1];
rz(-2.6527203) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.785381) q[0];
sx q[0];
rz(-1.4309034) q[0];
sx q[0];
rz(0.032462151) q[0];
rz(-1.7411918) q[2];
sx q[2];
rz(-2.6238447) q[2];
sx q[2];
rz(-2.2108159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1355537) q[1];
sx q[1];
rz(-0.054243739) q[1];
sx q[1];
rz(-2.7179621) q[1];
rz(3.0472894) q[3];
sx q[3];
rz(-0.62178388) q[3];
sx q[3];
rz(2.0981385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44498542) q[2];
sx q[2];
rz(-0.61034909) q[2];
sx q[2];
rz(-2.0390017) q[2];
rz(-1.6164814) q[3];
sx q[3];
rz(-0.81086719) q[3];
sx q[3];
rz(-1.471224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0082323) q[0];
sx q[0];
rz(-2.6120549) q[0];
sx q[0];
rz(1.3326921) q[0];
rz(2.7545199) q[1];
sx q[1];
rz(-1.8862855) q[1];
sx q[1];
rz(-2.0852087) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1833216) q[0];
sx q[0];
rz(-2.859451) q[0];
sx q[0];
rz(-1.0484379) q[0];
x q[1];
rz(-1.6956639) q[2];
sx q[2];
rz(-1.9859295) q[2];
sx q[2];
rz(-1.268569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.38195189) q[1];
sx q[1];
rz(-2.3084062) q[1];
sx q[1];
rz(0.52565378) q[1];
rz(-2.4293618) q[3];
sx q[3];
rz(-0.81347403) q[3];
sx q[3];
rz(1.2016736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1308088) q[2];
sx q[2];
rz(-2.2937412) q[2];
sx q[2];
rz(1.7424142) q[2];
rz(-1.0731267) q[3];
sx q[3];
rz(-1.9173887) q[3];
sx q[3];
rz(0.46260241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.90049967) q[0];
sx q[0];
rz(-1.4971965) q[0];
sx q[0];
rz(1.5370488) q[0];
rz(-0.72147876) q[1];
sx q[1];
rz(-2.571048) q[1];
sx q[1];
rz(-1.2068988) q[1];
rz(-2.4084694) q[2];
sx q[2];
rz(-2.5903268) q[2];
sx q[2];
rz(-3.0916284) q[2];
rz(3.0400757) q[3];
sx q[3];
rz(-2.3675167) q[3];
sx q[3];
rz(1.7270796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
