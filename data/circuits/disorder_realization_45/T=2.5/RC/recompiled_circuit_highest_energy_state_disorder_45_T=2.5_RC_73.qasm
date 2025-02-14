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
rz(-0.51195872) q[0];
sx q[0];
rz(6.6867642) q[0];
sx q[0];
rz(9.4288958) q[0];
rz(-2.4601958) q[1];
sx q[1];
rz(-1.8713142) q[1];
sx q[1];
rz(-1.4407925) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6928322) q[0];
sx q[0];
rz(-1.3997243) q[0];
sx q[0];
rz(-1.8525415) q[0];
rz(-pi) q[1];
rz(1.2295159) q[2];
sx q[2];
rz(-2.2926039) q[2];
sx q[2];
rz(2.316458) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5087252) q[1];
sx q[1];
rz(-0.77287107) q[1];
sx q[1];
rz(-0.51846026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.639569) q[3];
sx q[3];
rz(-2.9675238) q[3];
sx q[3];
rz(0.93910142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19930856) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(0.47041565) q[2];
rz(-2.2453902) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(-1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.50297058) q[0];
sx q[0];
rz(-0.48415411) q[0];
sx q[0];
rz(-0.36889398) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.2838485) q[1];
sx q[1];
rz(2.8772112) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97130972) q[0];
sx q[0];
rz(-1.7242887) q[0];
sx q[0];
rz(0.92043368) q[0];
rz(-2.8910341) q[2];
sx q[2];
rz(-1.7224285) q[2];
sx q[2];
rz(2.3506468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7355613) q[1];
sx q[1];
rz(-1.8623061) q[1];
sx q[1];
rz(0.31707615) q[1];
rz(-pi) q[2];
rz(0.78231298) q[3];
sx q[3];
rz(-1.4851598) q[3];
sx q[3];
rz(-0.2667564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.6139063) q[2];
sx q[2];
rz(-0.1656342) q[2];
sx q[2];
rz(1.8370834) q[2];
rz(0.3197318) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(-0.50362292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17427915) q[0];
sx q[0];
rz(-2.9496084) q[0];
sx q[0];
rz(-1.1235896) q[0];
rz(-0.38791052) q[1];
sx q[1];
rz(-1.8759517) q[1];
sx q[1];
rz(0.33043114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0200121) q[0];
sx q[0];
rz(-0.3191443) q[0];
sx q[0];
rz(0.97346755) q[0];
rz(-pi) q[1];
rz(0.60744384) q[2];
sx q[2];
rz(-1.4943549) q[2];
sx q[2];
rz(0.7668524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3302686) q[1];
sx q[1];
rz(-1.9243334) q[1];
sx q[1];
rz(-1.2297022) q[1];
rz(-pi) q[2];
rz(2.4027545) q[3];
sx q[3];
rz(-0.78966245) q[3];
sx q[3];
rz(-1.9146321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.419751) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(1.9640131) q[2];
rz(1.8925331) q[3];
sx q[3];
rz(-1.8355337) q[3];
sx q[3];
rz(-3.1031928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053442001) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(-3.0528659) q[0];
rz(-0.42452043) q[1];
sx q[1];
rz(-0.71966925) q[1];
sx q[1];
rz(1.4094062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30520327) q[0];
sx q[0];
rz(-2.7557748) q[0];
sx q[0];
rz(2.0842537) q[0];
x q[1];
rz(-0.48894162) q[2];
sx q[2];
rz(-1.0283111) q[2];
sx q[2];
rz(-0.20636339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3864435) q[1];
sx q[1];
rz(-0.99011043) q[1];
sx q[1];
rz(-2.0232852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0921246) q[3];
sx q[3];
rz(-1.2884669) q[3];
sx q[3];
rz(-1.8518333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38922629) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-0.43238315) q[2];
rz(-2.6835175) q[3];
sx q[3];
rz(-1.658344) q[3];
sx q[3];
rz(-2.1336011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78419375) q[0];
sx q[0];
rz(-0.14506871) q[0];
sx q[0];
rz(2.8302637) q[0];
rz(-1.9642824) q[1];
sx q[1];
rz(-1.177634) q[1];
sx q[1];
rz(-0.48431531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949849) q[0];
sx q[0];
rz(-1.7292882) q[0];
sx q[0];
rz(-3.0452749) q[0];
x q[1];
rz(-1.863345) q[2];
sx q[2];
rz(-1.496409) q[2];
sx q[2];
rz(1.1400346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1664409) q[1];
sx q[1];
rz(-1.7130392) q[1];
sx q[1];
rz(1.1935545) q[1];
rz(-pi) q[2];
rz(-0.70607263) q[3];
sx q[3];
rz(-2.9459402) q[3];
sx q[3];
rz(-2.5062989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3207265) q[2];
sx q[2];
rz(-1.4297239) q[2];
sx q[2];
rz(-2.7166264) q[2];
rz(-0.10415569) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(0.66515508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7729618) q[0];
sx q[0];
rz(-0.63987982) q[0];
sx q[0];
rz(-2.9848918) q[0];
rz(-2.8796097) q[1];
sx q[1];
rz(-1.9888839) q[1];
sx q[1];
rz(-1.6798457) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7279434) q[0];
sx q[0];
rz(-2.5143412) q[0];
sx q[0];
rz(1.2722737) q[0];
x q[1];
rz(0.22905519) q[2];
sx q[2];
rz(-2.2543) q[2];
sx q[2];
rz(-2.7102269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81002767) q[1];
sx q[1];
rz(-1.3910783) q[1];
sx q[1];
rz(1.6941638) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3912638) q[3];
sx q[3];
rz(-1.273386) q[3];
sx q[3];
rz(-0.70575037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4713952) q[2];
sx q[2];
rz(-1.6452226) q[2];
sx q[2];
rz(2.4246598) q[2];
rz(0.21643058) q[3];
sx q[3];
rz(-1.2856893) q[3];
sx q[3];
rz(1.9560248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7316498) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(-1.8753847) q[0];
rz(-0.75824291) q[1];
sx q[1];
rz(-1.9634602) q[1];
sx q[1];
rz(1.0479124) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3558725) q[0];
sx q[0];
rz(-1.2730755) q[0];
sx q[0];
rz(2.1773318) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0133466) q[2];
sx q[2];
rz(-1.9614904) q[2];
sx q[2];
rz(1.3732131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0430773) q[1];
sx q[1];
rz(-1.2081573) q[1];
sx q[1];
rz(-3.0671544) q[1];
rz(1.326845) q[3];
sx q[3];
rz(-2.0097964) q[3];
sx q[3];
rz(-1.148996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3332112) q[2];
sx q[2];
rz(-2.1344678) q[2];
sx q[2];
rz(-3.132931) q[2];
rz(1.0041142) q[3];
sx q[3];
rz(-0.4929556) q[3];
sx q[3];
rz(-2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36326161) q[0];
sx q[0];
rz(-2.9982428) q[0];
sx q[0];
rz(-2.1891201) q[0];
rz(1.5337503) q[1];
sx q[1];
rz(-1.2865103) q[1];
sx q[1];
rz(-1.0362157) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8995953) q[0];
sx q[0];
rz(-1.138666) q[0];
sx q[0];
rz(2.5431387) q[0];
x q[1];
rz(1.9373193) q[2];
sx q[2];
rz(-0.96666217) q[2];
sx q[2];
rz(1.095739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10401512) q[1];
sx q[1];
rz(-1.0733782) q[1];
sx q[1];
rz(0.90105547) q[1];
x q[2];
rz(-2.4149706) q[3];
sx q[3];
rz(-1.9189827) q[3];
sx q[3];
rz(-0.49648703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2967534) q[2];
sx q[2];
rz(-2.4256458) q[2];
sx q[2];
rz(-2.1412264) q[2];
rz(1.2760466) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(1.067151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.100383) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(-0.25417438) q[0];
rz(1.8036448) q[1];
sx q[1];
rz(-0.86054069) q[1];
sx q[1];
rz(-0.77802229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26296534) q[0];
sx q[0];
rz(-2.3926982) q[0];
sx q[0];
rz(-0.48567943) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4073952) q[2];
sx q[2];
rz(-1.7694574) q[2];
sx q[2];
rz(-0.29115788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6584423) q[1];
sx q[1];
rz(-1.1643693) q[1];
sx q[1];
rz(0.61950018) q[1];
x q[2];
rz(1.7860402) q[3];
sx q[3];
rz(-1.6674893) q[3];
sx q[3];
rz(-1.4002864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.079387) q[2];
sx q[2];
rz(-1.6239245) q[2];
sx q[2];
rz(-3.0089231) q[2];
rz(-1.1072055) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(-1.6975105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4001813) q[0];
sx q[0];
rz(-1.9453229) q[0];
sx q[0];
rz(-2.3924526) q[0];
rz(0.49508849) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(-1.3912158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86319727) q[0];
sx q[0];
rz(-2.1957198) q[0];
sx q[0];
rz(2.1369336) q[0];
rz(-2.3074564) q[2];
sx q[2];
rz(-2.5747402) q[2];
sx q[2];
rz(1.2421654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5500945) q[1];
sx q[1];
rz(-1.6143394) q[1];
sx q[1];
rz(-2.7344879) q[1];
rz(-pi) q[2];
rz(1.0416609) q[3];
sx q[3];
rz(-1.1253998) q[3];
sx q[3];
rz(2.9187849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5007925) q[2];
sx q[2];
rz(-1.6668789) q[2];
sx q[2];
rz(2.2643209) q[2];
rz(2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(-1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40182879) q[0];
sx q[0];
rz(-2.5943828) q[0];
sx q[0];
rz(1.7964969) q[0];
rz(-1.8632035) q[1];
sx q[1];
rz(-2.8291193) q[1];
sx q[1];
rz(3.1183174) q[1];
rz(1.6017492) q[2];
sx q[2];
rz(-2.2014115) q[2];
sx q[2];
rz(1.9104107) q[2];
rz(1.18289) q[3];
sx q[3];
rz(-1.9961052) q[3];
sx q[3];
rz(0.75710184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
