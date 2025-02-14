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
rz(2.6296339) q[0];
sx q[0];
rz(-0.40357885) q[0];
sx q[0];
rz(-3.1374748) q[0];
rz(0.68139684) q[1];
sx q[1];
rz(1.8713142) q[1];
sx q[1];
rz(11.125578) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487605) q[0];
sx q[0];
rz(-1.3997243) q[0];
sx q[0];
rz(1.2890512) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3902293) q[2];
sx q[2];
rz(-1.3168502) q[2];
sx q[2];
rz(0.51515173) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3057249) q[1];
sx q[1];
rz(-2.222371) q[1];
sx q[1];
rz(2.0209794) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.639569) q[3];
sx q[3];
rz(-2.9675238) q[3];
sx q[3];
rz(-0.93910142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19930856) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(2.671177) q[2];
rz(-0.89620245) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50297058) q[0];
sx q[0];
rz(-0.48415411) q[0];
sx q[0];
rz(0.36889398) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.8577441) q[1];
sx q[1];
rz(0.2643815) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1702829) q[0];
sx q[0];
rz(-1.7242887) q[0];
sx q[0];
rz(0.92043368) q[0];
rz(2.5892841) q[2];
sx q[2];
rz(-2.8495516) q[2];
sx q[2];
rz(2.8949182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0708587) q[1];
sx q[1];
rz(-1.2675335) q[1];
sx q[1];
rz(1.8766848) q[1];
x q[2];
rz(1.4503497) q[3];
sx q[3];
rz(-2.3494738) q[3];
sx q[3];
rz(1.9223546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5276864) q[2];
sx q[2];
rz(-2.9759585) q[2];
sx q[2];
rz(1.8370834) q[2];
rz(-0.3197318) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(-2.6379697) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9673135) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(2.018003) q[0];
rz(0.38791052) q[1];
sx q[1];
rz(-1.265641) q[1];
sx q[1];
rz(0.33043114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39836568) q[0];
sx q[0];
rz(-1.8332229) q[0];
sx q[0];
rz(2.9578379) q[0];
rz(-0.13339146) q[2];
sx q[2];
rz(-0.61163227) q[2];
sx q[2];
rz(0.69452121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6092005) q[1];
sx q[1];
rz(-0.48621854) q[1];
sx q[1];
rz(2.4052252) q[1];
rz(-pi) q[2];
rz(2.4027545) q[3];
sx q[3];
rz(-2.3519302) q[3];
sx q[3];
rz(1.9146321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.419751) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(1.9640131) q[2];
rz(1.2490595) q[3];
sx q[3];
rz(-1.8355337) q[3];
sx q[3];
rz(-0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881507) q[0];
sx q[0];
rz(-1.4295239) q[0];
sx q[0];
rz(-3.0528659) q[0];
rz(0.42452043) q[1];
sx q[1];
rz(-2.4219234) q[1];
sx q[1];
rz(-1.7321865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78414375) q[0];
sx q[0];
rz(-1.7567092) q[0];
sx q[0];
rz(-1.2307412) q[0];
x q[1];
rz(0.97169496) q[2];
sx q[2];
rz(-1.1568152) q[2];
sx q[2];
rz(1.5091015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.7551491) q[1];
sx q[1];
rz(-2.1514822) q[1];
sx q[1];
rz(1.1183075) q[1];
rz(-pi) q[2];
rz(2.8187784) q[3];
sx q[3];
rz(-1.0720616) q[3];
sx q[3];
rz(-0.43969595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7523664) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-0.43238315) q[2];
rz(0.45807517) q[3];
sx q[3];
rz(-1.4832486) q[3];
sx q[3];
rz(-1.0079916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78419375) q[0];
sx q[0];
rz(-2.9965239) q[0];
sx q[0];
rz(0.31132895) q[0];
rz(1.1773102) q[1];
sx q[1];
rz(-1.9639587) q[1];
sx q[1];
rz(-2.6572773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949849) q[0];
sx q[0];
rz(-1.4123045) q[0];
sx q[0];
rz(-0.096317795) q[0];
rz(-1.3179146) q[2];
sx q[2];
rz(-0.30159471) q[2];
sx q[2];
rz(0.67267928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97515177) q[1];
sx q[1];
rz(-1.7130392) q[1];
sx q[1];
rz(-1.1935545) q[1];
x q[2];
rz(2.43552) q[3];
sx q[3];
rz(-0.19565249) q[3];
sx q[3];
rz(2.5062989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3207265) q[2];
sx q[2];
rz(-1.7118688) q[2];
sx q[2];
rz(0.42496625) q[2];
rz(0.10415569) q[3];
sx q[3];
rz(-0.36915532) q[3];
sx q[3];
rz(0.66515508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7729618) q[0];
sx q[0];
rz(-2.5017128) q[0];
sx q[0];
rz(-2.9848918) q[0];
rz(2.8796097) q[1];
sx q[1];
rz(-1.9888839) q[1];
sx q[1];
rz(1.6798457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3647385) q[0];
sx q[0];
rz(-0.97524736) q[0];
sx q[0];
rz(0.2100581) q[0];
rz(-pi) q[1];
rz(0.22905519) q[2];
sx q[2];
rz(-0.88729268) q[2];
sx q[2];
rz(-0.4313658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3586629) q[1];
sx q[1];
rz(-1.6921669) q[1];
sx q[1];
rz(-0.18106457) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1739985) q[3];
sx q[3];
rz(-2.2809416) q[3];
sx q[3];
rz(2.0098734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67019749) q[2];
sx q[2];
rz(-1.4963701) q[2];
sx q[2];
rz(-2.4246598) q[2];
rz(2.9251621) q[3];
sx q[3];
rz(-1.8559034) q[3];
sx q[3];
rz(-1.1855679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4099429) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(1.2662079) q[0];
rz(2.3833497) q[1];
sx q[1];
rz(-1.1781324) q[1];
sx q[1];
rz(-1.0479124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78572014) q[0];
sx q[0];
rz(-1.8685172) q[0];
sx q[0];
rz(-0.96426086) q[0];
rz(-pi) q[1];
rz(-2.2323147) q[2];
sx q[2];
rz(-0.66864852) q[2];
sx q[2];
rz(2.3956211) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3057237) q[1];
sx q[1];
rz(-2.7717238) q[1];
sx q[1];
rz(1.7643516) q[1];
rz(-0.47509681) q[3];
sx q[3];
rz(-0.49834278) q[3];
sx q[3];
rz(1.4628177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3332112) q[2];
sx q[2];
rz(-1.0071249) q[2];
sx q[2];
rz(3.132931) q[2];
rz(1.0041142) q[3];
sx q[3];
rz(-0.4929556) q[3];
sx q[3];
rz(1.0021817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.778331) q[0];
sx q[0];
rz(-2.9982428) q[0];
sx q[0];
rz(-2.1891201) q[0];
rz(-1.5337503) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(2.105377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0506315) q[0];
sx q[0];
rz(-1.0337752) q[0];
sx q[0];
rz(2.0799252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2042734) q[2];
sx q[2];
rz(-0.96666217) q[2];
sx q[2];
rz(1.095739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1324824) q[1];
sx q[1];
rz(-0.81071893) q[1];
sx q[1];
rz(-2.2894163) q[1];
rz(0.72662205) q[3];
sx q[3];
rz(-1.22261) q[3];
sx q[3];
rz(-2.6451056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2967534) q[2];
sx q[2];
rz(-0.71594683) q[2];
sx q[2];
rz(2.1412264) q[2];
rz(1.865546) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100383) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(-2.8874183) q[0];
rz(-1.3379478) q[1];
sx q[1];
rz(-2.281052) q[1];
sx q[1];
rz(0.77802229) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36150071) q[0];
sx q[0];
rz(-0.92467148) q[0];
sx q[0];
rz(-1.9801937) q[0];
rz(-pi) q[1];
rz(-0.67963375) q[2];
sx q[2];
rz(-0.25654116) q[2];
sx q[2];
rz(2.7367965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.48315037) q[1];
sx q[1];
rz(-1.1643693) q[1];
sx q[1];
rz(2.5220925) q[1];
rz(1.7860402) q[3];
sx q[3];
rz(-1.6674893) q[3];
sx q[3];
rz(1.7413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.079387) q[2];
sx q[2];
rz(-1.6239245) q[2];
sx q[2];
rz(0.13266955) q[2];
rz(2.0343871) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(-1.6975105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4001813) q[0];
sx q[0];
rz(-1.1962698) q[0];
sx q[0];
rz(-0.74914002) q[0];
rz(2.6465042) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(-1.7503768) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86319727) q[0];
sx q[0];
rz(-0.94587284) q[0];
sx q[0];
rz(1.004659) q[0];
rz(2.3074564) q[2];
sx q[2];
rz(-0.5668525) q[2];
sx q[2];
rz(-1.8994272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0616845) q[1];
sx q[1];
rz(-2.7322953) q[1];
sx q[1];
rz(-0.10959919) q[1];
x q[2];
rz(-0.81328525) q[3];
sx q[3];
rz(-0.67768598) q[3];
sx q[3];
rz(-0.71302524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64080015) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(0.87727171) q[2];
rz(-2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40182879) q[0];
sx q[0];
rz(-2.5943828) q[0];
sx q[0];
rz(1.7964969) q[0];
rz(1.8632035) q[1];
sx q[1];
rz(-0.31247333) q[1];
sx q[1];
rz(-0.023275274) q[1];
rz(1.5398434) q[2];
sx q[2];
rz(-0.94018117) q[2];
sx q[2];
rz(-1.2311819) q[2];
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
