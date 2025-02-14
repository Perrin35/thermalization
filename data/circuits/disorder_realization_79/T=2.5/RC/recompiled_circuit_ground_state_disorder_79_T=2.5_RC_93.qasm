OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5918936) q[0];
sx q[0];
rz(-2.0285719) q[0];
sx q[0];
rz(-0.23166238) q[0];
rz(-2.7910233) q[1];
sx q[1];
rz(-2.9313593) q[1];
sx q[1];
rz(-2.3793013) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9462546) q[0];
sx q[0];
rz(-2.1911977) q[0];
sx q[0];
rz(2.0928732) q[0];
rz(0.76114391) q[2];
sx q[2];
rz(-2.5147438) q[2];
sx q[2];
rz(0.055920211) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8232441) q[1];
sx q[1];
rz(-1.6260514) q[1];
sx q[1];
rz(-0.18624185) q[1];
x q[2];
rz(-2.0559156) q[3];
sx q[3];
rz(-1.9645683) q[3];
sx q[3];
rz(-2.6504248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98213696) q[2];
sx q[2];
rz(-2.5348713) q[2];
sx q[2];
rz(2.7098126) q[2];
rz(-2.2744956) q[3];
sx q[3];
rz(-2.1853787) q[3];
sx q[3];
rz(1.506327) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0302439) q[0];
sx q[0];
rz(-0.83006492) q[0];
sx q[0];
rz(1.9507971) q[0];
rz(-1.8973154) q[1];
sx q[1];
rz(-1.1294653) q[1];
sx q[1];
rz(2.5608889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55739072) q[0];
sx q[0];
rz(-1.709849) q[0];
sx q[0];
rz(1.1625421) q[0];
rz(-pi) q[1];
rz(-2.7218616) q[2];
sx q[2];
rz(-1.0786622) q[2];
sx q[2];
rz(-2.555814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3327707) q[1];
sx q[1];
rz(-0.94478196) q[1];
sx q[1];
rz(-1.1762192) q[1];
rz(-pi) q[2];
rz(-2.1070117) q[3];
sx q[3];
rz(-0.76047069) q[3];
sx q[3];
rz(1.5602668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1397436) q[2];
sx q[2];
rz(-1.7364343) q[2];
sx q[2];
rz(-1.639036) q[2];
rz(-0.20770811) q[3];
sx q[3];
rz(-1.3108871) q[3];
sx q[3];
rz(-2.1343855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19596066) q[0];
sx q[0];
rz(-1.0769083) q[0];
sx q[0];
rz(-0.14863241) q[0];
rz(1.0736939) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(0.78280848) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4769094) q[0];
sx q[0];
rz(-2.9626155) q[0];
sx q[0];
rz(-1.8838691) q[0];
rz(-pi) q[1];
x q[1];
rz(0.071657091) q[2];
sx q[2];
rz(-0.53402099) q[2];
sx q[2];
rz(-1.7657745) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69082123) q[1];
sx q[1];
rz(-1.3714002) q[1];
sx q[1];
rz(2.6197207) q[1];
rz(1.820676) q[3];
sx q[3];
rz(-0.80554038) q[3];
sx q[3];
rz(2.4791634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17915501) q[2];
sx q[2];
rz(-2.16733) q[2];
sx q[2];
rz(-1.1243189) q[2];
rz(3.1149241) q[3];
sx q[3];
rz(-2.4494438) q[3];
sx q[3];
rz(-0.28287014) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7836595) q[0];
sx q[0];
rz(-1.1555576) q[0];
sx q[0];
rz(-1.577277) q[0];
rz(1.0464926) q[1];
sx q[1];
rz(-1.0341897) q[1];
sx q[1];
rz(1.1276721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30450102) q[0];
sx q[0];
rz(-0.9795042) q[0];
sx q[0];
rz(-1.558377) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0617969) q[2];
sx q[2];
rz(-1.2376806) q[2];
sx q[2];
rz(0.12399697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8279523) q[1];
sx q[1];
rz(-1.2544187) q[1];
sx q[1];
rz(-2.4183941) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2644482) q[3];
sx q[3];
rz(-0.1538642) q[3];
sx q[3];
rz(2.2779704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6005738) q[2];
sx q[2];
rz(-1.8709196) q[2];
sx q[2];
rz(-0.18826558) q[2];
rz(2.6514734) q[3];
sx q[3];
rz(-1.4508338) q[3];
sx q[3];
rz(-1.274444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1566496) q[0];
sx q[0];
rz(-1.5771414) q[0];
sx q[0];
rz(1.6749325) q[0];
rz(-2.6690392) q[1];
sx q[1];
rz(-0.87635374) q[1];
sx q[1];
rz(-0.98365012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2547823) q[0];
sx q[0];
rz(-1.7521825) q[0];
sx q[0];
rz(2.7104293) q[0];
rz(-2.6783597) q[2];
sx q[2];
rz(-1.7432799) q[2];
sx q[2];
rz(-1.1197937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77351219) q[1];
sx q[1];
rz(-1.0607914) q[1];
sx q[1];
rz(-3.0762005) q[1];
rz(-3.0346924) q[3];
sx q[3];
rz(-1.8351115) q[3];
sx q[3];
rz(-2.3472465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5094362) q[2];
sx q[2];
rz(-1.7551273) q[2];
sx q[2];
rz(0.40718386) q[2];
rz(-1.2716278) q[3];
sx q[3];
rz(-0.41912246) q[3];
sx q[3];
rz(0.65740383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42095175) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(-1.9049013) q[0];
rz(1.1300794) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(2.0225661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7216199) q[0];
sx q[0];
rz(-0.16789745) q[0];
sx q[0];
rz(-0.41746087) q[0];
rz(-pi) q[1];
rz(0.58216146) q[2];
sx q[2];
rz(-2.2560439) q[2];
sx q[2];
rz(-0.2174046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.055370959) q[1];
sx q[1];
rz(-1.1495591) q[1];
sx q[1];
rz(1.7293386) q[1];
rz(-pi) q[2];
rz(0.72596999) q[3];
sx q[3];
rz(-2.5790102) q[3];
sx q[3];
rz(-1.9051304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5469024) q[2];
sx q[2];
rz(-1.12744) q[2];
sx q[2];
rz(2.8315869) q[2];
rz(-1.7075432) q[3];
sx q[3];
rz(-1.6238345) q[3];
sx q[3];
rz(-1.8967418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86842787) q[0];
sx q[0];
rz(-2.7719066) q[0];
sx q[0];
rz(1.0721068) q[0];
rz(1.3605236) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(-1.2881813) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5679703) q[0];
sx q[0];
rz(-2.0814271) q[0];
sx q[0];
rz(0.88167067) q[0];
rz(-2.200618) q[2];
sx q[2];
rz(-0.67747203) q[2];
sx q[2];
rz(1.3969473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2350504) q[1];
sx q[1];
rz(-1.8120011) q[1];
sx q[1];
rz(1.3969087) q[1];
x q[2];
rz(2.6374972) q[3];
sx q[3];
rz(-2.305784) q[3];
sx q[3];
rz(-0.1288165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0897022) q[2];
sx q[2];
rz(-1.6334198) q[2];
sx q[2];
rz(-1.4349597) q[2];
rz(-0.47576225) q[3];
sx q[3];
rz(-2.4475554) q[3];
sx q[3];
rz(2.7580875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7504904) q[0];
sx q[0];
rz(-2.7848211) q[0];
sx q[0];
rz(1.179689) q[0];
rz(-2.7791595) q[1];
sx q[1];
rz(-1.06203) q[1];
sx q[1];
rz(1.5865883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060136) q[0];
sx q[0];
rz(-1.4188674) q[0];
sx q[0];
rz(-2.8169027) q[0];
rz(-2.4612975) q[2];
sx q[2];
rz(-1.1914502) q[2];
sx q[2];
rz(0.85000402) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52380864) q[1];
sx q[1];
rz(-2.0062399) q[1];
sx q[1];
rz(-0.82651386) q[1];
x q[2];
rz(0.15896564) q[3];
sx q[3];
rz(-2.1329125) q[3];
sx q[3];
rz(1.5288439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4156551) q[2];
sx q[2];
rz(-2.026181) q[2];
sx q[2];
rz(-1.1768781) q[2];
rz(-2.7365909) q[3];
sx q[3];
rz(-0.98086762) q[3];
sx q[3];
rz(-0.39308959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28615752) q[0];
sx q[0];
rz(-3.0431008) q[0];
sx q[0];
rz(-1.806102) q[0];
rz(-0.91241854) q[1];
sx q[1];
rz(-1.4211979) q[1];
sx q[1];
rz(3.1013427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7351609) q[0];
sx q[0];
rz(-1.4022151) q[0];
sx q[0];
rz(3.1030036) q[0];
x q[1];
rz(0.98699498) q[2];
sx q[2];
rz(-1.0539471) q[2];
sx q[2];
rz(-1.3911846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1564538) q[1];
sx q[1];
rz(-0.6985526) q[1];
sx q[1];
rz(-0.88628873) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76816948) q[3];
sx q[3];
rz(-1.4171367) q[3];
sx q[3];
rz(0.43936554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5361629) q[2];
sx q[2];
rz(-1.9998113) q[2];
sx q[2];
rz(2.6119168) q[2];
rz(-1.4179432) q[3];
sx q[3];
rz(-0.46897408) q[3];
sx q[3];
rz(2.6403707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1390851) q[0];
sx q[0];
rz(-1.3676099) q[0];
sx q[0];
rz(-3.1097155) q[0];
rz(1.2079976) q[1];
sx q[1];
rz(-0.54787689) q[1];
sx q[1];
rz(0.30778232) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986711) q[0];
sx q[0];
rz(-1.0772155) q[0];
sx q[0];
rz(0.30360766) q[0];
x q[1];
rz(-1.8372907) q[2];
sx q[2];
rz(-1.4478168) q[2];
sx q[2];
rz(0.83713712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5468502) q[1];
sx q[1];
rz(-1.062927) q[1];
sx q[1];
rz(-2.2617729) q[1];
x q[2];
rz(-0.73460893) q[3];
sx q[3];
rz(-2.9501868) q[3];
sx q[3];
rz(-0.020343971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6098183) q[2];
sx q[2];
rz(-2.4419624) q[2];
sx q[2];
rz(0.043665234) q[2];
rz(1.4788491) q[3];
sx q[3];
rz(-0.50373977) q[3];
sx q[3];
rz(1.3879205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8762348) q[0];
sx q[0];
rz(-1.6725412) q[0];
sx q[0];
rz(-1.5681736) q[0];
rz(0.028989446) q[1];
sx q[1];
rz(-2.034076) q[1];
sx q[1];
rz(-2.1709002) q[1];
rz(2.0736442) q[2];
sx q[2];
rz(-2.89123) q[2];
sx q[2];
rz(2.8889026) q[2];
rz(2.0652931) q[3];
sx q[3];
rz(-2.2702379) q[3];
sx q[3];
rz(-0.092580295) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
