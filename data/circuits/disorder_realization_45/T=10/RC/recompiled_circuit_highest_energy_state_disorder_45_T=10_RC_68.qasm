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
rz(-1.2019914) q[0];
sx q[0];
rz(-2.658598) q[0];
sx q[0];
rz(1.510409) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2886292) q[0];
sx q[0];
rz(-1.4242111) q[0];
sx q[0];
rz(0.96340553) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32663871) q[2];
sx q[2];
rz(-2.5212598) q[2];
sx q[2];
rz(0.49119821) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9485907) q[1];
sx q[1];
rz(-0.99094245) q[1];
sx q[1];
rz(-1.5327044) q[1];
rz(-pi) q[2];
rz(-1.8584941) q[3];
sx q[3];
rz(-1.7904864) q[3];
sx q[3];
rz(-0.9385329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0240747) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-0.5994125) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-2.3960466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9052542) q[0];
sx q[0];
rz(-0.31762341) q[0];
sx q[0];
rz(0.41877094) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7600766) q[2];
sx q[2];
rz(-1.809263) q[2];
sx q[2];
rz(-0.60205215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88953253) q[1];
sx q[1];
rz(-2.4870858) q[1];
sx q[1];
rz(0.26328523) q[1];
rz(0.50314869) q[3];
sx q[3];
rz(-2.1981578) q[3];
sx q[3];
rz(-0.37093758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(-1.7342742) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5582964) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(2.5335675) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(-1.8720522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8931173) q[0];
sx q[0];
rz(-0.60113827) q[0];
sx q[0];
rz(-2.2918743) q[0];
x q[1];
rz(-2.24176) q[2];
sx q[2];
rz(-1.6353893) q[2];
sx q[2];
rz(-3.10499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3976058) q[1];
sx q[1];
rz(-2.3922046) q[1];
sx q[1];
rz(-2.8515062) q[1];
rz(-pi) q[2];
rz(0.57165159) q[3];
sx q[3];
rz(-2.5749675) q[3];
sx q[3];
rz(1.0039312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.515392) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(-0.0083222566) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(2.9279809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(-1.7648765) q[0];
rz(-1.6142913) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462858) q[0];
sx q[0];
rz(-2.5881672) q[0];
sx q[0];
rz(0.85352104) q[0];
rz(0.3600271) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(1.479666) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69428378) q[1];
sx q[1];
rz(-2.1226624) q[1];
sx q[1];
rz(2.0815316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52526955) q[3];
sx q[3];
rz(-2.161918) q[3];
sx q[3];
rz(2.1256465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.3516124) q[2];
rz(-1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(-3.0426262) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(2.6720572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3833475) q[0];
sx q[0];
rz(-1.4329785) q[0];
sx q[0];
rz(0.024017781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37920538) q[2];
sx q[2];
rz(-2.0742886) q[2];
sx q[2];
rz(0.45631726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27891913) q[1];
sx q[1];
rz(-1.9544744) q[1];
sx q[1];
rz(2.5088599) q[1];
x q[2];
rz(2.9088777) q[3];
sx q[3];
rz(-2.5732627) q[3];
sx q[3];
rz(1.1956904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(1.8355231) q[2];
rz(-0.4246873) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(-2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.3816396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2974325) q[0];
sx q[0];
rz(-1.4844262) q[0];
sx q[0];
rz(2.3303836) q[0];
x q[1];
rz(0.53333111) q[2];
sx q[2];
rz(-1.0600277) q[2];
sx q[2];
rz(-3.0549218) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4037212) q[1];
sx q[1];
rz(-2.3624237) q[1];
sx q[1];
rz(2.8761151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9074677) q[3];
sx q[3];
rz(-0.30933274) q[3];
sx q[3];
rz(-2.7721289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7606925) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(-0.5091269) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(0.38207644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(-0.089381889) q[0];
rz(-1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(0.99348974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29879323) q[0];
sx q[0];
rz(-1.5325108) q[0];
sx q[0];
rz(0.78193112) q[0];
rz(-2.2682796) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(2.8755434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9536994) q[1];
sx q[1];
rz(-0.39361289) q[1];
sx q[1];
rz(-1.962349) q[1];
rz(2.089608) q[3];
sx q[3];
rz(-2.1362274) q[3];
sx q[3];
rz(-2.4861002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79954687) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(-0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.8079669) q[0];
rz(-2.2857621) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(2.2241101) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79865341) q[0];
sx q[0];
rz(-2.9173033) q[0];
sx q[0];
rz(2.0464315) q[0];
rz(-0.94709227) q[2];
sx q[2];
rz(-1.4741885) q[2];
sx q[2];
rz(1.4132702) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92261295) q[1];
sx q[1];
rz(-0.41480468) q[1];
sx q[1];
rz(0.99402438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64669426) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(-0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.8514006) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(-0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23583394) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(1.2241036) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.3714429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4215721) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(0.30465845) q[0];
rz(2.0061699) q[2];
sx q[2];
rz(-1.0658385) q[2];
sx q[2];
rz(1.4471042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4491357) q[1];
sx q[1];
rz(-0.68329158) q[1];
sx q[1];
rz(1.2399252) q[1];
rz(2.9670466) q[3];
sx q[3];
rz(-2.2851599) q[3];
sx q[3];
rz(-2.3256231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.203043) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(-2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(1.8800053) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7670583) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(0.64819711) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2637973) q[0];
sx q[0];
rz(-1.5938984) q[0];
sx q[0];
rz(-0.013851555) q[0];
x q[1];
rz(0.83309116) q[2];
sx q[2];
rz(-2.0948185) q[2];
sx q[2];
rz(-1.810871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5773864) q[1];
sx q[1];
rz(-1.5930452) q[1];
sx q[1];
rz(2.4425227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51697124) q[3];
sx q[3];
rz(-1.7423769) q[3];
sx q[3];
rz(-1.7652546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.466922) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(0.35015097) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(2.1351142) q[3];
sx q[3];
rz(-1.1075533) q[3];
sx q[3];
rz(1.5816734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
