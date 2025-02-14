OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6634231) q[0];
sx q[0];
rz(-2.2278251) q[0];
sx q[0];
rz(-1.0650286) q[0];
rz(-1.5364667) q[1];
sx q[1];
rz(-2.0361587) q[1];
sx q[1];
rz(3.0450568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9912579) q[0];
sx q[0];
rz(-1.152395) q[0];
sx q[0];
rz(2.9166448) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41806721) q[2];
sx q[2];
rz(-0.56116784) q[2];
sx q[2];
rz(0.18218606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21743821) q[1];
sx q[1];
rz(-0.70637843) q[1];
sx q[1];
rz(2.6732754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7126585) q[3];
sx q[3];
rz(-0.58462287) q[3];
sx q[3];
rz(-1.7138002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1221293) q[2];
sx q[2];
rz(-0.18522842) q[2];
sx q[2];
rz(1.8786001) q[2];
rz(0.39877912) q[3];
sx q[3];
rz(-1.140927) q[3];
sx q[3];
rz(-1.0293452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1102092) q[0];
sx q[0];
rz(-2.5094014) q[0];
sx q[0];
rz(-1.6803918) q[0];
rz(3.0714495) q[1];
sx q[1];
rz(-1.5488011) q[1];
sx q[1];
rz(-0.4313012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9330641) q[0];
sx q[0];
rz(-1.5613436) q[0];
sx q[0];
rz(1.509311) q[0];
x q[1];
rz(0.29745086) q[2];
sx q[2];
rz(-1.1701726) q[2];
sx q[2];
rz(1.144415) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.166242) q[1];
sx q[1];
rz(-1.0570475) q[1];
sx q[1];
rz(-0.789289) q[1];
rz(2.7106695) q[3];
sx q[3];
rz(-0.049073372) q[3];
sx q[3];
rz(1.2653093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9554319) q[2];
sx q[2];
rz(-1.7308851) q[2];
sx q[2];
rz(-0.53039941) q[2];
rz(2.0576599) q[3];
sx q[3];
rz(-1.089774) q[3];
sx q[3];
rz(-1.2497905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1029866) q[0];
sx q[0];
rz(-0.57612053) q[0];
sx q[0];
rz(-1.2023793) q[0];
rz(2.1309958) q[1];
sx q[1];
rz(-1.8162411) q[1];
sx q[1];
rz(0.030166322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3836244) q[0];
sx q[0];
rz(-0.080648184) q[0];
sx q[0];
rz(2.5190802) q[0];
x q[1];
rz(-1.6232477) q[2];
sx q[2];
rz(-1.5546335) q[2];
sx q[2];
rz(2.9132622) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2304014) q[1];
sx q[1];
rz(-1.4433858) q[1];
sx q[1];
rz(-0.34613737) q[1];
rz(-pi) q[2];
rz(1.7272337) q[3];
sx q[3];
rz(-0.57830252) q[3];
sx q[3];
rz(2.3012637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61115894) q[2];
sx q[2];
rz(-2.9035089) q[2];
sx q[2];
rz(0.37453026) q[2];
rz(1.1206867) q[3];
sx q[3];
rz(-1.4846385) q[3];
sx q[3];
rz(-0.69168958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.197914) q[0];
sx q[0];
rz(-1.8489842) q[0];
sx q[0];
rz(-1.9669272) q[0];
rz(-1.8265751) q[1];
sx q[1];
rz(-2.0349793) q[1];
sx q[1];
rz(0.16474251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4602866) q[0];
sx q[0];
rz(-2.3577655) q[0];
sx q[0];
rz(0.71390217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4061178) q[2];
sx q[2];
rz(-2.0782172) q[2];
sx q[2];
rz(-2.8230482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.03962) q[1];
sx q[1];
rz(-2.0594739) q[1];
sx q[1];
rz(0.38445977) q[1];
rz(-pi) q[2];
rz(-3.1026671) q[3];
sx q[3];
rz(-1.4994326) q[3];
sx q[3];
rz(1.4369278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69079196) q[2];
sx q[2];
rz(-1.6054634) q[2];
sx q[2];
rz(2.6353321) q[2];
rz(2.6311503) q[3];
sx q[3];
rz(-2.3531395) q[3];
sx q[3];
rz(2.3300664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6798169) q[0];
sx q[0];
rz(-0.33648574) q[0];
sx q[0];
rz(2.3315954) q[0];
rz(0.098377146) q[1];
sx q[1];
rz(-2.7521303) q[1];
sx q[1];
rz(-2.817165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8194991) q[0];
sx q[0];
rz(-3.0190912) q[0];
sx q[0];
rz(1.5550645) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7067028) q[2];
sx q[2];
rz(-1.7909484) q[2];
sx q[2];
rz(0.94583257) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.788841) q[1];
sx q[1];
rz(-0.60539027) q[1];
sx q[1];
rz(2.8245138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39119519) q[3];
sx q[3];
rz(-2.6325691) q[3];
sx q[3];
rz(2.5768821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35679945) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(-2.5830833) q[2];
rz(1.2950581) q[3];
sx q[3];
rz(-1.7241155) q[3];
sx q[3];
rz(-0.3524802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.9117987) q[0];
sx q[0];
rz(-2.7138382) q[0];
sx q[0];
rz(1.7708923) q[0];
rz(-1.3937344) q[1];
sx q[1];
rz(-1.0153898) q[1];
sx q[1];
rz(0.14399354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8407366) q[0];
sx q[0];
rz(-1.6413781) q[0];
sx q[0];
rz(-1.6750122) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8990895) q[2];
sx q[2];
rz(-1.1272301) q[2];
sx q[2];
rz(-2.7672845) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.6274144) q[1];
sx q[1];
rz(-2.5251221) q[1];
sx q[1];
rz(1.7656839) q[1];
rz(-0.64231974) q[3];
sx q[3];
rz(-1.2313456) q[3];
sx q[3];
rz(-2.385115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7250942) q[2];
sx q[2];
rz(-2.7855253) q[2];
sx q[2];
rz(0.98102513) q[2];
rz(-0.52870005) q[3];
sx q[3];
rz(-1.7873535) q[3];
sx q[3];
rz(-2.5813812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.8037146) q[0];
sx q[0];
rz(-0.80334544) q[0];
sx q[0];
rz(-1.2038318) q[0];
rz(-3.1375258) q[1];
sx q[1];
rz(-1.7151567) q[1];
sx q[1];
rz(0.34122658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2409397) q[0];
sx q[0];
rz(-1.9918858) q[0];
sx q[0];
rz(-2.3809041) q[0];
rz(-pi) q[1];
rz(3.1401743) q[2];
sx q[2];
rz(-2.2621415) q[2];
sx q[2];
rz(-1.6306016) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9014185) q[1];
sx q[1];
rz(-1.56159) q[1];
sx q[1];
rz(-2.9926278) q[1];
rz(-1.4248203) q[3];
sx q[3];
rz(-1.2314704) q[3];
sx q[3];
rz(1.4682352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5228806) q[2];
sx q[2];
rz(-1.2423542) q[2];
sx q[2];
rz(2.8704571) q[2];
rz(1.6297657) q[3];
sx q[3];
rz(-0.96704331) q[3];
sx q[3];
rz(-0.40201521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790344) q[0];
sx q[0];
rz(-2.5738578) q[0];
sx q[0];
rz(-2.973279) q[0];
rz(2.1315101) q[1];
sx q[1];
rz(-1.1770266) q[1];
sx q[1];
rz(-3.041306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71609488) q[0];
sx q[0];
rz(-1.8779813) q[0];
sx q[0];
rz(1.2745122) q[0];
rz(-pi) q[1];
rz(-2.9434154) q[2];
sx q[2];
rz(-0.74286425) q[2];
sx q[2];
rz(-2.4695549) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6339019) q[1];
sx q[1];
rz(-1.7724438) q[1];
sx q[1];
rz(1.545791) q[1];
x q[2];
rz(2.4504205) q[3];
sx q[3];
rz(-0.6414957) q[3];
sx q[3];
rz(-2.3506899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68891305) q[2];
sx q[2];
rz(-2.6426688) q[2];
sx q[2];
rz(-1.4208687) q[2];
rz(-1.0639327) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(-3.0498144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7559779) q[0];
sx q[0];
rz(-1.4571964) q[0];
sx q[0];
rz(-3.1047367) q[0];
rz(1.256975) q[1];
sx q[1];
rz(-1.6391862) q[1];
sx q[1];
rz(-1.3704376) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4697594) q[0];
sx q[0];
rz(-1.2984833) q[0];
sx q[0];
rz(0.16118523) q[0];
rz(-1.7212935) q[2];
sx q[2];
rz(-1.9353107) q[2];
sx q[2];
rz(0.30752674) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3210813) q[1];
sx q[1];
rz(-0.91225636) q[1];
sx q[1];
rz(0.35595004) q[1];
x q[2];
rz(-1.4716205) q[3];
sx q[3];
rz(-0.41891019) q[3];
sx q[3];
rz(0.3813972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2799985) q[2];
sx q[2];
rz(-1.7656606) q[2];
sx q[2];
rz(-0.17246788) q[2];
rz(-0.55781588) q[3];
sx q[3];
rz(-2.6004801) q[3];
sx q[3];
rz(-1.7012168) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43011343) q[0];
sx q[0];
rz(-2.5028296) q[0];
sx q[0];
rz(2.8016256) q[0];
rz(-2.4760447) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(-1.8284214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42178828) q[0];
sx q[0];
rz(-1.7504842) q[0];
sx q[0];
rz(1.522904) q[0];
rz(-pi) q[1];
rz(-2.8808241) q[2];
sx q[2];
rz(-1.4968902) q[2];
sx q[2];
rz(-2.7735965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39496702) q[1];
sx q[1];
rz(-1.2685742) q[1];
sx q[1];
rz(0.24924596) q[1];
rz(-0.35452133) q[3];
sx q[3];
rz(-0.42663048) q[3];
sx q[3];
rz(-1.420295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7790935) q[2];
sx q[2];
rz(-2.7069147) q[2];
sx q[2];
rz(0.67081007) q[2];
rz(-0.32495156) q[3];
sx q[3];
rz(-0.81348014) q[3];
sx q[3];
rz(1.0130829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83298564) q[0];
sx q[0];
rz(-1.8600464) q[0];
sx q[0];
rz(-1.6639584) q[0];
rz(2.1113405) q[1];
sx q[1];
rz(-1.4719084) q[1];
sx q[1];
rz(1.7869064) q[1];
rz(0.75922913) q[2];
sx q[2];
rz(-2.4674131) q[2];
sx q[2];
rz(-2.7684733) q[2];
rz(-0.53976157) q[3];
sx q[3];
rz(-0.65315795) q[3];
sx q[3];
rz(-1.5140189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
