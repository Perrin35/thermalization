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
rz(-1.7614814) q[0];
sx q[0];
rz(7.8742134) q[0];
sx q[0];
rz(5.9013517) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5474553) q[0];
sx q[0];
rz(-1.0411123) q[0];
sx q[0];
rz(-2.269777) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.051048267) q[2];
sx q[2];
rz(-2.7709392) q[2];
sx q[2];
rz(-1.6459344) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3081835) q[1];
sx q[1];
rz(-1.7004968) q[1];
sx q[1];
rz(0.42713487) q[1];
x q[2];
rz(2.1006868) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.4640946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(-1.5748242) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(0.69275698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67398706) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.1760733) q[0];
rz(-0.067497079) q[1];
sx q[1];
rz(-0.57756966) q[1];
sx q[1];
rz(-0.31879058) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8628937) q[0];
sx q[0];
rz(-2.3434533) q[0];
sx q[0];
rz(-0.44095914) q[0];
rz(-0.484157) q[2];
sx q[2];
rz(-2.6645497) q[2];
sx q[2];
rz(-2.8271528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76004878) q[1];
sx q[1];
rz(-2.4439619) q[1];
sx q[1];
rz(0.95894869) q[1];
rz(-1.8570921) q[3];
sx q[3];
rz(-1.4501374) q[3];
sx q[3];
rz(0.21721043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(0.47743615) q[2];
rz(1.7834974) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(-1.1392449) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(2.4635945) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292873) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(-1.4112006) q[0];
x q[1];
rz(0.94495541) q[2];
sx q[2];
rz(-1.2765108) q[2];
sx q[2];
rz(-1.4465894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4809847) q[1];
sx q[1];
rz(-1.0502195) q[1];
sx q[1];
rz(2.1064227) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0150681) q[3];
sx q[3];
rz(-0.58150269) q[3];
sx q[3];
rz(-0.27669981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.5072301) q[2];
sx q[2];
rz(-2.0646084) q[2];
rz(1.8814603) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-2.8857359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28427163) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(-0.76048365) q[0];
rz(-0.35176945) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(0.15028353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1635123) q[0];
sx q[0];
rz(-0.0038359782) q[0];
sx q[0];
rz(1.0413961) q[0];
rz(-1.9530832) q[2];
sx q[2];
rz(-1.4103544) q[2];
sx q[2];
rz(0.40942243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80881883) q[1];
sx q[1];
rz(-1.9017694) q[1];
sx q[1];
rz(-0.95167758) q[1];
x q[2];
rz(-2.6500449) q[3];
sx q[3];
rz(-0.50954362) q[3];
sx q[3];
rz(-2.5423637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1111697) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-0.58491582) q[2];
rz(-0.065718204) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(-2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718662) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(2.4140893) q[0];
rz(1.8456521) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(2.4268699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22435221) q[0];
sx q[0];
rz(-0.2324129) q[0];
sx q[0];
rz(-0.31347855) q[0];
rz(0.3008414) q[2];
sx q[2];
rz(-1.46035) q[2];
sx q[2];
rz(-0.12953239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2421869) q[1];
sx q[1];
rz(-1.1922298) q[1];
sx q[1];
rz(0.19902163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0017702) q[3];
sx q[3];
rz(-2.4845893) q[3];
sx q[3];
rz(-0.13323076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456445) q[2];
sx q[2];
rz(-1.5278634) q[2];
sx q[2];
rz(1.9963416) q[2];
rz(0.20259419) q[3];
sx q[3];
rz(-0.069772094) q[3];
sx q[3];
rz(-0.32874671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458493) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(2.5766418) q[1];
sx q[1];
rz(-2.0940557) q[1];
sx q[1];
rz(-1.3289183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.57083) q[0];
sx q[0];
rz(-1.098363) q[0];
sx q[0];
rz(0.60586141) q[0];
rz(-pi) q[1];
rz(-0.33514993) q[2];
sx q[2];
rz(-2.9932635) q[2];
sx q[2];
rz(0.89491612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79949483) q[1];
sx q[1];
rz(-1.3159385) q[1];
sx q[1];
rz(1.8138769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9181913) q[3];
sx q[3];
rz(-1.3376682) q[3];
sx q[3];
rz(-2.4883771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38587511) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(0.34783777) q[2];
rz(0.31473413) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9921853) q[0];
sx q[0];
rz(-1.6353761) q[0];
sx q[0];
rz(-0.95112479) q[0];
rz(0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(0.94400418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2601449) q[0];
sx q[0];
rz(-1.2627739) q[0];
sx q[0];
rz(-0.016168895) q[0];
rz(-pi) q[1];
rz(-1.0877785) q[2];
sx q[2];
rz(-1.1406787) q[2];
sx q[2];
rz(-1.6295691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2846825) q[1];
sx q[1];
rz(-1.3976685) q[1];
sx q[1];
rz(-0.022927479) q[1];
x q[2];
rz(1.0349501) q[3];
sx q[3];
rz(-2.6831045) q[3];
sx q[3];
rz(-1.637984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(1.7000807) q[2];
rz(3.1169543) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(2.1520481) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80855075) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(-3.0035875) q[0];
rz(-1.0778614) q[1];
sx q[1];
rz(-1.4638008) q[1];
sx q[1];
rz(0.93625751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2168343) q[0];
sx q[0];
rz(-0.48473293) q[0];
sx q[0];
rz(-0.68883552) q[0];
rz(-0.8291275) q[2];
sx q[2];
rz(-1.1489604) q[2];
sx q[2];
rz(1.2527996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.028658) q[1];
sx q[1];
rz(-1.3880265) q[1];
sx q[1];
rz(0.014928047) q[1];
x q[2];
rz(-2.4560087) q[3];
sx q[3];
rz(-1.778025) q[3];
sx q[3];
rz(2.384546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1428947) q[2];
sx q[2];
rz(-1.2724718) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(3.0190492) q[3];
sx q[3];
rz(-0.69869852) q[3];
sx q[3];
rz(2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97813022) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(-0.36460707) q[0];
rz(1.7568024) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(-2.2273831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62637882) q[0];
sx q[0];
rz(-2.6914586) q[0];
sx q[0];
rz(1.2594957) q[0];
rz(-pi) q[1];
rz(1.6169102) q[2];
sx q[2];
rz(-0.83674016) q[2];
sx q[2];
rz(-3.0457599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57371512) q[1];
sx q[1];
rz(-2.2063631) q[1];
sx q[1];
rz(-1.741841) q[1];
rz(-pi) q[2];
rz(1.7201603) q[3];
sx q[3];
rz(-1.0348667) q[3];
sx q[3];
rz(-1.5045297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8387575) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(-1.6287104) q[2];
rz(2.5527939) q[3];
sx q[3];
rz(-1.9353119) q[3];
sx q[3];
rz(-0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673679) q[0];
sx q[0];
rz(-2.5741757) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(0.65471634) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(-2.7213352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3493372) q[0];
sx q[0];
rz(-0.47264034) q[0];
sx q[0];
rz(-2.0562754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7905529) q[2];
sx q[2];
rz(-1.9091064) q[2];
sx q[2];
rz(-1.4684791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6984562) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(3.0522947) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5758366) q[3];
sx q[3];
rz(-1.2671736) q[3];
sx q[3];
rz(0.88241258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(2.1743656) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(-0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16348542) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(-0.12374395) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(-1.4370357) q[2];
sx q[2];
rz(-1.8885713) q[2];
sx q[2];
rz(1.0071913) q[2];
rz(2.8295201) q[3];
sx q[3];
rz(-1.5371116) q[3];
sx q[3];
rz(-0.1826382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
