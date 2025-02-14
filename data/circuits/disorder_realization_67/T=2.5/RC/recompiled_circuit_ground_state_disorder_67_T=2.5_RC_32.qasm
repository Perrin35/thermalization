OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8118892) q[0];
sx q[0];
rz(5.9726297) q[0];
sx q[0];
rz(7.4818727) q[0];
rz(1.1341473) q[1];
sx q[1];
rz(-2.3448047) q[1];
sx q[1];
rz(-0.30153433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3350007) q[0];
sx q[0];
rz(-1.7552688) q[0];
sx q[0];
rz(-0.16452275) q[0];
x q[1];
rz(2.9542175) q[2];
sx q[2];
rz(-1.2760386) q[2];
sx q[2];
rz(0.32887884) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1323213) q[1];
sx q[1];
rz(-1.5030716) q[1];
sx q[1];
rz(0.81791116) q[1];
rz(-1.1490378) q[3];
sx q[3];
rz(-0.43614498) q[3];
sx q[3];
rz(-0.25806043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(-3.0196043) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(0.056099135) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695456) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(3.1392414) q[0];
rz(0.40070847) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(0.8561264) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1416505) q[0];
sx q[0];
rz(-2.460833) q[0];
sx q[0];
rz(1.9826243) q[0];
rz(-pi) q[1];
rz(-0.28601525) q[2];
sx q[2];
rz(-2.0692424) q[2];
sx q[2];
rz(1.1666544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0217738) q[1];
sx q[1];
rz(-2.0975862) q[1];
sx q[1];
rz(-0.057838126) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5991531) q[3];
sx q[3];
rz(-2.503201) q[3];
sx q[3];
rz(0.29104376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1043642) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(0.47002235) q[2];
rz(0.59703279) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3700579) q[0];
sx q[0];
rz(-2.6486588) q[0];
sx q[0];
rz(2.9135627) q[0];
rz(-0.44949624) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(-0.94625783) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3920604) q[0];
sx q[0];
rz(-2.6596053) q[0];
sx q[0];
rz(-2.1493692) q[0];
x q[1];
rz(-2.4837355) q[2];
sx q[2];
rz(-2.0924753) q[2];
sx q[2];
rz(1.1945832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5214649) q[1];
sx q[1];
rz(-1.9072272) q[1];
sx q[1];
rz(-1.7640339) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15839307) q[3];
sx q[3];
rz(-1.2906646) q[3];
sx q[3];
rz(1.9060077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2366644) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(1.6620592) q[2];
rz(0.18105257) q[3];
sx q[3];
rz(-0.79020399) q[3];
sx q[3];
rz(2.2553867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(2.2430578) q[0];
rz(-2.6354375) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285883) q[0];
sx q[0];
rz(-0.53697907) q[0];
sx q[0];
rz(0.12852886) q[0];
rz(0.47414549) q[2];
sx q[2];
rz(-2.3398826) q[2];
sx q[2];
rz(-0.8566044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6502225) q[1];
sx q[1];
rz(-1.6852324) q[1];
sx q[1];
rz(3.070862) q[1];
rz(-0.060486501) q[3];
sx q[3];
rz(-0.91676676) q[3];
sx q[3];
rz(-2.5791283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8755181) q[2];
sx q[2];
rz(-2.6660599) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(-1.255704) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(-3.0953395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2332377) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(-2.9610942) q[0];
rz(1.3795229) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(-2.0464121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22893208) q[0];
sx q[0];
rz(-1.4353509) q[0];
sx q[0];
rz(1.5082466) q[0];
x q[1];
rz(-1.6542475) q[2];
sx q[2];
rz(-0.7755643) q[2];
sx q[2];
rz(-0.47093876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74873456) q[1];
sx q[1];
rz(-1.7104516) q[1];
sx q[1];
rz(0.77844324) q[1];
rz(-pi) q[2];
rz(-2.2593656) q[3];
sx q[3];
rz(-1.467657) q[3];
sx q[3];
rz(-2.8571199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77202648) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(1.0687211) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1123947) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(-0.29119626) q[0];
rz(1.4578106) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(-2.2529032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68316685) q[0];
sx q[0];
rz(-2.1447276) q[0];
sx q[0];
rz(2.5794585) q[0];
x q[1];
rz(-0.027300731) q[2];
sx q[2];
rz(-2.3841287) q[2];
sx q[2];
rz(0.94099076) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7900438) q[1];
sx q[1];
rz(-2.1218532) q[1];
sx q[1];
rz(-0.95981564) q[1];
rz(-pi) q[2];
rz(-1.2567886) q[3];
sx q[3];
rz(-1.8591789) q[3];
sx q[3];
rz(2.7183804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39151829) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(-2.1938426) q[2];
rz(0.1968955) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(-0.31015629) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31203684) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(0.46992508) q[0];
rz(1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(1.7376815) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.2228773) q[0];
sx q[0];
rz(0.97396429) q[0];
rz(-pi) q[1];
x q[1];
rz(1.092488) q[2];
sx q[2];
rz(-1.385139) q[2];
sx q[2];
rz(0.80362475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.95937377) q[1];
sx q[1];
rz(-1.9314018) q[1];
sx q[1];
rz(1.006516) q[1];
rz(-pi) q[2];
rz(-2.6440355) q[3];
sx q[3];
rz(-0.64421875) q[3];
sx q[3];
rz(1.8807172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5002084) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(3.0339962) q[2];
rz(2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951931) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-0.25303823) q[0];
rz(0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(-0.94892445) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2491566) q[0];
sx q[0];
rz(-1.2741486) q[0];
sx q[0];
rz(-0.15693024) q[0];
x q[1];
rz(-0.19472943) q[2];
sx q[2];
rz(-1.1695332) q[2];
sx q[2];
rz(-1.3347609) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6151877) q[1];
sx q[1];
rz(-1.5171836) q[1];
sx q[1];
rz(-2.0508906) q[1];
rz(-pi) q[2];
rz(2.6213245) q[3];
sx q[3];
rz(-1.0661849) q[3];
sx q[3];
rz(0.59719539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(-0.31069791) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(1.0827304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-0.40626353) q[0];
sx q[0];
rz(1.9061506) q[0];
rz(0.48108092) q[1];
sx q[1];
rz(-0.95293871) q[1];
sx q[1];
rz(1.4280041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9202703) q[0];
sx q[0];
rz(-2.3820665) q[0];
sx q[0];
rz(-3.1253472) q[0];
x q[1];
rz(2.253654) q[2];
sx q[2];
rz(-0.35338923) q[2];
sx q[2];
rz(-1.8977752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7117611) q[1];
sx q[1];
rz(-1.8277543) q[1];
sx q[1];
rz(2.8896232) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8974182) q[3];
sx q[3];
rz(-2.2664968) q[3];
sx q[3];
rz(0.43750924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1741751) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(-3.1136801) q[2];
rz(0.86999718) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025539909) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(-2.0822339) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.4780412) q[1];
sx q[1];
rz(0.47763225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896509) q[0];
sx q[0];
rz(-2.383259) q[0];
sx q[0];
rz(-2.9248021) q[0];
rz(-1.9066493) q[2];
sx q[2];
rz(-1.5689465) q[2];
sx q[2];
rz(-2.2950878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79523008) q[1];
sx q[1];
rz(-1.6246965) q[1];
sx q[1];
rz(1.4583711) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82179339) q[3];
sx q[3];
rz(-2.0374277) q[3];
sx q[3];
rz(-0.94576437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9719505) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(0.36894813) q[2];
rz(1.87489) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(-0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030768) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(2.6970462) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(1.8989904) q[2];
sx q[2];
rz(-1.1482308) q[2];
sx q[2];
rz(-2.2238735) q[2];
rz(0.038158554) q[3];
sx q[3];
rz(-0.46783075) q[3];
sx q[3];
rz(-0.086188407) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
