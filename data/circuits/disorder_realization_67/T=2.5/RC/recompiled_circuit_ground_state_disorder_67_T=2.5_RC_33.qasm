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
rz(-0.31055561) q[0];
sx q[0];
rz(-1.9429053) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(0.30153433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3350007) q[0];
sx q[0];
rz(-1.7552688) q[0];
sx q[0];
rz(0.16452275) q[0];
rz(-pi) q[1];
rz(2.9542175) q[2];
sx q[2];
rz(-1.865554) q[2];
sx q[2];
rz(2.8127138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1323213) q[1];
sx q[1];
rz(-1.6385211) q[1];
sx q[1];
rz(-2.3236815) q[1];
x q[2];
rz(1.1687247) q[3];
sx q[3];
rz(-1.744606) q[3];
sx q[3];
rz(-2.2151092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4528759) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(2.0067046) q[2];
rz(0.12198837) q[3];
sx q[3];
rz(-0.9361836) q[3];
sx q[3];
rz(3.0854935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695456) q[0];
sx q[0];
rz(-0.24709728) q[0];
sx q[0];
rz(0.0023512996) q[0];
rz(-0.40070847) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(-2.2854663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8979748) q[0];
sx q[0];
rz(-1.3161193) q[0];
sx q[0];
rz(-0.93230482) q[0];
x q[1];
rz(2.0868504) q[2];
sx q[2];
rz(-1.8212089) q[2];
sx q[2];
rz(2.5977787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90710956) q[1];
sx q[1];
rz(-2.6119348) q[1];
sx q[1];
rz(-1.6698599) q[1];
rz(-pi) q[2];
rz(-2.5754588) q[3];
sx q[3];
rz(-1.8834891) q[3];
sx q[3];
rz(1.7305525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1043642) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(-0.59703279) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(-1.6980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(0.22802995) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(-0.94625783) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11433521) q[0];
sx q[0];
rz(-1.9693621) q[0];
sx q[0];
rz(2.8629659) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1990132) q[2];
sx q[2];
rz(-1.0120856) q[2];
sx q[2];
rz(-3.1329346) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0277532) q[1];
sx q[1];
rz(-1.7530754) q[1];
sx q[1];
rz(-0.34232847) q[1];
rz(-pi) q[2];
rz(-2.9831996) q[3];
sx q[3];
rz(-1.2906646) q[3];
sx q[3];
rz(1.235585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9049282) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(1.4795335) q[2];
rz(-0.18105257) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(-0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(-2.2430578) q[0];
rz(-0.50615519) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(-1.4599962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.579297) q[0];
sx q[0];
rz(-2.102872) q[0];
sx q[0];
rz(1.4946372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0115667) q[2];
sx q[2];
rz(-0.87723644) q[2];
sx q[2];
rz(0.22089411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0467028) q[1];
sx q[1];
rz(-3.0071435) q[1];
sx q[1];
rz(2.1220757) q[1];
rz(-pi) q[2];
rz(-1.4921032) q[3];
sx q[3];
rz(-2.4851812) q[3];
sx q[3];
rz(2.479913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(1.6391222) q[2];
rz(1.255704) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(3.0953395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083549) q[0];
sx q[0];
rz(-0.73035705) q[0];
sx q[0];
rz(-2.9610942) q[0];
rz(-1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(1.0951805) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22893208) q[0];
sx q[0];
rz(-1.7062418) q[0];
sx q[0];
rz(1.6333461) q[0];
x q[1];
rz(-1.4873452) q[2];
sx q[2];
rz(-0.7755643) q[2];
sx q[2];
rz(-2.6706539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95848786) q[1];
sx q[1];
rz(-2.339683) q[1];
sx q[1];
rz(1.7657178) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4093231) q[3];
sx q[3];
rz(-0.69500178) q[3];
sx q[3];
rz(1.7308066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77202648) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(-1.0687211) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0291979) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(-2.8503964) q[0];
rz(-1.4578106) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(0.88868946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55808961) q[0];
sx q[0];
rz(-2.0348685) q[0];
sx q[0];
rz(-2.2232673) q[0];
x q[1];
rz(-2.3843147) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(-0.64964408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7900438) q[1];
sx q[1];
rz(-2.1218532) q[1];
sx q[1];
rz(-0.95981564) q[1];
rz(-1.2567886) q[3];
sx q[3];
rz(-1.8591789) q[3];
sx q[3];
rz(-0.42321229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7500744) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(-0.94775003) q[2];
rz(-2.9446972) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(2.8314364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31203684) q[0];
sx q[0];
rz(-1.1126248) q[0];
sx q[0];
rz(-0.46992508) q[0];
rz(-1.6795109) q[1];
sx q[1];
rz(-1.8406248) q[1];
sx q[1];
rz(-1.7376815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71588281) q[0];
sx q[0];
rz(-0.68001594) q[0];
sx q[0];
rz(-2.143857) q[0];
rz(-pi) q[1];
rz(1.1833625) q[2];
sx q[2];
rz(-0.51046267) q[2];
sx q[2];
rz(2.032377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39170489) q[1];
sx q[1];
rz(-1.0467231) q[1];
sx q[1];
rz(-0.4197555) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49755712) q[3];
sx q[3];
rz(-0.64421875) q[3];
sx q[3];
rz(-1.2608755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64138428) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(-0.10759648) q[2];
rz(2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951931) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-2.8885544) q[0];
rz(-3.053275) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(2.1926682) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7530156) q[0];
sx q[0];
rz(-0.33451329) q[0];
sx q[0];
rz(1.0981513) q[0];
rz(-pi) q[1];
rz(1.1429092) q[2];
sx q[2];
rz(-2.6978931) q[2];
sx q[2];
rz(2.2744409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.526405) q[1];
sx q[1];
rz(-1.5171836) q[1];
sx q[1];
rz(-1.0907021) q[1];
rz(-2.3036495) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(-1.6745245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(-0.31069791) q[2];
rz(0.24142309) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(1.235442) q[0];
rz(-2.6605117) q[1];
sx q[1];
rz(-0.95293871) q[1];
sx q[1];
rz(1.4280041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7803376) q[0];
sx q[0];
rz(-1.5819823) q[0];
sx q[0];
rz(-0.7594603) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9129006) q[2];
sx q[2];
rz(-1.2989825) q[2];
sx q[2];
rz(-0.52955176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9199099) q[1];
sx q[1];
rz(-2.7836728) q[1];
sx q[1];
rz(0.81166761) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4190524) q[3];
sx q[3];
rz(-1.3219537) q[3];
sx q[3];
rz(-1.3470847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(3.1136801) q[2];
rz(0.86999718) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025539909) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(-1.0593587) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.4780412) q[1];
sx q[1];
rz(-2.6639604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5641812) q[0];
sx q[0];
rz(-1.7192657) q[0];
sx q[0];
rz(-0.74651697) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.576409) q[2];
sx q[2];
rz(-2.8057348) q[2];
sx q[2];
rz(0.72959057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79523008) q[1];
sx q[1];
rz(-1.5168961) q[1];
sx q[1];
rz(-1.4583711) q[1];
x q[2];
rz(-2.5390901) q[3];
sx q[3];
rz(-0.91711603) q[3];
sx q[3];
rz(-0.22887558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9719505) q[2];
sx q[2];
rz(-0.98196882) q[2];
sx q[2];
rz(2.7726445) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(2.1784311) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030768) q[0];
sx q[0];
rz(-1.2110447) q[0];
sx q[0];
rz(-1.3883653) q[0];
rz(2.6970462) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(0.44345345) q[2];
sx q[2];
rz(-1.8691861) q[2];
sx q[2];
rz(-0.791823) q[2];
rz(0.46753771) q[3];
sx q[3];
rz(-1.553592) q[3];
sx q[3];
rz(-1.6910465) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
