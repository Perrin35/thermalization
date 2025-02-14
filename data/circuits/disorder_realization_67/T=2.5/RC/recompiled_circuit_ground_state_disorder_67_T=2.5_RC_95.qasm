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
rz(-2.8400583) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26623959) q[0];
sx q[0];
rz(-1.7325028) q[0];
sx q[0];
rz(-1.7577359) q[0];
rz(-pi) q[1];
rz(-2.1211336) q[2];
sx q[2];
rz(-2.7937903) q[2];
sx q[2];
rz(0.90711275) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1323213) q[1];
sx q[1];
rz(-1.6385211) q[1];
sx q[1];
rz(2.3236815) q[1];
rz(-2.9530617) q[3];
sx q[3];
rz(-1.9664696) q[3];
sx q[3];
rz(-0.71771948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4528759) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(-3.0196043) q[3];
sx q[3];
rz(-0.9361836) q[3];
sx q[3];
rz(3.0854935) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720471) q[0];
sx q[0];
rz(-0.24709728) q[0];
sx q[0];
rz(-0.0023512996) q[0];
rz(-2.7408842) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(0.8561264) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24361783) q[0];
sx q[0];
rz(-1.8254733) q[0];
sx q[0];
rz(-0.93230482) q[0];
x q[1];
rz(-2.0868504) q[2];
sx q[2];
rz(-1.8212089) q[2];
sx q[2];
rz(0.543814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2344831) q[1];
sx q[1];
rz(-0.52965783) q[1];
sx q[1];
rz(1.4717327) q[1];
rz(0.56613381) q[3];
sx q[3];
rz(-1.8834891) q[3];
sx q[3];
rz(1.7305525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.037228435) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(0.47002235) q[2];
rz(2.5445599) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(1.6980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-2.6486588) q[0];
sx q[0];
rz(-2.9135627) q[0];
rz(0.44949624) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(2.1953348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957009) q[0];
sx q[0];
rz(-1.8270703) q[0];
sx q[0];
rz(1.9836224) q[0];
rz(0.7545289) q[2];
sx q[2];
rz(-0.81475329) q[2];
sx q[2];
rz(-2.192564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1552004) q[1];
sx q[1];
rz(-0.3861392) q[1];
sx q[1];
rz(-0.50215747) q[1];
rz(-1.2872996) q[3];
sx q[3];
rz(-1.4186267) q[3];
sx q[3];
rz(0.37934549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9049282) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(-1.4795335) q[2];
rz(-2.9605401) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36412305) q[0];
sx q[0];
rz(-1.5622219) q[0];
sx q[0];
rz(-0.89853483) q[0];
rz(2.6354375) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130044) q[0];
sx q[0];
rz(-2.6046136) q[0];
sx q[0];
rz(3.0130638) q[0];
rz(-2.6674472) q[2];
sx q[2];
rz(-0.80171004) q[2];
sx q[2];
rz(-2.2849883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0702563) q[1];
sx q[1];
rz(-1.500529) q[1];
sx q[1];
rz(1.6855168) q[1];
rz(-pi) q[2];
rz(2.2257099) q[3];
sx q[3];
rz(-1.5228027) q[3];
sx q[3];
rz(-0.97150436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(-1.255704) q[3];
sx q[3];
rz(-1.4016822) q[3];
sx q[3];
rz(3.0953395) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2332377) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(0.18049845) q[0];
rz(1.7620697) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(-2.0464121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9126606) q[0];
sx q[0];
rz(-1.4353509) q[0];
sx q[0];
rz(1.5082466) q[0];
rz(0.081549598) q[2];
sx q[2];
rz(-0.79863859) q[2];
sx q[2];
rz(-0.35432409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68183652) q[1];
sx q[1];
rz(-2.3533258) q[1];
sx q[1];
rz(-0.19758176) q[1];
x q[2];
rz(-1.7322695) q[3];
sx q[3];
rz(-0.69500178) q[3];
sx q[3];
rz(1.7308066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77202648) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(1.0687211) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(-1.2515742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1123947) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(2.8503964) q[0];
rz(-1.4578106) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(-2.2529032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5989953) q[0];
sx q[0];
rz(-0.7805191) q[0];
sx q[0];
rz(0.88135834) q[0];
rz(0.027300731) q[2];
sx q[2];
rz(-0.75746398) q[2];
sx q[2];
rz(0.94099076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7900438) q[1];
sx q[1];
rz(-1.0197395) q[1];
sx q[1];
rz(2.181777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2567886) q[3];
sx q[3];
rz(-1.8591789) q[3];
sx q[3];
rz(2.7183804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7500744) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(-2.1938426) q[2];
rz(0.1968955) q[3];
sx q[3];
rz(-2.9010549) q[3];
sx q[3];
rz(0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31203684) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(0.46992508) q[0];
rz(-1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(1.4039111) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71588281) q[0];
sx q[0];
rz(-2.4615767) q[0];
sx q[0];
rz(-2.143857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9582301) q[2];
sx q[2];
rz(-0.51046267) q[2];
sx q[2];
rz(-1.1092157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39170489) q[1];
sx q[1];
rz(-2.0948695) q[1];
sx q[1];
rz(2.7218372) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2265731) q[3];
sx q[3];
rz(-1.0148409) q[3];
sx q[3];
rz(1.2840301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5002084) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(3.0339962) q[2];
rz(-0.61947668) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-0.25303823) q[0];
rz(-0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(0.94892445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8924361) q[0];
sx q[0];
rz(-1.2741486) q[0];
sx q[0];
rz(-2.9846624) q[0];
x q[1];
rz(-1.1429092) q[2];
sx q[2];
rz(-2.6978931) q[2];
sx q[2];
rz(0.86715172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9946549) q[1];
sx q[1];
rz(-2.6587464) q[1];
sx q[1];
rz(-1.4551244) q[1];
rz(0.83794318) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(-1.6745245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(0.31069791) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0017515) q[0];
sx q[0];
rz(-0.40626353) q[0];
sx q[0];
rz(-1.9061506) q[0];
rz(2.6605117) q[1];
sx q[1];
rz(-0.95293871) q[1];
sx q[1];
rz(1.7135886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9426711) q[0];
sx q[0];
rz(-2.3301972) q[0];
sx q[0];
rz(1.5553724) q[0];
rz(-pi) q[1];
rz(0.88793869) q[2];
sx q[2];
rz(-2.7882034) q[2];
sx q[2];
rz(-1.8977752) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2216827) q[1];
sx q[1];
rz(-2.7836728) q[1];
sx q[1];
rz(2.329925) q[1];
x q[2];
rz(1.8974182) q[3];
sx q[3];
rz(-0.87509586) q[3];
sx q[3];
rz(0.43750924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9674176) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(-0.027912557) q[2];
rz(0.86999718) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(1.1869259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025539909) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(1.0593587) q[0];
rz(1.4147883) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(2.6639604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57741149) q[0];
sx q[0];
rz(-1.7192657) q[0];
sx q[0];
rz(-0.74651697) q[0];
x q[1];
rz(-1.576409) q[2];
sx q[2];
rz(-0.33585784) q[2];
sx q[2];
rz(2.4120021) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.6832215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82179339) q[3];
sx q[3];
rz(-1.104165) q[3];
sx q[3];
rz(0.94576437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9719505) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(-0.36894813) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
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
rz(0.44454642) q[1];
sx q[1];
rz(-2.3311756) q[1];
sx q[1];
rz(-0.12745007) q[1];
rz(2.6981392) q[2];
sx q[2];
rz(-1.2724066) q[2];
sx q[2];
rz(2.3497697) q[2];
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
