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
rz(1.1062082) q[0];
sx q[0];
rz(-2.455403) q[0];
sx q[0];
rz(-1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(-2.6931813) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657239) q[0];
sx q[0];
rz(-2.4910036) q[0];
sx q[0];
rz(-3.0537729) q[0];
x q[1];
rz(-1.6719867) q[2];
sx q[2];
rz(-1.8253583) q[2];
sx q[2];
rz(0.86092608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7192711) q[1];
sx q[1];
rz(-1.0485253) q[1];
sx q[1];
rz(2.2838255) q[1];
rz(-pi) q[2];
rz(-0.5011933) q[3];
sx q[3];
rz(-1.5287011) q[3];
sx q[3];
rz(0.30748366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6260234) q[2];
sx q[2];
rz(-1.7531351) q[2];
sx q[2];
rz(2.5173729) q[2];
rz(2.7624687) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(-2.8930801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17495951) q[0];
sx q[0];
rz(-2.3269854) q[0];
sx q[0];
rz(1.0354743) q[0];
rz(-1.8677208) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(-0.48286352) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5855756) q[0];
sx q[0];
rz(-0.10679467) q[0];
sx q[0];
rz(2.6983745) q[0];
rz(0.16629433) q[2];
sx q[2];
rz(-1.6876843) q[2];
sx q[2];
rz(2.6720195) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.00512) q[1];
sx q[1];
rz(-2.0025952) q[1];
sx q[1];
rz(3.1312607) q[1];
x q[2];
rz(-0.24077529) q[3];
sx q[3];
rz(-2.8117883) q[3];
sx q[3];
rz(0.93798104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1172993) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(-0.67000669) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(0.89699927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4033177) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(2.4025412) q[0];
rz(0.96145472) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(0.75698537) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805646) q[0];
sx q[0];
rz(-1.8317199) q[0];
sx q[0];
rz(-0.83103128) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4503373) q[2];
sx q[2];
rz(-1.2757841) q[2];
sx q[2];
rz(1.6694836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5685421) q[1];
sx q[1];
rz(-1.4366163) q[1];
sx q[1];
rz(2.5106984) q[1];
rz(1.2631031) q[3];
sx q[3];
rz(-1.6374) q[3];
sx q[3];
rz(2.8860725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5028533) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(-2.4364831) q[2];
rz(-2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(1.9158844) q[0];
rz(-1.907584) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(3.0753678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787251) q[0];
sx q[0];
rz(-1.6604794) q[0];
sx q[0];
rz(1.8171726) q[0];
x q[1];
rz(0.80133665) q[2];
sx q[2];
rz(-1.3532012) q[2];
sx q[2];
rz(-0.8720397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45715047) q[1];
sx q[1];
rz(-1.252414) q[1];
sx q[1];
rz(2.0778947) q[1];
rz(-pi) q[2];
rz(-2.8230571) q[3];
sx q[3];
rz(-0.38357601) q[3];
sx q[3];
rz(2.2936402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-0.44060102) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(-2.3073176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202268) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(-1.45586) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5680926) q[0];
sx q[0];
rz(-1.4254036) q[0];
sx q[0];
rz(2.306434) q[0];
x q[1];
rz(-0.46874502) q[2];
sx q[2];
rz(-2.8879117) q[2];
sx q[2];
rz(1.6799334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33259049) q[1];
sx q[1];
rz(-2.2223516) q[1];
sx q[1];
rz(1.5452191) q[1];
x q[2];
rz(-0.60815706) q[3];
sx q[3];
rz(-1.4440315) q[3];
sx q[3];
rz(2.7931282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8684034) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(2.2301162) q[2];
rz(0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-3.1349283) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8580496) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(-0.48769543) q[1];
sx q[1];
rz(-0.54134381) q[1];
sx q[1];
rz(-2.3977051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(1.6532142) q[0];
rz(-pi) q[1];
rz(3.0024372) q[2];
sx q[2];
rz(-2.1916813) q[2];
sx q[2];
rz(2.5992952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53092693) q[1];
sx q[1];
rz(-1.0333038) q[1];
sx q[1];
rz(-0.33157886) q[1];
x q[2];
rz(2.8091431) q[3];
sx q[3];
rz(-0.5088734) q[3];
sx q[3];
rz(2.8116941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0222212) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(-2.1807097) q[2];
rz(2.7458701) q[3];
sx q[3];
rz(-1.3362249) q[3];
sx q[3];
rz(3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845067) q[0];
sx q[0];
rz(-2.0261903) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(2.537435) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7121838) q[0];
sx q[0];
rz(-1.3954221) q[0];
sx q[0];
rz(-0.87923572) q[0];
rz(-pi) q[1];
rz(3.0839594) q[2];
sx q[2];
rz(-0.95256348) q[2];
sx q[2];
rz(1.0317486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8128374) q[1];
sx q[1];
rz(-0.27662524) q[1];
sx q[1];
rz(-3.0429716) q[1];
rz(-pi) q[2];
rz(-1.8363399) q[3];
sx q[3];
rz(-1.8887541) q[3];
sx q[3];
rz(-1.4972655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0543694) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(1.3537815) q[3];
sx q[3];
rz(-1.9446707) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488309) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(0.30717474) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-2.2408392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9884315) q[0];
sx q[0];
rz(-1.5876731) q[0];
sx q[0];
rz(-0.26217006) q[0];
rz(-pi) q[1];
rz(2.6866954) q[2];
sx q[2];
rz(-0.96856801) q[2];
sx q[2];
rz(1.0554316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14054243) q[1];
sx q[1];
rz(-1.8866354) q[1];
sx q[1];
rz(1.9204101) q[1];
rz(-1.4291271) q[3];
sx q[3];
rz(-2.1974539) q[3];
sx q[3];
rz(0.10654813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33425346) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(1.0963415) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39356247) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(2.6322741) q[0];
rz(-0.18109426) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(-2.1720355) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3345393) q[0];
sx q[0];
rz(-1.4167804) q[0];
sx q[0];
rz(-2.3898983) q[0];
rz(-pi) q[1];
rz(2.0779586) q[2];
sx q[2];
rz(-1.615287) q[2];
sx q[2];
rz(-2.8756093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0721832) q[1];
sx q[1];
rz(-1.2320215) q[1];
sx q[1];
rz(-1.6852001) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2990555) q[3];
sx q[3];
rz(-1.304783) q[3];
sx q[3];
rz(2.0913948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8934882) q[2];
sx q[2];
rz(-2.9073145) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.76081) q[0];
sx q[0];
rz(-0.2247227) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(2.6864247) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(-0.761935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2357016) q[0];
sx q[0];
rz(-2.3755223) q[0];
sx q[0];
rz(2.7188042) q[0];
rz(-pi) q[1];
rz(2.5013431) q[2];
sx q[2];
rz(-2.2138322) q[2];
sx q[2];
rz(-1.4852448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1208204) q[1];
sx q[1];
rz(-1.6050395) q[1];
sx q[1];
rz(-2.4832151) q[1];
rz(-pi) q[2];
rz(1.1033789) q[3];
sx q[3];
rz(-1.5776792) q[3];
sx q[3];
rz(0.21871834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(0.90395149) q[2];
rz(-1.5229185) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.9797549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.224613) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(2.0595111) q[2];
sx q[2];
rz(-0.50177028) q[2];
sx q[2];
rz(-0.72015464) q[2];
rz(0.0062777304) q[3];
sx q[3];
rz(-0.43000392) q[3];
sx q[3];
rz(3.1109711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
