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
rz(-0.65547216) q[0];
sx q[0];
rz(-2.1894426) q[0];
sx q[0];
rz(0.093753554) q[0];
rz(1.3682415) q[1];
sx q[1];
rz(4.7028766) q[1];
sx q[1];
rz(8.75755) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3190368) q[0];
sx q[0];
rz(-2.0635491) q[0];
sx q[0];
rz(-0.42903729) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7736499) q[2];
sx q[2];
rz(-1.8241992) q[2];
sx q[2];
rz(-2.5784166) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2591483) q[1];
sx q[1];
rz(-1.8516774) q[1];
sx q[1];
rz(0.85390635) q[1];
rz(2.2556858) q[3];
sx q[3];
rz(-2.0768171) q[3];
sx q[3];
rz(0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61624709) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(0.023539143) q[2];
rz(2.8871138) q[3];
sx q[3];
rz(-1.1602217) q[3];
sx q[3];
rz(0.75828534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94754058) q[0];
sx q[0];
rz(-1.3964615) q[0];
sx q[0];
rz(-0.61554712) q[0];
rz(2.5415892) q[1];
sx q[1];
rz(-1.871385) q[1];
sx q[1];
rz(-2.792865) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6541512) q[0];
sx q[0];
rz(-0.92986996) q[0];
sx q[0];
rz(1.4750255) q[0];
rz(-pi) q[1];
rz(-0.87375529) q[2];
sx q[2];
rz(-0.43461576) q[2];
sx q[2];
rz(-0.86948621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24838994) q[1];
sx q[1];
rz(-0.55596184) q[1];
sx q[1];
rz(0.979579) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28517752) q[3];
sx q[3];
rz(-1.0833519) q[3];
sx q[3];
rz(0.91966682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6500924) q[2];
sx q[2];
rz(-1.3050175) q[2];
sx q[2];
rz(-0.65518641) q[2];
rz(0.16264597) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(0.20774016) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1732037) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(0.094245687) q[0];
rz(-1.3000129) q[1];
sx q[1];
rz(-2.2424707) q[1];
sx q[1];
rz(1.1945486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.081095) q[0];
sx q[0];
rz(-2.0682242) q[0];
sx q[0];
rz(1.917385) q[0];
x q[1];
rz(-1.934951) q[2];
sx q[2];
rz(-1.6192163) q[2];
sx q[2];
rz(-0.044755699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.647741) q[1];
sx q[1];
rz(-1.5268699) q[1];
sx q[1];
rz(-1.8274587) q[1];
rz(0.76285513) q[3];
sx q[3];
rz(-1.8486946) q[3];
sx q[3];
rz(-2.7369473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3208348) q[2];
sx q[2];
rz(-1.4471549) q[2];
sx q[2];
rz(-2.7336332) q[2];
rz(1.7631433) q[3];
sx q[3];
rz(-0.89865509) q[3];
sx q[3];
rz(-0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5648062) q[0];
sx q[0];
rz(-1.3917568) q[0];
sx q[0];
rz(-0.80192649) q[0];
rz(0.39464125) q[1];
sx q[1];
rz(-2.092974) q[1];
sx q[1];
rz(-0.035042979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76709439) q[0];
sx q[0];
rz(-1.2001925) q[0];
sx q[0];
rz(0.9824533) q[0];
x q[1];
rz(0.50582992) q[2];
sx q[2];
rz(-2.7625045) q[2];
sx q[2];
rz(2.068678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9600507) q[1];
sx q[1];
rz(-1.6142577) q[1];
sx q[1];
rz(0.04219136) q[1];
rz(-pi) q[2];
rz(2.6529516) q[3];
sx q[3];
rz(-0.39791574) q[3];
sx q[3];
rz(-0.15109381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0617712) q[2];
sx q[2];
rz(-2.0786736) q[2];
sx q[2];
rz(-1.4351832) q[2];
rz(-1.7363413) q[3];
sx q[3];
rz(-2.0515714) q[3];
sx q[3];
rz(2.0315571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46471304) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(-1.9997464) q[0];
rz(-2.6308718) q[1];
sx q[1];
rz(-1.2341713) q[1];
sx q[1];
rz(1.4265192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1349484) q[0];
sx q[0];
rz(-1.0492968) q[0];
sx q[0];
rz(-2.5000076) q[0];
rz(-pi) q[1];
rz(0.22682206) q[2];
sx q[2];
rz(-1.2661627) q[2];
sx q[2];
rz(-0.90639988) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.63663) q[1];
sx q[1];
rz(-0.98082029) q[1];
sx q[1];
rz(-3.0061199) q[1];
x q[2];
rz(2.9742091) q[3];
sx q[3];
rz(-1.9731083) q[3];
sx q[3];
rz(-1.2237807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.3072689) q[2];
sx q[2];
rz(0.2571787) q[2];
rz(0.87279618) q[3];
sx q[3];
rz(-1.3215439) q[3];
sx q[3];
rz(2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207101) q[0];
sx q[0];
rz(-0.45419422) q[0];
sx q[0];
rz(2.0632451) q[0];
rz(-2.2348166) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(0.041898601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7920162) q[0];
sx q[0];
rz(-1.0984165) q[0];
sx q[0];
rz(-3.1012606) q[0];
rz(2.2380377) q[2];
sx q[2];
rz(-1.6674526) q[2];
sx q[2];
rz(-0.54723155) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4087569) q[1];
sx q[1];
rz(-2.6799767) q[1];
sx q[1];
rz(-2.2240586) q[1];
rz(-0.53421212) q[3];
sx q[3];
rz(-2.231866) q[3];
sx q[3];
rz(-2.7345757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8580253) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(-2.2502327) q[2];
rz(-0.71409613) q[3];
sx q[3];
rz(-0.73914206) q[3];
sx q[3];
rz(-1.063063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080268) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(-2.4666393) q[0];
rz(-1.630111) q[1];
sx q[1];
rz(-2.5210896) q[1];
sx q[1];
rz(-0.57704467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9723608) q[0];
sx q[0];
rz(-1.0857538) q[0];
sx q[0];
rz(0.26974704) q[0];
rz(-2.7084025) q[2];
sx q[2];
rz(-1.4775039) q[2];
sx q[2];
rz(-2.6733638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.1516421) q[1];
sx q[1];
rz(-0.90817876) q[1];
sx q[1];
rz(3.0376787) q[1];
rz(-pi) q[2];
rz(-2.2636599) q[3];
sx q[3];
rz(-2.0725277) q[3];
sx q[3];
rz(-0.89661613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.453489) q[2];
sx q[2];
rz(-2.0192396) q[2];
sx q[2];
rz(0.92448676) q[2];
rz(1.2201355) q[3];
sx q[3];
rz(-2.852738) q[3];
sx q[3];
rz(0.060062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7500551) q[0];
sx q[0];
rz(-0.67111641) q[0];
sx q[0];
rz(1.7975988) q[0];
rz(-1.2229819) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(-1.5917684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5179199) q[0];
sx q[0];
rz(-2.8123283) q[0];
sx q[0];
rz(0.98401208) q[0];
rz(-pi) q[1];
rz(-0.37495592) q[2];
sx q[2];
rz(-1.6757312) q[2];
sx q[2];
rz(-2.8255759) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38993714) q[1];
sx q[1];
rz(-0.90607535) q[1];
sx q[1];
rz(1.0996363) q[1];
rz(-3.0430085) q[3];
sx q[3];
rz(-2.3752593) q[3];
sx q[3];
rz(0.51829547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.082077114) q[2];
sx q[2];
rz(-1.0249219) q[2];
sx q[2];
rz(-0.22641851) q[2];
rz(-0.039479937) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(1.9421019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3969642) q[0];
sx q[0];
rz(-1.2122943) q[0];
sx q[0];
rz(2.3126171) q[0];
rz(0.12116155) q[1];
sx q[1];
rz(-0.96210259) q[1];
sx q[1];
rz(-1.4519579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71428107) q[0];
sx q[0];
rz(-1.7085008) q[0];
sx q[0];
rz(2.4192823) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.317886) q[2];
sx q[2];
rz(-1.2531157) q[2];
sx q[2];
rz(1.6967271) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.02579298) q[1];
sx q[1];
rz(-2.3069972) q[1];
sx q[1];
rz(-0.92206149) q[1];
rz(2.7424654) q[3];
sx q[3];
rz(-2.0889971) q[3];
sx q[3];
rz(-0.46405989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58763233) q[2];
sx q[2];
rz(-0.82531896) q[2];
sx q[2];
rz(0.86453214) q[2];
rz(-0.65166059) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(0.017875044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.637735) q[0];
sx q[0];
rz(-0.88541579) q[0];
sx q[0];
rz(-2.6813685) q[0];
rz(1.3453311) q[1];
sx q[1];
rz(-2.5049152) q[1];
sx q[1];
rz(-1.771079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84708285) q[0];
sx q[0];
rz(-2.6381603) q[0];
sx q[0];
rz(2.2398021) q[0];
rz(-1.2865169) q[2];
sx q[2];
rz(-2.1889969) q[2];
sx q[2];
rz(2.6524891) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7975841) q[1];
sx q[1];
rz(-1.3329778) q[1];
sx q[1];
rz(1.8348376) q[1];
x q[2];
rz(1.7779782) q[3];
sx q[3];
rz(-1.6231114) q[3];
sx q[3];
rz(-0.68488831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7728077) q[2];
sx q[2];
rz(-1.7926755) q[2];
sx q[2];
rz(0.14300145) q[2];
rz(-2.9755106) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(-2.3486229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76846692) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(-1.985818) q[1];
sx q[1];
rz(-1.118569) q[1];
sx q[1];
rz(-2.7973693) q[1];
rz(-2.106582) q[2];
sx q[2];
rz(-2.2314318) q[2];
sx q[2];
rz(-0.4690276) q[2];
rz(-1.0977911) q[3];
sx q[3];
rz(-2.4840631) q[3];
sx q[3];
rz(-2.4854131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
