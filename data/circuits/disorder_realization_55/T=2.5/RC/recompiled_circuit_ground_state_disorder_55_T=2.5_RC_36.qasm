OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3352873) q[0];
sx q[0];
rz(-1.0399613) q[0];
sx q[0];
rz(-0.34997532) q[0];
rz(1.5179874) q[1];
sx q[1];
rz(-2.1972158) q[1];
sx q[1];
rz(1.1854393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957592) q[0];
sx q[0];
rz(-1.1544219) q[0];
sx q[0];
rz(1.9365947) q[0];
x q[1];
rz(0.41857206) q[2];
sx q[2];
rz(-1.3579733) q[2];
sx q[2];
rz(-0.71453071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6373316) q[1];
sx q[1];
rz(-2.0954014) q[1];
sx q[1];
rz(2.8693512) q[1];
rz(0.18326944) q[3];
sx q[3];
rz(-2.1784867) q[3];
sx q[3];
rz(2.303108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27652937) q[2];
sx q[2];
rz(-2.7294366) q[2];
sx q[2];
rz(0.40974799) q[2];
rz(-0.77541882) q[3];
sx q[3];
rz(-1.5284458) q[3];
sx q[3];
rz(-2.8342136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85768098) q[0];
sx q[0];
rz(-2.6429521) q[0];
sx q[0];
rz(0.45251244) q[0];
rz(0.2287989) q[1];
sx q[1];
rz(-0.50140536) q[1];
sx q[1];
rz(2.9948044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75232601) q[0];
sx q[0];
rz(-1.379836) q[0];
sx q[0];
rz(1.2021078) q[0];
rz(-2.7640259) q[2];
sx q[2];
rz(-1.7930328) q[2];
sx q[2];
rz(1.3106032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2816946) q[1];
sx q[1];
rz(-1.7051776) q[1];
sx q[1];
rz(-1.5308389) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5654752) q[3];
sx q[3];
rz(-0.74339044) q[3];
sx q[3];
rz(-2.2593123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9734834) q[2];
sx q[2];
rz(-0.3349458) q[2];
sx q[2];
rz(2.2701263) q[2];
rz(2.787309) q[3];
sx q[3];
rz(-0.93577093) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5062434) q[0];
sx q[0];
rz(-0.62017089) q[0];
sx q[0];
rz(2.035602) q[0];
rz(-2.1982819) q[1];
sx q[1];
rz(-2.3425075) q[1];
sx q[1];
rz(1.4897289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2882799) q[0];
sx q[0];
rz(-1.5810471) q[0];
sx q[0];
rz(1.6985189) q[0];
rz(-pi) q[1];
rz(-0.42647894) q[2];
sx q[2];
rz(-0.22828776) q[2];
sx q[2];
rz(-0.84289614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4909648) q[1];
sx q[1];
rz(-2.5450383) q[1];
sx q[1];
rz(-3.0741398) q[1];
x q[2];
rz(-0.23814985) q[3];
sx q[3];
rz(-1.643073) q[3];
sx q[3];
rz(-1.0840462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51365435) q[2];
sx q[2];
rz(-1.1955806) q[2];
sx q[2];
rz(0.99620831) q[2];
rz(0.57322383) q[3];
sx q[3];
rz(-0.51155353) q[3];
sx q[3];
rz(-0.63428026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298249) q[0];
sx q[0];
rz(-1.8695762) q[0];
sx q[0];
rz(-1.5856532) q[0];
rz(-0.32444435) q[1];
sx q[1];
rz(-2.3093846) q[1];
sx q[1];
rz(-1.9437887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2749589) q[0];
sx q[0];
rz(-2.500545) q[0];
sx q[0];
rz(2.7085371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3487071) q[2];
sx q[2];
rz(-0.59704268) q[2];
sx q[2];
rz(2.3131903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19732319) q[1];
sx q[1];
rz(-1.0900647) q[1];
sx q[1];
rz(-0.14667796) q[1];
rz(-pi) q[2];
rz(2.1230202) q[3];
sx q[3];
rz(-2.8365144) q[3];
sx q[3];
rz(0.48459241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.045362554) q[2];
sx q[2];
rz(-2.5924293) q[2];
sx q[2];
rz(1.9722923) q[2];
rz(0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(-0.72871488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.8600334) q[0];
sx q[0];
rz(-0.39880729) q[0];
sx q[0];
rz(-1.1997724) q[0];
rz(-1.3072321) q[1];
sx q[1];
rz(-2.2016826) q[1];
sx q[1];
rz(-1.4365546) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2970151) q[0];
sx q[0];
rz(-2.7488193) q[0];
sx q[0];
rz(-1.7295087) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3398814) q[2];
sx q[2];
rz(-1.4428939) q[2];
sx q[2];
rz(-2.2890039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.349166) q[1];
sx q[1];
rz(-2.6712121) q[1];
sx q[1];
rz(-0.70655148) q[1];
rz(-pi) q[2];
rz(-0.9229696) q[3];
sx q[3];
rz(-1.9119091) q[3];
sx q[3];
rz(0.34527147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3522219) q[2];
sx q[2];
rz(-1.4381961) q[2];
sx q[2];
rz(-2.5686725) q[2];
rz(-2.2574183) q[3];
sx q[3];
rz(-0.97777647) q[3];
sx q[3];
rz(2.6463215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5837412) q[0];
sx q[0];
rz(-1.8889677) q[0];
sx q[0];
rz(1.0766693) q[0];
rz(2.63511) q[1];
sx q[1];
rz(-2.5264637) q[1];
sx q[1];
rz(-1.3004251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87430743) q[0];
sx q[0];
rz(-1.6194664) q[0];
sx q[0];
rz(-3.1278473) q[0];
x q[1];
rz(0.22407786) q[2];
sx q[2];
rz(-0.59106088) q[2];
sx q[2];
rz(0.056211744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36339736) q[1];
sx q[1];
rz(-1.4397238) q[1];
sx q[1];
rz(1.3568727) q[1];
rz(-0.20008101) q[3];
sx q[3];
rz(-1.225138) q[3];
sx q[3];
rz(2.8576771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41636813) q[2];
sx q[2];
rz(-0.84364265) q[2];
sx q[2];
rz(-1.9492524) q[2];
rz(-0.013817712) q[3];
sx q[3];
rz(-0.71607029) q[3];
sx q[3];
rz(0.97318399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8927807) q[0];
sx q[0];
rz(-0.92614323) q[0];
sx q[0];
rz(-2.7698351) q[0];
rz(-0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(0.54381347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1756729) q[0];
sx q[0];
rz(-1.9750711) q[0];
sx q[0];
rz(2.5258397) q[0];
rz(2.6863348) q[2];
sx q[2];
rz(-0.90418679) q[2];
sx q[2];
rz(-2.8316759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84219599) q[1];
sx q[1];
rz(-2.7666515) q[1];
sx q[1];
rz(0.018126651) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4187846) q[3];
sx q[3];
rz(-1.6789888) q[3];
sx q[3];
rz(-0.044807981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1112993) q[2];
sx q[2];
rz(-1.0617504) q[2];
sx q[2];
rz(0.57478762) q[2];
rz(0.8542257) q[3];
sx q[3];
rz(-1.6653019) q[3];
sx q[3];
rz(2.0265719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92693555) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(-0.32447234) q[0];
rz(-2.532161) q[1];
sx q[1];
rz(-2.29988) q[1];
sx q[1];
rz(-2.6268974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97159779) q[0];
sx q[0];
rz(-2.7063402) q[0];
sx q[0];
rz(1.5488141) q[0];
rz(-pi) q[1];
rz(-1.1616522) q[2];
sx q[2];
rz(-2.3858641) q[2];
sx q[2];
rz(0.31015304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60694164) q[1];
sx q[1];
rz(-1.3365627) q[1];
sx q[1];
rz(-2.5638084) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2949766) q[3];
sx q[3];
rz(-1.6758741) q[3];
sx q[3];
rz(0.17227473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9784166) q[2];
sx q[2];
rz(-1.6360444) q[2];
sx q[2];
rz(2.8669538) q[2];
rz(2.2091673) q[3];
sx q[3];
rz(-0.43007389) q[3];
sx q[3];
rz(0.29621223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68359971) q[0];
sx q[0];
rz(-1.4608811) q[0];
sx q[0];
rz(0.13548166) q[0];
rz(-0.97700351) q[1];
sx q[1];
rz(-2.3833387) q[1];
sx q[1];
rz(3.1405084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8086804) q[0];
sx q[0];
rz(-1.7263733) q[0];
sx q[0];
rz(-2.6434487) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7747701) q[2];
sx q[2];
rz(-2.6088031) q[2];
sx q[2];
rz(0.29292695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3848968) q[1];
sx q[1];
rz(-2.2731509) q[1];
sx q[1];
rz(0.76354247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3850609) q[3];
sx q[3];
rz(-1.0660785) q[3];
sx q[3];
rz(1.0738508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5929426) q[2];
sx q[2];
rz(-1.0385916) q[2];
sx q[2];
rz(2.9366117) q[2];
rz(-0.18260469) q[3];
sx q[3];
rz(-0.20668106) q[3];
sx q[3];
rz(-0.73777795) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36271998) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(-2.6642098) q[0];
rz(0.1167156) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(0.83888549) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6434833) q[0];
sx q[0];
rz(-2.2543397) q[0];
sx q[0];
rz(-0.34434005) q[0];
rz(-1.2195385) q[2];
sx q[2];
rz(-2.2645586) q[2];
sx q[2];
rz(1.1379682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9899163) q[1];
sx q[1];
rz(-2.4258759) q[1];
sx q[1];
rz(2.8552961) q[1];
rz(-1.4492744) q[3];
sx q[3];
rz(-0.24059939) q[3];
sx q[3];
rz(-0.46328637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3597151) q[2];
sx q[2];
rz(-0.30656591) q[2];
sx q[2];
rz(0.26620418) q[2];
rz(-0.5075469) q[3];
sx q[3];
rz(-0.72269732) q[3];
sx q[3];
rz(-0.44107309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97892852) q[0];
sx q[0];
rz(-1.6829818) q[0];
sx q[0];
rz(-0.80243954) q[0];
rz(-2.2294527) q[1];
sx q[1];
rz(-1.9504539) q[1];
sx q[1];
rz(0.84074195) q[1];
rz(1.4100762) q[2];
sx q[2];
rz(-1.3257324) q[2];
sx q[2];
rz(-0.64468862) q[2];
rz(2.38337) q[3];
sx q[3];
rz(-2.6144386) q[3];
sx q[3];
rz(2.2948882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
