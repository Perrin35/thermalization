OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.469874) q[0];
sx q[0];
rz(4.2966312) q[0];
sx q[0];
rz(11.136461) q[0];
rz(-2.6861796) q[1];
sx q[1];
rz(-2.5502584) q[1];
sx q[1];
rz(1.1270181) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15412946) q[0];
sx q[0];
rz(-2.9107981) q[0];
sx q[0];
rz(1.5202281) q[0];
rz(-pi) q[1];
rz(-0.34268171) q[2];
sx q[2];
rz(-1.3200054) q[2];
sx q[2];
rz(-0.33736526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0309567) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(2.8085676) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5529049) q[3];
sx q[3];
rz(-1.3800637) q[3];
sx q[3];
rz(2.408556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0226125) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(1.653778) q[2];
rz(-2.7704499) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(0.77905542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30598518) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(0.66407472) q[0];
rz(0.56655073) q[1];
sx q[1];
rz(-1.605875) q[1];
sx q[1];
rz(-2.9552592) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16435745) q[0];
sx q[0];
rz(-2.1002227) q[0];
sx q[0];
rz(0.25821547) q[0];
rz(1.1900224) q[2];
sx q[2];
rz(-1.3212034) q[2];
sx q[2];
rz(-1.6935503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9606321) q[1];
sx q[1];
rz(-2.4419028) q[1];
sx q[1];
rz(-1.8262499) q[1];
rz(-1.4545031) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8191007) q[2];
sx q[2];
rz(-1.2995316) q[2];
sx q[2];
rz(-1.5619649) q[2];
rz(1.1473848) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(-1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086562432) q[0];
sx q[0];
rz(-1.6494305) q[0];
sx q[0];
rz(-2.8360039) q[0];
rz(-1.2809523) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(1.618128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9394835) q[0];
sx q[0];
rz(-1.8022707) q[0];
sx q[0];
rz(-0.37136308) q[0];
rz(-pi) q[1];
rz(-0.32714897) q[2];
sx q[2];
rz(-0.047657813) q[2];
sx q[2];
rz(2.7313155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.87464) q[1];
sx q[1];
rz(-2.2494509) q[1];
sx q[1];
rz(2.0372422) q[1];
rz(2.9039246) q[3];
sx q[3];
rz(-2.1690627) q[3];
sx q[3];
rz(-3.0868343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4898701) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(0.22001246) q[2];
rz(3.0618727) q[3];
sx q[3];
rz(-2.1646175) q[3];
sx q[3];
rz(-0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43059573) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(0.0084477607) q[0];
rz(-1.1620109) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(1.1147503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24766391) q[0];
sx q[0];
rz(-1.8801483) q[0];
sx q[0];
rz(0.19750144) q[0];
rz(-pi) q[1];
rz(-2.5774245) q[2];
sx q[2];
rz(-2.0551066) q[2];
sx q[2];
rz(-0.51646215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47059688) q[1];
sx q[1];
rz(-0.45491114) q[1];
sx q[1];
rz(0.90683337) q[1];
x q[2];
rz(-0.77372257) q[3];
sx q[3];
rz(-2.5829909) q[3];
sx q[3];
rz(-0.62880558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4343425) q[2];
sx q[2];
rz(-2.12314) q[2];
sx q[2];
rz(-0.84558359) q[2];
rz(-0.25104684) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(-0.74560753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276328) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(2.1115671) q[0];
rz(2.5469942) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(1.3535708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6973482) q[0];
sx q[0];
rz(-1.5500907) q[0];
sx q[0];
rz(0.022788825) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3859904) q[2];
sx q[2];
rz(-1.2535742) q[2];
sx q[2];
rz(0.76390195) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.374267) q[1];
sx q[1];
rz(-2.0637636) q[1];
sx q[1];
rz(2.3460991) q[1];
x q[2];
rz(2.0849801) q[3];
sx q[3];
rz(-1.0407036) q[3];
sx q[3];
rz(-0.3584273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9109965) q[2];
sx q[2];
rz(-1.6350919) q[2];
sx q[2];
rz(-0.050749151) q[2];
rz(1.1132318) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-2.7509403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0473061) q[0];
sx q[0];
rz(-1.0599437) q[0];
sx q[0];
rz(0.37288368) q[0];
rz(1.8376384) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(2.1162927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4639125) q[0];
sx q[0];
rz(-0.45651528) q[0];
sx q[0];
rz(0.17386659) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0228902) q[2];
sx q[2];
rz(-2.7774924) q[2];
sx q[2];
rz(-0.61227476) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48580468) q[1];
sx q[1];
rz(-0.88690573) q[1];
sx q[1];
rz(-2.8832316) q[1];
rz(-pi) q[2];
rz(1.147109) q[3];
sx q[3];
rz(-1.7762587) q[3];
sx q[3];
rz(0.68483554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.23124) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.7115889) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-0.68683306) q[3];
sx q[3];
rz(1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1533399) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(2.7100995) q[0];
rz(-1.8776114) q[1];
sx q[1];
rz(-1.6004205) q[1];
sx q[1];
rz(2.4914609) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7509814) q[0];
sx q[0];
rz(-1.6529875) q[0];
sx q[0];
rz(2.85336) q[0];
rz(0.88754295) q[2];
sx q[2];
rz(-1.6517963) q[2];
sx q[2];
rz(2.0320867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0228867) q[1];
sx q[1];
rz(-2.431015) q[1];
sx q[1];
rz(-1.7479244) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5863083) q[3];
sx q[3];
rz(-2.6836694) q[3];
sx q[3];
rz(2.0934888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5905137) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(3.1371269) q[2];
rz(0.0095857754) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(-2.7981304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15653217) q[0];
sx q[0];
rz(-0.28230202) q[0];
sx q[0];
rz(0.16047934) q[0];
rz(1.7566682) q[1];
sx q[1];
rz(-0.79912186) q[1];
sx q[1];
rz(-2.8327732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35058103) q[0];
sx q[0];
rz(-1.6948943) q[0];
sx q[0];
rz(1.3826293) q[0];
rz(-pi) q[1];
rz(-0.71134348) q[2];
sx q[2];
rz(-1.2738196) q[2];
sx q[2];
rz(-2.1893196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2330197) q[1];
sx q[1];
rz(-1.9601788) q[1];
sx q[1];
rz(1.1638457) q[1];
rz(-pi) q[2];
rz(-1.2587955) q[3];
sx q[3];
rz(-1.9083169) q[3];
sx q[3];
rz(2.5520293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34362346) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(-1.5416175) q[2];
rz(0.68495098) q[3];
sx q[3];
rz(-1.5870321) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(0.55139971) q[0];
rz(-0.4772056) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(-2.3068857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4906368) q[0];
sx q[0];
rz(-1.0792562) q[0];
sx q[0];
rz(1.0440582) q[0];
x q[1];
rz(-1.536252) q[2];
sx q[2];
rz(-0.87022034) q[2];
sx q[2];
rz(0.81734428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0352201) q[1];
sx q[1];
rz(-2.1823931) q[1];
sx q[1];
rz(-0.04208438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3605462) q[3];
sx q[3];
rz(-2.3608942) q[3];
sx q[3];
rz(2.4311284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4386374) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(2.9023602) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(-1.5391866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5156373) q[0];
sx q[0];
rz(-3.0037168) q[0];
sx q[0];
rz(-2.4249518) q[0];
rz(0.58249885) q[1];
sx q[1];
rz(-1.7889675) q[1];
sx q[1];
rz(2.8315721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48381708) q[0];
sx q[0];
rz(-0.80379009) q[0];
sx q[0];
rz(-1.0763542) q[0];
x q[1];
rz(-1.9642682) q[2];
sx q[2];
rz(-1.9542017) q[2];
sx q[2];
rz(-2.368106) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88921949) q[1];
sx q[1];
rz(-2.3205726) q[1];
sx q[1];
rz(-2.0682425) q[1];
rz(-1.0729372) q[3];
sx q[3];
rz(-1.080435) q[3];
sx q[3];
rz(2.2435202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(1.4198111) q[2];
rz(1.4451197) q[3];
sx q[3];
rz(-2.4308379) q[3];
sx q[3];
rz(-2.3783309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2390908) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(1.6613962) q[1];
sx q[1];
rz(-1.7330685) q[1];
sx q[1];
rz(0.20191244) q[1];
rz(0.26915941) q[2];
sx q[2];
rz(-2.8648389) q[2];
sx q[2];
rz(1.7493389) q[2];
rz(-2.2524184) q[3];
sx q[3];
rz(-1.9445264) q[3];
sx q[3];
rz(1.0497937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
