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
rz(1.7215913) q[0];
sx q[0];
rz(-2.1821332) q[0];
sx q[0];
rz(1.0406915) q[0];
rz(2.8031082) q[1];
sx q[1];
rz(-1.6522633) q[1];
sx q[1];
rz(-1.7117865) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8393134) q[0];
sx q[0];
rz(-0.29713085) q[0];
sx q[0];
rz(-0.41584797) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0770887) q[2];
sx q[2];
rz(-2.4790451) q[2];
sx q[2];
rz(2.8665989) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46619848) q[1];
sx q[1];
rz(-2.3643092) q[1];
sx q[1];
rz(0.40730469) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8398714) q[3];
sx q[3];
rz(-2.1198049) q[3];
sx q[3];
rz(-1.3807856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0040943) q[2];
sx q[2];
rz(-1.7137004) q[2];
sx q[2];
rz(0.049169866) q[2];
rz(-2.5947425) q[3];
sx q[3];
rz(-2.2873736) q[3];
sx q[3];
rz(1.2949519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26286724) q[0];
sx q[0];
rz(-0.85064369) q[0];
sx q[0];
rz(-3.1021297) q[0];
rz(-0.70611686) q[1];
sx q[1];
rz(-1.1742274) q[1];
sx q[1];
rz(1.2605234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.880917) q[0];
sx q[0];
rz(-1.5124784) q[0];
sx q[0];
rz(-1.9888982) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4917212) q[2];
sx q[2];
rz(-0.77929493) q[2];
sx q[2];
rz(0.036606006) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83585868) q[1];
sx q[1];
rz(-1.6129588) q[1];
sx q[1];
rz(-3.107454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6373487) q[3];
sx q[3];
rz(-1.7643098) q[3];
sx q[3];
rz(-2.6849417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4827106) q[2];
sx q[2];
rz(-2.1363246) q[2];
sx q[2];
rz(0.91747326) q[2];
rz(-0.42258036) q[3];
sx q[3];
rz(-1.8924507) q[3];
sx q[3];
rz(1.0379855) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4611886) q[0];
sx q[0];
rz(-0.55484158) q[0];
sx q[0];
rz(0.94938266) q[0];
rz(1.2481015) q[1];
sx q[1];
rz(-0.58369842) q[1];
sx q[1];
rz(1.6277574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064977) q[0];
sx q[0];
rz(-1.2076006) q[0];
sx q[0];
rz(1.3083544) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1749985) q[2];
sx q[2];
rz(-0.51766073) q[2];
sx q[2];
rz(0.84830059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8366295) q[1];
sx q[1];
rz(-0.29759544) q[1];
sx q[1];
rz(-2.8862428) q[1];
rz(-pi) q[2];
rz(1.2825427) q[3];
sx q[3];
rz(-1.6950399) q[3];
sx q[3];
rz(-1.2584723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44024399) q[2];
sx q[2];
rz(-2.0155719) q[2];
sx q[2];
rz(1.6181642) q[2];
rz(-2.8273888) q[3];
sx q[3];
rz(-1.0256297) q[3];
sx q[3];
rz(1.5020802) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358129) q[0];
sx q[0];
rz(-2.950225) q[0];
sx q[0];
rz(-2.9845003) q[0];
rz(-3.0063903) q[1];
sx q[1];
rz(-2.356485) q[1];
sx q[1];
rz(-2.8083727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5557809) q[0];
sx q[0];
rz(-2.4436722) q[0];
sx q[0];
rz(0.74549739) q[0];
x q[1];
rz(-0.84246796) q[2];
sx q[2];
rz(-1.4204475) q[2];
sx q[2];
rz(-1.5068693) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1318479) q[1];
sx q[1];
rz(-1.8341843) q[1];
sx q[1];
rz(1.9326278) q[1];
rz(-pi) q[2];
rz(0.85835056) q[3];
sx q[3];
rz(-1.8110523) q[3];
sx q[3];
rz(-0.95360707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3335483) q[2];
sx q[2];
rz(-0.6490038) q[2];
sx q[2];
rz(1.6928847) q[2];
rz(0.73271218) q[3];
sx q[3];
rz(-1.2196187) q[3];
sx q[3];
rz(-1.609751) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919959) q[0];
sx q[0];
rz(-1.6549598) q[0];
sx q[0];
rz(2.9100371) q[0];
rz(-0.72136503) q[1];
sx q[1];
rz(-0.8911348) q[1];
sx q[1];
rz(0.84842938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285632) q[0];
sx q[0];
rz(-0.60270488) q[0];
sx q[0];
rz(-2.187285) q[0];
rz(-pi) q[1];
rz(-1.1662899) q[2];
sx q[2];
rz(-1.0644039) q[2];
sx q[2];
rz(0.17932349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2347772) q[1];
sx q[1];
rz(-0.39420745) q[1];
sx q[1];
rz(-0.18313198) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7293398) q[3];
sx q[3];
rz(-2.8839211) q[3];
sx q[3];
rz(0.6850971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1344177) q[2];
sx q[2];
rz(-0.60322064) q[2];
sx q[2];
rz(-1.2009386) q[2];
rz(-2.4522771) q[3];
sx q[3];
rz(-1.1700234) q[3];
sx q[3];
rz(-2.5758666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597964) q[0];
sx q[0];
rz(-1.5317651) q[0];
sx q[0];
rz(1.7933886) q[0];
rz(-0.7181522) q[1];
sx q[1];
rz(-0.57128692) q[1];
sx q[1];
rz(-1.2917554) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51707089) q[0];
sx q[0];
rz(-2.7291131) q[0];
sx q[0];
rz(-0.58858354) q[0];
rz(-pi) q[1];
rz(2.767105) q[2];
sx q[2];
rz(-1.0051703) q[2];
sx q[2];
rz(-2.781812) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2709658) q[1];
sx q[1];
rz(-0.18947225) q[1];
sx q[1];
rz(0.47196526) q[1];
rz(0.060543493) q[3];
sx q[3];
rz(-0.39602867) q[3];
sx q[3];
rz(2.372449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9888088) q[2];
sx q[2];
rz(-0.96502105) q[2];
sx q[2];
rz(-2.7394845) q[2];
rz(1.0058588) q[3];
sx q[3];
rz(-0.22336762) q[3];
sx q[3];
rz(-1.4373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724029) q[0];
sx q[0];
rz(-2.993354) q[0];
sx q[0];
rz(-1.8120026) q[0];
rz(-1.8136464) q[1];
sx q[1];
rz(-1.4168394) q[1];
sx q[1];
rz(-0.45164576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62147776) q[0];
sx q[0];
rz(-1.4796487) q[0];
sx q[0];
rz(-2.6976748) q[0];
rz(-pi) q[1];
rz(-0.30173413) q[2];
sx q[2];
rz(-0.44281755) q[2];
sx q[2];
rz(0.91609611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.627915) q[1];
sx q[1];
rz(-1.6213525) q[1];
sx q[1];
rz(-1.7611658) q[1];
rz(-pi) q[2];
rz(-2.3450646) q[3];
sx q[3];
rz(-2.3083271) q[3];
sx q[3];
rz(-3.0298373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20246501) q[2];
sx q[2];
rz(-1.3718601) q[2];
sx q[2];
rz(-2.8537214) q[2];
rz(1.6599844) q[3];
sx q[3];
rz(-2.2899254) q[3];
sx q[3];
rz(-1.456858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.562029) q[0];
sx q[0];
rz(-2.1810739) q[0];
sx q[0];
rz(0.55291837) q[0];
rz(1.2506073) q[1];
sx q[1];
rz(-0.41025531) q[1];
sx q[1];
rz(1.7611354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.834873) q[0];
sx q[0];
rz(-2.0836897) q[0];
sx q[0];
rz(1.1514649) q[0];
x q[1];
rz(0.12166656) q[2];
sx q[2];
rz(-0.98699283) q[2];
sx q[2];
rz(2.9505299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1637818) q[1];
sx q[1];
rz(-1.3339143) q[1];
sx q[1];
rz(1.6131748) q[1];
rz(0.80069009) q[3];
sx q[3];
rz(-1.691406) q[3];
sx q[3];
rz(-2.1757954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55014253) q[2];
sx q[2];
rz(-1.7937135) q[2];
sx q[2];
rz(2.4264917) q[2];
rz(3.0969369) q[3];
sx q[3];
rz(-0.5831334) q[3];
sx q[3];
rz(-2.513212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1735246) q[0];
sx q[0];
rz(-2.3723497) q[0];
sx q[0];
rz(-0.64088696) q[0];
rz(-1.6454654) q[1];
sx q[1];
rz(-0.38495266) q[1];
sx q[1];
rz(2.4900751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2007569) q[0];
sx q[0];
rz(-0.88229942) q[0];
sx q[0];
rz(-2.8066638) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8046494) q[2];
sx q[2];
rz(-2.049198) q[2];
sx q[2];
rz(-3.0870147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87043283) q[1];
sx q[1];
rz(-1.9385797) q[1];
sx q[1];
rz(-0.8720542) q[1];
rz(-pi) q[2];
rz(1.2226168) q[3];
sx q[3];
rz(-1.3908252) q[3];
sx q[3];
rz(-2.6258385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4064101) q[2];
sx q[2];
rz(-1.0446905) q[2];
sx q[2];
rz(-0.45808211) q[2];
rz(-3.1297019) q[3];
sx q[3];
rz(-1.0606822) q[3];
sx q[3];
rz(-0.021765821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7387725) q[0];
sx q[0];
rz(-0.194855) q[0];
sx q[0];
rz(3.1126157) q[0];
rz(0.92652357) q[1];
sx q[1];
rz(-1.69311) q[1];
sx q[1];
rz(0.40625939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0615912) q[0];
sx q[0];
rz(-1.5793945) q[0];
sx q[0];
rz(-2.8614869) q[0];
rz(-1.1154956) q[2];
sx q[2];
rz(-2.4433854) q[2];
sx q[2];
rz(-0.9309665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0613289) q[1];
sx q[1];
rz(-2.8483263) q[1];
sx q[1];
rz(1.6786511) q[1];
rz(-0.79979898) q[3];
sx q[3];
rz(-0.97808981) q[3];
sx q[3];
rz(-2.0519902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3709092) q[2];
sx q[2];
rz(-1.6434881) q[2];
sx q[2];
rz(-2.5468199) q[2];
rz(-0.70991436) q[3];
sx q[3];
rz(-2.1845332) q[3];
sx q[3];
rz(2.1962568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5668673) q[0];
sx q[0];
rz(-2.18676) q[0];
sx q[0];
rz(-2.2796897) q[0];
rz(-2.8361539) q[1];
sx q[1];
rz(-1.8157235) q[1];
sx q[1];
rz(-2.9095412) q[1];
rz(2.7965056) q[2];
sx q[2];
rz(-1.9660334) q[2];
sx q[2];
rz(2.6031969) q[2];
rz(-1.7954682) q[3];
sx q[3];
rz(-1.9763038) q[3];
sx q[3];
rz(1.4122813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
