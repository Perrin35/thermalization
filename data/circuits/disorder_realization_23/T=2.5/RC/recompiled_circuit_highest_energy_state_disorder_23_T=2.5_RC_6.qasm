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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(1.1516655) q[0];
rz(-2.0517853) q[1];
sx q[1];
rz(-1.6648219) q[1];
sx q[1];
rz(2.9906315) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.359084) q[0];
sx q[0];
rz(-1.5182966) q[0];
sx q[0];
rz(1.2760602) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8279575) q[2];
sx q[2];
rz(-1.6534963) q[2];
sx q[2];
rz(-1.5717779) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6358984) q[1];
sx q[1];
rz(-1.5855967) q[1];
sx q[1];
rz(3.1379051) q[1];
rz(2.742953) q[3];
sx q[3];
rz(-0.10484914) q[3];
sx q[3];
rz(2.9440126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0207409) q[2];
sx q[2];
rz(-3.1181702) q[2];
sx q[2];
rz(2.6502996) q[2];
rz(1.8190207) q[3];
sx q[3];
rz(-1.5343851) q[3];
sx q[3];
rz(-1.8137431) q[3];
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
rz(0.78111929) q[0];
sx q[0];
rz(-2.5871215) q[0];
sx q[0];
rz(2.4419899) q[0];
rz(-1.5505002) q[1];
sx q[1];
rz(-2.6467549) q[1];
sx q[1];
rz(-2.9717305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9014329) q[0];
sx q[0];
rz(-0.086726464) q[0];
sx q[0];
rz(-2.6804352) q[0];
rz(-pi) q[1];
rz(-1.4361977) q[2];
sx q[2];
rz(-1.2685809) q[2];
sx q[2];
rz(2.5562191) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8123229) q[1];
sx q[1];
rz(-1.0877177) q[1];
sx q[1];
rz(-0.4050576) q[1];
rz(-1.1137257) q[3];
sx q[3];
rz(-1.4997484) q[3];
sx q[3];
rz(0.20807264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26213172) q[2];
sx q[2];
rz(-0.36236557) q[2];
sx q[2];
rz(1.7246838) q[2];
rz(-2.0587685) q[3];
sx q[3];
rz(-0.88734752) q[3];
sx q[3];
rz(-0.82720238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5603492) q[0];
sx q[0];
rz(-0.92009783) q[0];
sx q[0];
rz(1.5592519) q[0];
rz(2.4371367) q[1];
sx q[1];
rz(-1.9943941) q[1];
sx q[1];
rz(0.53680435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671331) q[0];
sx q[0];
rz(-1.6170814) q[0];
sx q[0];
rz(0.68741902) q[0];
rz(-1.6073999) q[2];
sx q[2];
rz(-1.6609909) q[2];
sx q[2];
rz(-0.48675234) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0276913) q[1];
sx q[1];
rz(-1.4639106) q[1];
sx q[1];
rz(2.9902469) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2800673) q[3];
sx q[3];
rz(-2.42413) q[3];
sx q[3];
rz(-2.2253401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.74865666) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(0.7788457) q[2];
rz(0.084176453) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(-1.4028153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.2445143) q[0];
sx q[0];
rz(-0.92105138) q[0];
sx q[0];
rz(0.47065863) q[0];
rz(1.6828407) q[1];
sx q[1];
rz(-3.1367446) q[1];
sx q[1];
rz(2.2768314) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7047459) q[0];
sx q[0];
rz(-0.517874) q[0];
sx q[0];
rz(-3.0240866) q[0];
x q[1];
rz(-1.7030348) q[2];
sx q[2];
rz(-1.3515499) q[2];
sx q[2];
rz(-1.8869149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7219077) q[1];
sx q[1];
rz(-0.083217155) q[1];
sx q[1];
rz(-2.2305942) q[1];
rz(2.6506977) q[3];
sx q[3];
rz(-2.5449736) q[3];
sx q[3];
rz(0.20363775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1137769) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(-1.1569542) q[2];
rz(2.786934) q[3];
sx q[3];
rz(-1.9054474) q[3];
sx q[3];
rz(3.0310596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.81195867) q[0];
sx q[0];
rz(-0.93557731) q[0];
sx q[0];
rz(2.6079566) q[0];
rz(1.2136906) q[1];
sx q[1];
rz(-3.1184989) q[1];
sx q[1];
rz(-1.9603221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034485949) q[0];
sx q[0];
rz(-1.6167322) q[0];
sx q[0];
rz(2.9784543) q[0];
x q[1];
rz(-0.65106647) q[2];
sx q[2];
rz(-1.2322556) q[2];
sx q[2];
rz(-0.19442788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86210862) q[1];
sx q[1];
rz(-0.92593926) q[1];
sx q[1];
rz(2.9456861) q[1];
rz(-1.3245246) q[3];
sx q[3];
rz(-2.6239873) q[3];
sx q[3];
rz(-2.8235265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.033325) q[2];
sx q[2];
rz(-0.51563087) q[2];
sx q[2];
rz(0.92833129) q[2];
rz(2.8632274) q[3];
sx q[3];
rz(-1.3185578) q[3];
sx q[3];
rz(-3.0729455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5334897) q[0];
sx q[0];
rz(-2.6187496) q[0];
sx q[0];
rz(-2.9485517) q[0];
rz(-0.67735425) q[1];
sx q[1];
rz(-0.028379863) q[1];
sx q[1];
rz(-1.5100286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0832779) q[0];
sx q[0];
rz(-1.711635) q[0];
sx q[0];
rz(1.5910287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6608428) q[2];
sx q[2];
rz(-2.6777732) q[2];
sx q[2];
rz(-3.0141789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2360797) q[1];
sx q[1];
rz(-1.326099) q[1];
sx q[1];
rz(-0.83531816) q[1];
rz(-0.55431788) q[3];
sx q[3];
rz(-2.333967) q[3];
sx q[3];
rz(0.19632158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3626927) q[2];
sx q[2];
rz(-2.429246) q[2];
sx q[2];
rz(0.98711291) q[2];
rz(1.0846064) q[3];
sx q[3];
rz(-0.54607138) q[3];
sx q[3];
rz(0.62091056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5793594) q[0];
sx q[0];
rz(-0.99246228) q[0];
sx q[0];
rz(1.5480504) q[0];
rz(-0.69001895) q[1];
sx q[1];
rz(-3.1179805) q[1];
sx q[1];
rz(-1.7049559) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48955917) q[0];
sx q[0];
rz(-1.0512182) q[0];
sx q[0];
rz(1.5205617) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1146554) q[2];
sx q[2];
rz(-2.4566894) q[2];
sx q[2];
rz(1.2482582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7014536) q[1];
sx q[1];
rz(-1.5112644) q[1];
sx q[1];
rz(1.0284222) q[1];
rz(-2.6855647) q[3];
sx q[3];
rz(-1.8397962) q[3];
sx q[3];
rz(1.8331127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.37316608) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(-2.5567059) q[2];
rz(-1.5960426) q[3];
sx q[3];
rz(-1.7187748) q[3];
sx q[3];
rz(0.91442951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5698513) q[0];
sx q[0];
rz(-2.2241156) q[0];
sx q[0];
rz(-1.5745987) q[0];
rz(-0.51708108) q[1];
sx q[1];
rz(-0.93544975) q[1];
sx q[1];
rz(1.8749974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5735949) q[0];
sx q[0];
rz(-1.5758152) q[0];
sx q[0];
rz(-1.5673593) q[0];
x q[1];
rz(1.0158407) q[2];
sx q[2];
rz(-0.84677831) q[2];
sx q[2];
rz(-2.1479549) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8116906) q[1];
sx q[1];
rz(-1.2214097) q[1];
sx q[1];
rz(2.1495649) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.35584) q[3];
sx q[3];
rz(-2.9693797) q[3];
sx q[3];
rz(0.39628753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37303415) q[2];
sx q[2];
rz(-0.80058241) q[2];
sx q[2];
rz(-1.8689092) q[2];
rz(-0.4396762) q[3];
sx q[3];
rz(-1.9821854) q[3];
sx q[3];
rz(-0.064597733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.30176297) q[0];
sx q[0];
rz(-3.1253212) q[0];
sx q[0];
rz(-2.8453258) q[0];
rz(-3.0773194) q[1];
sx q[1];
rz(-1.4646894) q[1];
sx q[1];
rz(1.6305264) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8730174) q[0];
sx q[0];
rz(-0.35726705) q[0];
sx q[0];
rz(1.9505469) q[0];
x q[1];
rz(-0.6728386) q[2];
sx q[2];
rz(-1.4110067) q[2];
sx q[2];
rz(2.9304913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0546186) q[1];
sx q[1];
rz(-1.6999632) q[1];
sx q[1];
rz(1.8968326) q[1];
rz(0.011476403) q[3];
sx q[3];
rz(-1.6919129) q[3];
sx q[3];
rz(0.83601421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5955547) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(1.1240553) q[2];
rz(-1.3042287) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(-0.058163253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8703363) q[0];
sx q[0];
rz(-0.79126343) q[0];
sx q[0];
rz(1.1668209) q[0];
rz(1.600949) q[1];
sx q[1];
rz(-0.29778844) q[1];
sx q[1];
rz(-1.3265532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8542503) q[0];
sx q[0];
rz(-2.4357987) q[0];
sx q[0];
rz(0.38278929) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4197056) q[2];
sx q[2];
rz(-0.52209243) q[2];
sx q[2];
rz(2.4185857) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9611813) q[1];
sx q[1];
rz(-1.3588393) q[1];
sx q[1];
rz(2.4418264) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3438026) q[3];
sx q[3];
rz(-0.41068893) q[3];
sx q[3];
rz(-2.3410377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88663179) q[2];
sx q[2];
rz(-2.9780264) q[2];
sx q[2];
rz(1.6939885) q[2];
rz(0.12668954) q[3];
sx q[3];
rz(-1.4796175) q[3];
sx q[3];
rz(-1.1160342) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806725) q[0];
sx q[0];
rz(-1.3074449) q[0];
sx q[0];
rz(1.773651) q[0];
rz(1.548829) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(1.7573865) q[2];
sx q[2];
rz(-2.0108847) q[2];
sx q[2];
rz(-0.080367676) q[2];
rz(2.1550989) q[3];
sx q[3];
rz(-1.7135847) q[3];
sx q[3];
rz(1.8691487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
