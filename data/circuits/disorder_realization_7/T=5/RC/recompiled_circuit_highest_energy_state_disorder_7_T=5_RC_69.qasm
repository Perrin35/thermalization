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
rz(0.12662521) q[0];
sx q[0];
rz(1.5746483) q[0];
sx q[0];
rz(9.9404542) q[0];
rz(-2.9134143) q[1];
sx q[1];
rz(-2.394634) q[1];
sx q[1];
rz(2.7153314) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92862064) q[0];
sx q[0];
rz(-2.3317695) q[0];
sx q[0];
rz(-2.9024603) q[0];
x q[1];
rz(-0.49246712) q[2];
sx q[2];
rz(-0.55869192) q[2];
sx q[2];
rz(-0.44670263) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8350153) q[1];
sx q[1];
rz(-1.6284429) q[1];
sx q[1];
rz(-0.60204864) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98661749) q[3];
sx q[3];
rz(-1.8446577) q[3];
sx q[3];
rz(-0.010771839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8969741) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(0.2429602) q[2];
rz(0.36863676) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826137) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(-0.36732236) q[0];
rz(2.3452554) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(0.64250362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655661) q[0];
sx q[0];
rz(-2.4867704) q[0];
sx q[0];
rz(-0.73237082) q[0];
rz(1.9733866) q[2];
sx q[2];
rz(-1.0559096) q[2];
sx q[2];
rz(1.3001315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8290303) q[1];
sx q[1];
rz(-2.1305269) q[1];
sx q[1];
rz(-2.1055431) q[1];
x q[2];
rz(3.0152937) q[3];
sx q[3];
rz(-1.0923315) q[3];
sx q[3];
rz(-1.0431521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(1.7712234) q[2];
rz(0.227452) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(1.1915709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(-1.9434209) q[0];
rz(2.0612969) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(0.094873039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8137275) q[0];
sx q[0];
rz(-1.828232) q[0];
sx q[0];
rz(-1.0125005) q[0];
rz(2.2302365) q[2];
sx q[2];
rz(-1.9007287) q[2];
sx q[2];
rz(1.0418237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.171189) q[1];
sx q[1];
rz(-0.53817777) q[1];
sx q[1];
rz(1.6506877) q[1];
x q[2];
rz(2.8646144) q[3];
sx q[3];
rz(-1.2001032) q[3];
sx q[3];
rz(-2.1408368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0723116) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(2.7351232) q[2];
rz(-1.1786002) q[3];
sx q[3];
rz(-1.7013763) q[3];
sx q[3];
rz(1.9788205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793295) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(-2.80559) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(3.0373108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2662243) q[0];
sx q[0];
rz(-1.2422529) q[0];
sx q[0];
rz(-1.203241) q[0];
rz(-pi) q[1];
rz(3.0523446) q[2];
sx q[2];
rz(-0.1570905) q[2];
sx q[2];
rz(-2.2721827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9469918) q[1];
sx q[1];
rz(-0.90819383) q[1];
sx q[1];
rz(-2.6688982) q[1];
x q[2];
rz(-0.23747344) q[3];
sx q[3];
rz(-2.6769612) q[3];
sx q[3];
rz(0.97685087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.1445507) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(-1.087629) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8892141) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(-2.8401781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6431993) q[0];
sx q[0];
rz(-1.3424771) q[0];
sx q[0];
rz(-0.96781815) q[0];
rz(-pi) q[1];
rz(1.7784987) q[2];
sx q[2];
rz(-1.1900717) q[2];
sx q[2];
rz(2.8247339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6993048) q[1];
sx q[1];
rz(-1.3584474) q[1];
sx q[1];
rz(-1.9599171) q[1];
rz(0.88121342) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(2.415433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-3.1346698) q[2];
rz(2.210468) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(1.1091703) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329426) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(1.5337926) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(-2.839397) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9990197) q[0];
sx q[0];
rz(-0.66830615) q[0];
sx q[0];
rz(-1.4060941) q[0];
x q[1];
rz(2.4259461) q[2];
sx q[2];
rz(-1.3001912) q[2];
sx q[2];
rz(1.9220966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8540878) q[1];
sx q[1];
rz(-0.64784986) q[1];
sx q[1];
rz(-1.2811529) q[1];
rz(-pi) q[2];
rz(1.2421397) q[3];
sx q[3];
rz(-2.9813926) q[3];
sx q[3];
rz(-0.70806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89221421) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(0.92602473) q[2];
rz(1.9269491) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(-0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183384) q[0];
sx q[0];
rz(-0.76699081) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(2.8077937) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(2.0557859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1340116) q[0];
sx q[0];
rz(-1.7756878) q[0];
sx q[0];
rz(-1.6519068) q[0];
rz(-pi) q[1];
x q[1];
rz(1.169431) q[2];
sx q[2];
rz(-1.737672) q[2];
sx q[2];
rz(2.3292975) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9232417) q[1];
sx q[1];
rz(-2.3585547) q[1];
sx q[1];
rz(-2.3438575) q[1];
rz(1.7560739) q[3];
sx q[3];
rz(-2.4967125) q[3];
sx q[3];
rz(-1.4035937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51599017) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(-1.9773352) q[2];
rz(-0.26816756) q[3];
sx q[3];
rz(-1.8987013) q[3];
sx q[3];
rz(0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(1.9926158) q[0];
rz(0.2941429) q[1];
sx q[1];
rz(-1.5584757) q[1];
sx q[1];
rz(1.0424967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71871686) q[0];
sx q[0];
rz(-1.2407173) q[0];
sx q[0];
rz(-2.3407558) q[0];
rz(-pi) q[1];
rz(-2.9127321) q[2];
sx q[2];
rz(-2.3386152) q[2];
sx q[2];
rz(-0.91509089) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1139026) q[1];
sx q[1];
rz(-1.2053145) q[1];
sx q[1];
rz(-1.384907) q[1];
x q[2];
rz(-1.5339666) q[3];
sx q[3];
rz(-0.66434089) q[3];
sx q[3];
rz(-3.1337332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.7768804) q[2];
rz(2.4890066) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6259916) q[0];
sx q[0];
rz(-2.7643272) q[0];
sx q[0];
rz(-0.01734497) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(2.9877072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34129225) q[0];
sx q[0];
rz(-1.253739) q[0];
sx q[0];
rz(-1.4912259) q[0];
x q[1];
rz(-0.71435228) q[2];
sx q[2];
rz(-1.0754943) q[2];
sx q[2];
rz(2.60204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8887043) q[1];
sx q[1];
rz(-1.9333464) q[1];
sx q[1];
rz(1.8860555) q[1];
x q[2];
rz(3.1176223) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(2.6872203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-0.81622684) q[2];
sx q[2];
rz(0.44450644) q[2];
rz(2.9514173) q[3];
sx q[3];
rz(-1.5835652) q[3];
sx q[3];
rz(1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.001215) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(-2.0709399) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(0.79107034) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5673253) q[0];
sx q[0];
rz(-1.9057894) q[0];
sx q[0];
rz(0.092760249) q[0];
rz(-pi) q[1];
rz(2.5612381) q[2];
sx q[2];
rz(-2.6761746) q[2];
sx q[2];
rz(0.67829715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6562499) q[1];
sx q[1];
rz(-2.5963077) q[1];
sx q[1];
rz(2.4052909) q[1];
x q[2];
rz(-1.3338575) q[3];
sx q[3];
rz(-2.624369) q[3];
sx q[3];
rz(2.5631529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(1.6688639) q[2];
rz(-1.5276927) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(-2.6279502) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(2.6381941) q[2];
sx q[2];
rz(-1.361327) q[2];
sx q[2];
rz(0.8045902) q[2];
rz(-0.41889965) q[3];
sx q[3];
rz(-2.2875026) q[3];
sx q[3];
rz(0.84411375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
