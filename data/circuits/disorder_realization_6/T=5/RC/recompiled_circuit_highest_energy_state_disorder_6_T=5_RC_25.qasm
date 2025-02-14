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
rz(1.037984) q[0];
sx q[0];
rz(-1.4572287) q[0];
sx q[0];
rz(1.7888223) q[0];
rz(-5.0985131) q[1];
sx q[1];
rz(2.4689622) q[1];
sx q[1];
rz(11.050635) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0294581) q[0];
sx q[0];
rz(-1.7687651) q[0];
sx q[0];
rz(-2.9517118) q[0];
rz(-pi) q[1];
rz(2.6625161) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(0.75032633) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7616854) q[1];
sx q[1];
rz(-0.60098472) q[1];
sx q[1];
rz(-0.65324776) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58324849) q[3];
sx q[3];
rz(-2.6494559) q[3];
sx q[3];
rz(1.0517091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1379913) q[2];
sx q[2];
rz(-0.10819745) q[2];
sx q[2];
rz(2.9555964) q[2];
rz(0.018639175) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(-2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587392) q[0];
sx q[0];
rz(-0.56621972) q[0];
sx q[0];
rz(-0.3183611) q[0];
rz(0.63703498) q[1];
sx q[1];
rz(-0.77992264) q[1];
sx q[1];
rz(0.46924082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97142727) q[0];
sx q[0];
rz(-1.741647) q[0];
sx q[0];
rz(2.9396482) q[0];
rz(-pi) q[1];
rz(-0.083949371) q[2];
sx q[2];
rz(-1.1962039) q[2];
sx q[2];
rz(0.89886802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9904203) q[1];
sx q[1];
rz(-1.3440596) q[1];
sx q[1];
rz(2.8137035) q[1];
rz(1.5777228) q[3];
sx q[3];
rz(-2.4516461) q[3];
sx q[3];
rz(1.5010709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8239173) q[2];
sx q[2];
rz(-0.21024148) q[2];
sx q[2];
rz(-2.4483185) q[2];
rz(1.3200101) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(1.0953974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005118) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(-2.6620423) q[0];
rz(2.6374822) q[1];
sx q[1];
rz(-1.1481552) q[1];
sx q[1];
rz(0.058301059) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47051278) q[0];
sx q[0];
rz(-1.3127749) q[0];
sx q[0];
rz(-0.51768984) q[0];
x q[1];
rz(1.0487172) q[2];
sx q[2];
rz(-0.95013844) q[2];
sx q[2];
rz(1.716937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95851499) q[1];
sx q[1];
rz(-2.5766008) q[1];
sx q[1];
rz(1.7989669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8350868) q[3];
sx q[3];
rz(-2.0619146) q[3];
sx q[3];
rz(-0.49784218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38640675) q[2];
sx q[2];
rz(-1.1815973) q[2];
sx q[2];
rz(0.03726658) q[2];
rz(-1.6218328) q[3];
sx q[3];
rz(-2.683679) q[3];
sx q[3];
rz(0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94654361) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(-0.066548912) q[0];
rz(-1.6171803) q[1];
sx q[1];
rz(-1.2434554) q[1];
sx q[1];
rz(-0.7349416) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682869) q[0];
sx q[0];
rz(-2.5096008) q[0];
sx q[0];
rz(2.7422264) q[0];
rz(1.8040015) q[2];
sx q[2];
rz(-2.2873023) q[2];
sx q[2];
rz(-0.75762123) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0543136) q[1];
sx q[1];
rz(-0.4590408) q[1];
sx q[1];
rz(2.4695818) q[1];
x q[2];
rz(-1.4052584) q[3];
sx q[3];
rz(-2.41417) q[3];
sx q[3];
rz(-0.35044985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8526326) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(0.10870474) q[2];
rz(-2.2147801) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(-1.577781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0573334) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(-1.0908701) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(2.3172839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4949345) q[0];
sx q[0];
rz(-1.3368155) q[0];
sx q[0];
rz(1.9569546) q[0];
rz(2.213843) q[2];
sx q[2];
rz(-1.130419) q[2];
sx q[2];
rz(-1.7945031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5368176) q[1];
sx q[1];
rz(-2.1120434) q[1];
sx q[1];
rz(-0.46525915) q[1];
rz(-1.1371255) q[3];
sx q[3];
rz(-1.1992362) q[3];
sx q[3];
rz(0.56321689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.480964) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(-0.29120293) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(1.8891107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3785192) q[0];
sx q[0];
rz(-0.28994361) q[0];
sx q[0];
rz(0.13105233) q[0];
rz(-1.9746926) q[1];
sx q[1];
rz(-2.1578372) q[1];
sx q[1];
rz(3.1189721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9126836) q[0];
sx q[0];
rz(-1.0297518) q[0];
sx q[0];
rz(1.9846657) q[0];
rz(-0.51949595) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(-1.0102538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5675332) q[1];
sx q[1];
rz(-1.0242192) q[1];
sx q[1];
rz(2.778227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87225391) q[3];
sx q[3];
rz(-0.24980751) q[3];
sx q[3];
rz(1.8810617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.268198) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(2.0568636) q[2];
rz(1.1449413) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(-2.1196608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5893843) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(-0.4735221) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-2.7291606) q[1];
sx q[1];
rz(-1.8730877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354123) q[0];
sx q[0];
rz(-2.0783193) q[0];
sx q[0];
rz(2.1730459) q[0];
rz(-0.67737214) q[2];
sx q[2];
rz(-2.2808135) q[2];
sx q[2];
rz(-1.3959007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4422677) q[1];
sx q[1];
rz(-1.066477) q[1];
sx q[1];
rz(2.0513351) q[1];
rz(-pi) q[2];
rz(-1.394519) q[3];
sx q[3];
rz(-2.699614) q[3];
sx q[3];
rz(0.37261841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.97741693) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(2.8951728) q[2];
rz(2.1988846) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(-2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134269) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(-0.060591977) q[0];
rz(1.6391485) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(-1.5477808) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4386182) q[0];
sx q[0];
rz(-1.5718223) q[0];
sx q[0];
rz(-0.0038504168) q[0];
x q[1];
rz(-0.28569371) q[2];
sx q[2];
rz(-0.59425747) q[2];
sx q[2];
rz(2.2069634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2912078) q[1];
sx q[1];
rz(-2.3602848) q[1];
sx q[1];
rz(-0.92618295) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71560685) q[3];
sx q[3];
rz(-1.9247492) q[3];
sx q[3];
rz(-0.41381124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8760406) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(-1.6843686) q[2];
rz(3.0217116) q[3];
sx q[3];
rz(-0.20457743) q[3];
sx q[3];
rz(-2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(-1.1353528) q[0];
rz(-0.10336939) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(-0.9476544) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1670939) q[0];
sx q[0];
rz(-1.5038953) q[0];
sx q[0];
rz(2.5029906) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41288167) q[2];
sx q[2];
rz(-0.92501174) q[2];
sx q[2];
rz(-1.4476763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7506489) q[1];
sx q[1];
rz(-2.237197) q[1];
sx q[1];
rz(-2.8007068) q[1];
rz(-pi) q[2];
rz(-0.89449785) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(1.2813527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4244298) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(2.2372645) q[2];
rz(1.2633911) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271027) q[0];
sx q[0];
rz(-0.75815433) q[0];
sx q[0];
rz(2.1562321) q[0];
rz(-2.789978) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(-0.34686372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7642286) q[0];
sx q[0];
rz(-0.94616156) q[0];
sx q[0];
rz(-0.2701811) q[0];
rz(-1.3386334) q[2];
sx q[2];
rz(-2.684991) q[2];
sx q[2];
rz(-0.62397623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.092791768) q[1];
sx q[1];
rz(-1.6597224) q[1];
sx q[1];
rz(2.4503178) q[1];
rz(1.9371766) q[3];
sx q[3];
rz(-1.9854331) q[3];
sx q[3];
rz(-0.36568991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88006192) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(-3.0066709) q[2];
rz(3.0123582) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(-1.038704) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018910949) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(0.46161721) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(2.2444637) q[2];
sx q[2];
rz(-1.2012325) q[2];
sx q[2];
rz(2.7129632) q[2];
rz(0.57030525) q[3];
sx q[3];
rz(-2.8371596) q[3];
sx q[3];
rz(0.53573487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
