OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4190726) q[0];
sx q[0];
rz(-2.4604586) q[0];
sx q[0];
rz(1.0406159) q[0];
rz(2.430727) q[1];
sx q[1];
rz(3.791888) q[1];
sx q[1];
rz(11.314582) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2501564) q[0];
sx q[0];
rz(-2.6134239) q[0];
sx q[0];
rz(0.78869606) q[0];
rz(-pi) q[1];
rz(2.0712713) q[2];
sx q[2];
rz(-1.0756167) q[2];
sx q[2];
rz(1.6409846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2179778) q[1];
sx q[1];
rz(-1.5461951) q[1];
sx q[1];
rz(-2.4776211) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4667418) q[3];
sx q[3];
rz(-0.94776692) q[3];
sx q[3];
rz(2.2365776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17994443) q[2];
sx q[2];
rz(-1.0893931) q[2];
sx q[2];
rz(-1.3570448) q[2];
rz(2.0876136) q[3];
sx q[3];
rz(-1.6097924) q[3];
sx q[3];
rz(1.0668782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2422975) q[0];
sx q[0];
rz(-2.5488148) q[0];
sx q[0];
rz(1.5574667) q[0];
rz(-1.2552931) q[1];
sx q[1];
rz(-0.86304945) q[1];
sx q[1];
rz(-3.0994298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145139) q[0];
sx q[0];
rz(-1.9111484) q[0];
sx q[0];
rz(2.447489) q[0];
rz(-pi) q[1];
rz(0.57664906) q[2];
sx q[2];
rz(-2.0248785) q[2];
sx q[2];
rz(-0.25694914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0706341) q[1];
sx q[1];
rz(-1.8629304) q[1];
sx q[1];
rz(-0.16865428) q[1];
rz(-pi) q[2];
rz(-0.20037074) q[3];
sx q[3];
rz(-1.0215825) q[3];
sx q[3];
rz(-2.7436243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2178847) q[2];
sx q[2];
rz(-1.8929241) q[2];
sx q[2];
rz(-0.62151796) q[2];
rz(-2.7593241) q[3];
sx q[3];
rz(-1.1416124) q[3];
sx q[3];
rz(2.3069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650448) q[0];
sx q[0];
rz(-1.6355729) q[0];
sx q[0];
rz(-1.7682834) q[0];
rz(-2.7690167) q[1];
sx q[1];
rz(-2.617045) q[1];
sx q[1];
rz(3.0212044) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62154462) q[0];
sx q[0];
rz(-2.1839474) q[0];
sx q[0];
rz(0.91213062) q[0];
rz(-pi) q[1];
rz(-1.8502813) q[2];
sx q[2];
rz(-1.6095543) q[2];
sx q[2];
rz(-1.9732158) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98657246) q[1];
sx q[1];
rz(-0.84671686) q[1];
sx q[1];
rz(-1.1447722) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27311886) q[3];
sx q[3];
rz(-1.0550753) q[3];
sx q[3];
rz(1.7777594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4176214) q[2];
sx q[2];
rz(-1.3105023) q[2];
sx q[2];
rz(-1.2869147) q[2];
rz(-0.78479615) q[3];
sx q[3];
rz(-2.0311821) q[3];
sx q[3];
rz(-2.3634461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50412905) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(2.8959287) q[0];
rz(-3.1371112) q[1];
sx q[1];
rz(-2.5878398) q[1];
sx q[1];
rz(3.0779238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.132936) q[0];
sx q[0];
rz(-1.8171165) q[0];
sx q[0];
rz(0.66605796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0374299) q[2];
sx q[2];
rz(-1.2793102) q[2];
sx q[2];
rz(-1.9574036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6323327) q[1];
sx q[1];
rz(-2.2508988) q[1];
sx q[1];
rz(-2.9010309) q[1];
rz(-pi) q[2];
rz(0.44594619) q[3];
sx q[3];
rz(-1.3896488) q[3];
sx q[3];
rz(1.1904448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8513546) q[2];
sx q[2];
rz(-1.3685702) q[2];
sx q[2];
rz(-1.9264539) q[2];
rz(-2.3859207) q[3];
sx q[3];
rz(-2.1561421) q[3];
sx q[3];
rz(0.22172609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.82472411) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(-2.6487937) q[0];
rz(-0.97152501) q[1];
sx q[1];
rz(-1.8172455) q[1];
sx q[1];
rz(-0.18878254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.550866) q[0];
sx q[0];
rz(-1.8456689) q[0];
sx q[0];
rz(-1.7455533) q[0];
rz(-pi) q[1];
rz(-2.9813779) q[2];
sx q[2];
rz(-3.1268178) q[2];
sx q[2];
rz(-2.229634) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35313126) q[1];
sx q[1];
rz(-1.6470688) q[1];
sx q[1];
rz(1.0623054) q[1];
rz(-pi) q[2];
rz(-1.3195511) q[3];
sx q[3];
rz(-0.94175168) q[3];
sx q[3];
rz(-0.97240868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8652953) q[2];
sx q[2];
rz(-1.6388288) q[2];
sx q[2];
rz(3.0740645) q[2];
rz(-1.4589795) q[3];
sx q[3];
rz(-2.429481) q[3];
sx q[3];
rz(-2.0290831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35906288) q[0];
sx q[0];
rz(-1.1787865) q[0];
sx q[0];
rz(-1.4558526) q[0];
rz(0.23446941) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(0.04105982) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053692) q[0];
sx q[0];
rz(-1.5991044) q[0];
sx q[0];
rz(-0.7014277) q[0];
x q[1];
rz(-1.7975989) q[2];
sx q[2];
rz(-0.58355809) q[2];
sx q[2];
rz(-0.22554071) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1225374) q[1];
sx q[1];
rz(-1.1047939) q[1];
sx q[1];
rz(-0.48961498) q[1];
rz(0.84223522) q[3];
sx q[3];
rz(-1.4596617) q[3];
sx q[3];
rz(-2.5499217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91512758) q[2];
sx q[2];
rz(-2.3751986) q[2];
sx q[2];
rz(0.29092947) q[2];
rz(-2.9229524) q[3];
sx q[3];
rz(-0.69059697) q[3];
sx q[3];
rz(2.2712928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17756322) q[0];
sx q[0];
rz(-0.95528269) q[0];
sx q[0];
rz(-2.5780504) q[0];
rz(0.82849416) q[1];
sx q[1];
rz(-0.68873134) q[1];
sx q[1];
rz(0.73961893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1391382) q[0];
sx q[0];
rz(-2.1147303) q[0];
sx q[0];
rz(0.55940658) q[0];
x q[1];
rz(1.8663667) q[2];
sx q[2];
rz(-2.6049714) q[2];
sx q[2];
rz(-1.581096) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1465211) q[1];
sx q[1];
rz(-2.5634829) q[1];
sx q[1];
rz(2.0721295) q[1];
rz(-pi) q[2];
rz(1.4007934) q[3];
sx q[3];
rz(-0.94906607) q[3];
sx q[3];
rz(-2.9431557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.171611) q[2];
sx q[2];
rz(-2.3889399) q[2];
sx q[2];
rz(-0.29772154) q[2];
rz(0.026084829) q[3];
sx q[3];
rz(-1.9484768) q[3];
sx q[3];
rz(-2.397876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.96935) q[0];
sx q[0];
rz(-1.5765215) q[0];
sx q[0];
rz(1.4132389) q[0];
rz(0.50624943) q[1];
sx q[1];
rz(-1.0289611) q[1];
sx q[1];
rz(1.9821573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6613061) q[0];
sx q[0];
rz(-1.6723858) q[0];
sx q[0];
rz(2.8364968) q[0];
rz(-pi) q[1];
rz(1.5250823) q[2];
sx q[2];
rz(-2.1518143) q[2];
sx q[2];
rz(-2.8948262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22190801) q[1];
sx q[1];
rz(-1.3338998) q[1];
sx q[1];
rz(-0.48478957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6990463) q[3];
sx q[3];
rz(-0.82144605) q[3];
sx q[3];
rz(-1.3155703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96782288) q[2];
sx q[2];
rz(-2.2161667) q[2];
sx q[2];
rz(-1.7769163) q[2];
rz(0.81554282) q[3];
sx q[3];
rz(-1.848685) q[3];
sx q[3];
rz(-1.3968141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57104617) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(1.1676769) q[0];
rz(2.9171464) q[1];
sx q[1];
rz(-0.75643221) q[1];
sx q[1];
rz(-0.44463739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1261869) q[0];
sx q[0];
rz(-3.0660136) q[0];
sx q[0];
rz(0.74962693) q[0];
rz(0.860487) q[2];
sx q[2];
rz(-2.103946) q[2];
sx q[2];
rz(1.8955961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9079355) q[1];
sx q[1];
rz(-2.7747215) q[1];
sx q[1];
rz(1.9625825) q[1];
rz(-2.8520582) q[3];
sx q[3];
rz(-2.6752895) q[3];
sx q[3];
rz(2.795855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.726752) q[2];
sx q[2];
rz(-1.8291992) q[2];
sx q[2];
rz(1.7337588) q[2];
rz(1.5358745) q[3];
sx q[3];
rz(-1.9460257) q[3];
sx q[3];
rz(-2.2685952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6294412) q[0];
sx q[0];
rz(-0.81708556) q[0];
sx q[0];
rz(-2.0508811) q[0];
rz(-1.0461944) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(-2.2549021) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8532904) q[0];
sx q[0];
rz(-1.9294943) q[0];
sx q[0];
rz(-0.38889287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1524296) q[2];
sx q[2];
rz(-2.5083087) q[2];
sx q[2];
rz(3.0208602) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1613933) q[1];
sx q[1];
rz(-2.7100013) q[1];
sx q[1];
rz(-0.031567911) q[1];
x q[2];
rz(-2.2374152) q[3];
sx q[3];
rz(-1.2632367) q[3];
sx q[3];
rz(-0.22790652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1297168) q[2];
sx q[2];
rz(-1.0150212) q[2];
sx q[2];
rz(1.2474308) q[2];
rz(-0.85793197) q[3];
sx q[3];
rz(-1.6777638) q[3];
sx q[3];
rz(-1.5425382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3457376) q[0];
sx q[0];
rz(-1.416774) q[0];
sx q[0];
rz(-0.10373535) q[0];
rz(0.70032447) q[1];
sx q[1];
rz(-1.843597) q[1];
sx q[1];
rz(1.2784943) q[1];
rz(0.1968443) q[2];
sx q[2];
rz(-1.4509401) q[2];
sx q[2];
rz(0.069139253) q[2];
rz(2.2195001) q[3];
sx q[3];
rz(-2.5917883) q[3];
sx q[3];
rz(-2.5592309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
