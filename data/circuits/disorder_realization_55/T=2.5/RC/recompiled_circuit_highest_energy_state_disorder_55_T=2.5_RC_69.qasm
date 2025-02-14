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
rz(-1.363938) q[0];
sx q[0];
rz(-2.7231556) q[0];
sx q[0];
rz(1.3016181) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(1.4459223) q[1];
sx q[1];
rz(12.583153) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7977475) q[0];
sx q[0];
rz(-0.36917403) q[0];
sx q[0];
rz(1.7789715) q[0];
x q[1];
rz(0.016525908) q[2];
sx q[2];
rz(-1.7503097) q[2];
sx q[2];
rz(1.5500039) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1396347) q[1];
sx q[1];
rz(-3.1345001) q[1];
sx q[1];
rz(-0.79801871) q[1];
rz(0.34256012) q[3];
sx q[3];
rz(-2.9800219) q[3];
sx q[3];
rz(1.2662966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5825384) q[2];
sx q[2];
rz(-0.6607008) q[2];
sx q[2];
rz(-1.5793229) q[2];
rz(2.9301379) q[3];
sx q[3];
rz(-3.1410757) q[3];
sx q[3];
rz(-2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.6120537) q[0];
sx q[0];
rz(-0.27612975) q[0];
sx q[0];
rz(1.3268693) q[0];
rz(0.57948411) q[1];
sx q[1];
rz(-3.1377628) q[1];
sx q[1];
rz(-0.63900596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746637) q[0];
sx q[0];
rz(-0.94024728) q[0];
sx q[0];
rz(0.85682822) q[0];
rz(-pi) q[1];
rz(-0.12095708) q[2];
sx q[2];
rz(-1.5541583) q[2];
sx q[2];
rz(-1.5543092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7182448) q[1];
sx q[1];
rz(-3.1210174) q[1];
sx q[1];
rz(-0.59845509) q[1];
rz(1.6577107) q[3];
sx q[3];
rz(-0.79223903) q[3];
sx q[3];
rz(-2.5681005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7063893) q[2];
sx q[2];
rz(-0.13627626) q[2];
sx q[2];
rz(1.603568) q[2];
rz(1.5609353) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8921709) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(-0.38145915) q[0];
rz(0.70746607) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(2.0170508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1046421) q[0];
sx q[0];
rz(-0.35250394) q[0];
sx q[0];
rz(-2.332325) q[0];
x q[1];
rz(-1.7922282) q[2];
sx q[2];
rz(-3.0237758) q[2];
sx q[2];
rz(-3.0750781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3055691) q[1];
sx q[1];
rz(-1.5587121) q[1];
sx q[1];
rz(3.0786242) q[1];
rz(-pi) q[2];
rz(-2.400983) q[3];
sx q[3];
rz(-2.4428074) q[3];
sx q[3];
rz(-1.5027588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4418929) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(-0.049467889) q[2];
rz(-0.60702819) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98168755) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(0.017177563) q[0];
rz(2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5471829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80205432) q[0];
sx q[0];
rz(-2.0996764) q[0];
sx q[0];
rz(2.2418029) q[0];
rz(1.3457005) q[2];
sx q[2];
rz(-1.0402816) q[2];
sx q[2];
rz(2.7599285) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3354825) q[1];
sx q[1];
rz(-1.6948189) q[1];
sx q[1];
rz(1.6308484) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4956216) q[3];
sx q[3];
rz(-1.5630018) q[3];
sx q[3];
rz(0.92168671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(-0.68224254) q[2];
rz(3.0864129) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76160112) q[0];
sx q[0];
rz(-2.70607) q[0];
sx q[0];
rz(-2.7165661) q[0];
rz(-1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(-2.3262598) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7777625) q[0];
sx q[0];
rz(-1.6319931) q[0];
sx q[0];
rz(1.5598925) q[0];
x q[1];
rz(-1.4311976) q[2];
sx q[2];
rz(-3.1181212) q[2];
sx q[2];
rz(-1.0271629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22742352) q[1];
sx q[1];
rz(-0.14761758) q[1];
sx q[1];
rz(-0.62995894) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32834239) q[3];
sx q[3];
rz(-0.36188618) q[3];
sx q[3];
rz(2.7992424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1347947) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(-1.6689782) q[2];
rz(-0.93004477) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(-2.3013733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(-1.7148788) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-2.5529824) q[1];
sx q[1];
rz(1.1013365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5752714) q[0];
sx q[0];
rz(-2.3744517) q[0];
sx q[0];
rz(1.1772035) q[0];
rz(-pi) q[1];
rz(2.1997994) q[2];
sx q[2];
rz(-2.859349) q[2];
sx q[2];
rz(1.6779225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7633865) q[1];
sx q[1];
rz(-0.15608938) q[1];
sx q[1];
rz(2.2264477) q[1];
x q[2];
rz(1.3188478) q[3];
sx q[3];
rz(-3.0130606) q[3];
sx q[3];
rz(2.4737308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71658984) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(1.8343743) q[2];
rz(-0.37846765) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-2.4881261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8398447) q[0];
sx q[0];
rz(-1.8740338) q[0];
sx q[0];
rz(2.2140455) q[0];
rz(-1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.5313139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5132234) q[0];
sx q[0];
rz(-1.5484865) q[0];
sx q[0];
rz(1.5324027) q[0];
x q[1];
rz(0.53801914) q[2];
sx q[2];
rz(-0.45390651) q[2];
sx q[2];
rz(1.5858142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5802637) q[1];
sx q[1];
rz(-0.12910832) q[1];
sx q[1];
rz(1.5498954) q[1];
rz(-pi) q[2];
rz(1.0449991) q[3];
sx q[3];
rz(-2.5657095) q[3];
sx q[3];
rz(-0.27920846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7808468) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(-1.8878262) q[2];
rz(-2.4474261) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(-2.9085801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6041782) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(1.0373212) q[0];
rz(-1.5327771) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(-1.6745837) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4983385) q[0];
sx q[0];
rz(-0.20659978) q[0];
sx q[0];
rz(-1.2363966) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0612512) q[2];
sx q[2];
rz(-2.8626277) q[2];
sx q[2];
rz(2.7921048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26997494) q[1];
sx q[1];
rz(-1.570357) q[1];
sx q[1];
rz(-1.5696708) q[1];
rz(-pi) q[2];
rz(1.476273) q[3];
sx q[3];
rz(-2.0243939) q[3];
sx q[3];
rz(1.9840553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1303225) q[2];
sx q[2];
rz(-2.932817) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-0.52055001) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(1.6827778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.7875824) q[0];
sx q[0];
rz(-0.0024604877) q[0];
sx q[0];
rz(1.3186697) q[0];
rz(-1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(-1.5444548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1653571) q[0];
sx q[0];
rz(-2.6476323) q[0];
sx q[0];
rz(2.8970092) q[0];
x q[1];
rz(-2.804432) q[2];
sx q[2];
rz(-0.68575194) q[2];
sx q[2];
rz(-2.8711988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.006598) q[1];
sx q[1];
rz(-1.2603575) q[1];
sx q[1];
rz(-0.18532345) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61267743) q[3];
sx q[3];
rz(-2.3542488) q[3];
sx q[3];
rz(1.2932216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80457193) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(2.9447832) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-2.9350023) q[3];
sx q[3];
rz(-2.9406252) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8404959) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(1.949973) q[0];
rz(-1.6169029) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5651388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7708262) q[0];
sx q[0];
rz(-1.8642555) q[0];
sx q[0];
rz(1.3956885) q[0];
rz(1.1168408) q[2];
sx q[2];
rz(-1.7773787) q[2];
sx q[2];
rz(1.1517186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2597771) q[1];
sx q[1];
rz(-1.5702973) q[1];
sx q[1];
rz(-3.1406443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11348806) q[3];
sx q[3];
rz(-1.4295414) q[3];
sx q[3];
rz(3.0761513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(-1.4584165) q[2];
rz(-3.1111187) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771582) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(1.5741813) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(-1.6363999) q[2];
sx q[2];
rz(-3.0658683) q[2];
sx q[2];
rz(0.22584596) q[2];
rz(-2.3248657) q[3];
sx q[3];
rz(-1.1882267) q[3];
sx q[3];
rz(-1.4878275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
