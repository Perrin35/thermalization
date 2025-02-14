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
rz(1.7776547) q[0];
sx q[0];
rz(-0.41843709) q[0];
sx q[0];
rz(1.8399746) q[0];
rz(0.16297451) q[1];
sx q[1];
rz(-1.4459223) q[1];
sx q[1];
rz(0.016782848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12110773) q[0];
sx q[0];
rz(-1.2099625) q[0];
sx q[0];
rz(0.079795436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7503337) q[2];
sx q[2];
rz(-1.5870567) q[2];
sx q[2];
rz(-0.017841466) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1396347) q[1];
sx q[1];
rz(-0.0070925709) q[1];
sx q[1];
rz(0.79801871) q[1];
x q[2];
rz(0.34256012) q[3];
sx q[3];
rz(-0.16157074) q[3];
sx q[3];
rz(-1.2662966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(-1.5622697) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(0.16099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(-1.3268693) q[0];
rz(-0.57948411) q[1];
sx q[1];
rz(-3.1377628) q[1];
sx q[1];
rz(0.63900596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746637) q[0];
sx q[0];
rz(-0.94024728) q[0];
sx q[0];
rz(0.85682822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5540358) q[2];
sx q[2];
rz(-1.6917366) q[2];
sx q[2];
rz(-0.018509381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.745805) q[1];
sx q[1];
rz(-1.5592054) q[1];
sx q[1];
rz(-0.017000217) q[1];
x q[2];
rz(1.483882) q[3];
sx q[3];
rz(-0.79223903) q[3];
sx q[3];
rz(-0.57349216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(1.603568) q[2];
rz(-1.5609353) q[3];
sx q[3];
rz(-0.014336421) q[3];
sx q[3];
rz(3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8921709) q[0];
sx q[0];
rz(-2.6264661) q[0];
sx q[0];
rz(-2.7601335) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(1.1245419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9455348) q[0];
sx q[0];
rz(-1.3302517) q[0];
sx q[0];
rz(-1.831016) q[0];
rz(-pi) q[1];
rz(1.7922282) q[2];
sx q[2];
rz(-0.11781684) q[2];
sx q[2];
rz(-3.0750781) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3055691) q[1];
sx q[1];
rz(-1.5587121) q[1];
sx q[1];
rz(3.0786242) q[1];
x q[2];
rz(-1.0550523) q[3];
sx q[3];
rz(-2.0655144) q[3];
sx q[3];
rz(-2.5124541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4418929) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(-0.049467889) q[2];
rz(-2.5345645) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(-1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.1599051) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(0.29302868) q[1];
sx q[1];
rz(-0.79048645) q[1];
sx q[1];
rz(-1.5944098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7538504) q[0];
sx q[0];
rz(-1.004129) q[0];
sx q[0];
rz(0.64105861) q[0];
x q[1];
rz(1.7958922) q[2];
sx q[2];
rz(-2.101311) q[2];
sx q[2];
rz(2.7599285) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7878927) q[1];
sx q[1];
rz(-0.13772923) q[1];
sx q[1];
rz(2.6929123) q[1];
rz(-pi) q[2];
rz(0.64597102) q[3];
sx q[3];
rz(-1.5630018) q[3];
sx q[3];
rz(2.2199059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25345099) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(-2.4593501) q[2];
rz(-3.0864129) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(-1.3050219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76160112) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(2.7165661) q[0];
rz(1.5399326) q[1];
sx q[1];
rz(-2.6585177) q[1];
sx q[1];
rz(-0.81533283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013215) q[0];
sx q[0];
rz(-0.062159408) q[0];
sx q[0];
rz(2.9654853) q[0];
x q[1];
rz(1.4311976) q[2];
sx q[2];
rz(-0.023471467) q[2];
sx q[2];
rz(2.1144298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9141691) q[1];
sx q[1];
rz(-2.9939751) q[1];
sx q[1];
rz(2.5116337) q[1];
x q[2];
rz(-0.34408042) q[3];
sx q[3];
rz(-1.4563784) q[3];
sx q[3];
rz(-0.92000577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0067979) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(1.6689782) q[2];
rz(-2.2115479) q[3];
sx q[3];
rz(-0.01447066) q[3];
sx q[3];
rz(-2.3013733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(1.4267138) q[0];
rz(-0.73000437) q[1];
sx q[1];
rz(-2.5529824) q[1];
sx q[1];
rz(1.1013365) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8555785) q[0];
sx q[0];
rz(-1.8402303) q[0];
sx q[0];
rz(-2.2982909) q[0];
rz(-0.16898245) q[2];
sx q[2];
rz(-1.3436396) q[2];
sx q[2];
rz(-2.3262466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8582272) q[1];
sx q[1];
rz(-1.4472597) q[1];
sx q[1];
rz(-3.0459411) q[1];
rz(-3.1093842) q[3];
sx q[3];
rz(-1.695249) q[3];
sx q[3];
rz(2.2197753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4250028) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(1.8343743) q[2];
rz(-0.37846765) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(0.65346658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8398447) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(0.92754716) q[0];
rz(-1.357366) q[1];
sx q[1];
rz(-2.3097242) q[1];
sx q[1];
rz(1.6102788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0831628) q[0];
sx q[0];
rz(-1.5324123) q[0];
sx q[0];
rz(-3.1192664) q[0];
rz(-0.39674098) q[2];
sx q[2];
rz(-1.3441685) q[2];
sx q[2];
rz(-2.6643348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.011259638) q[1];
sx q[1];
rz(-1.5734871) q[1];
sx q[1];
rz(-1.6998768) q[1];
rz(-pi) q[2];
x q[2];
rz(2.82656) q[3];
sx q[3];
rz(-1.0803534) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53741443) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(2.1042714) q[0];
rz(-1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(-1.467009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1572185) q[0];
sx q[0];
rz(-1.7658002) q[0];
sx q[0];
rz(0.068679811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.59378) q[2];
sx q[2];
rz(-1.8488374) q[2];
sx q[2];
rz(-2.8756623) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8407708) q[1];
sx q[1];
rz(-1.5719218) q[1];
sx q[1];
rz(-3.1411533) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45536228) q[3];
sx q[3];
rz(-1.6557367) q[3];
sx q[3];
rz(2.7698539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1303225) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-2.6210426) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540102) q[0];
sx q[0];
rz(-0.0024604877) q[0];
sx q[0];
rz(-1.3186697) q[0];
rz(1.7240546) q[1];
sx q[1];
rz(-0.28957614) q[1];
sx q[1];
rz(1.5971378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97623551) q[0];
sx q[0];
rz(-0.49396038) q[0];
sx q[0];
rz(2.8970092) q[0];
rz(-pi) q[1];
rz(0.6575281) q[2];
sx q[2];
rz(-1.7818461) q[2];
sx q[2];
rz(1.5653277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68543562) q[1];
sx q[1];
rz(-0.36000571) q[1];
sx q[1];
rz(2.092157) q[1];
x q[2];
rz(-0.61267743) q[3];
sx q[3];
rz(-0.78734382) q[3];
sx q[3];
rz(-1.8483711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3370207) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(-2.9447832) q[2];
rz(-1.192441) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(-0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8404959) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(1.949973) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(1.5651388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25116205) q[0];
sx q[0];
rz(-1.403247) q[0];
sx q[0];
rz(-0.29775374) q[0];
x q[1];
rz(-0.2290957) q[2];
sx q[2];
rz(-1.1271897) q[2];
sx q[2];
rz(-0.51883051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31101878) q[1];
sx q[1];
rz(-1.5717447) q[1];
sx q[1];
rz(1.5702973) q[1];
rz(2.2433167) q[3];
sx q[3];
rz(-2.9606323) q[3];
sx q[3];
rz(2.3955936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18371753) q[2];
sx q[2];
rz(-2.5536733) q[2];
sx q[2];
rz(1.6831762) q[2];
rz(-3.1111187) q[3];
sx q[3];
rz(-0.009549791) q[3];
sx q[3];
rz(-0.20235801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6644345) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(-1.5674113) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(1.6363999) q[2];
sx q[2];
rz(-0.075724307) q[2];
sx q[2];
rz(-2.9157467) q[2];
rz(2.3248657) q[3];
sx q[3];
rz(-1.9533659) q[3];
sx q[3];
rz(1.6537651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
