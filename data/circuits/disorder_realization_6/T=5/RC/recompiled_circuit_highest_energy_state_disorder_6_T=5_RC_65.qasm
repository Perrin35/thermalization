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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(-1.7888223) q[0];
rz(1.1846722) q[1];
sx q[1];
rz(-0.67263043) q[1];
sx q[1];
rz(1.5157359) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11213451) q[0];
sx q[0];
rz(-1.3728276) q[0];
sx q[0];
rz(2.9517118) q[0];
rz(-0.47907655) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(-2.3912663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3876311) q[1];
sx q[1];
rz(-1.921614) q[1];
sx q[1];
rz(-0.49855788) q[1];
rz(-pi) q[2];
rz(0.58324849) q[3];
sx q[3];
rz(-0.49213675) q[3];
sx q[3];
rz(-1.0517091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1379913) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(-2.9555964) q[2];
rz(-0.018639175) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587392) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(0.3183611) q[0];
rz(2.5045577) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(-2.6723518) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56457389) q[0];
sx q[0];
rz(-1.7697608) q[0];
sx q[0];
rz(-1.7451203) q[0];
rz(-pi) q[1];
rz(1.1950011) q[2];
sx q[2];
rz(-1.6489121) q[2];
sx q[2];
rz(2.5004435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7982895) q[1];
sx q[1];
rz(-1.8899916) q[1];
sx q[1];
rz(1.3317687) q[1];
rz(-1.5638698) q[3];
sx q[3];
rz(-2.4516461) q[3];
sx q[3];
rz(1.5010709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3176754) q[2];
sx q[2];
rz(-0.21024148) q[2];
sx q[2];
rz(-2.4483185) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(-2.0461953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410809) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(-0.47955036) q[0];
rz(-0.50411049) q[1];
sx q[1];
rz(-1.1481552) q[1];
sx q[1];
rz(0.058301059) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710799) q[0];
sx q[0];
rz(-1.3127749) q[0];
sx q[0];
rz(-2.6239028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68966663) q[2];
sx q[2];
rz(-1.9885157) q[2];
sx q[2];
rz(-2.6725685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2268277) q[1];
sx q[1];
rz(-2.119442) q[1];
sx q[1];
rz(2.9991722) q[1];
rz(-0.50595768) q[3];
sx q[3];
rz(-1.3383597) q[3];
sx q[3];
rz(0.94601226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38640675) q[2];
sx q[2];
rz(-1.1815973) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(-1.6218328) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(2.7601385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-1.5244124) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(2.4066511) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5552444) q[0];
sx q[0];
rz(-2.1463094) q[0];
sx q[0];
rz(1.8481514) q[0];
rz(-1.3375912) q[2];
sx q[2];
rz(-2.2873023) q[2];
sx q[2];
rz(-0.75762123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0385602) q[1];
sx q[1];
rz(-1.2913229) q[1];
sx q[1];
rz(-0.36906645) q[1];
x q[2];
rz(-2.2913886) q[3];
sx q[3];
rz(-1.461004) q[3];
sx q[3];
rz(1.3444855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8526326) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(-3.0328879) q[2];
rz(-2.2147801) q[3];
sx q[3];
rz(-0.97065297) q[3];
sx q[3];
rz(1.577781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0573334) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(2.0507226) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(2.3172839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4949345) q[0];
sx q[0];
rz(-1.8047771) q[0];
sx q[0];
rz(-1.9569546) q[0];
rz(-0.53218158) q[2];
sx q[2];
rz(-2.144118) q[2];
sx q[2];
rz(-2.6087201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9059772) q[1];
sx q[1];
rz(-2.4433377) q[1];
sx q[1];
rz(-2.2119616) q[1];
rz(-pi) q[2];
rz(2.3184954) q[3];
sx q[3];
rz(-0.56328008) q[3];
sx q[3];
rz(2.7989504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6606286) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(0.29120293) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(-1.252482) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3785192) q[0];
sx q[0];
rz(-0.28994361) q[0];
sx q[0];
rz(-0.13105233) q[0];
rz(1.1669) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(0.022620591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0221828) q[0];
sx q[0];
rz(-1.2188605) q[0];
sx q[0];
rz(-0.58076136) q[0];
rz(-1.6940704) q[2];
sx q[2];
rz(-1.7825812) q[2];
sx q[2];
rz(-0.47763863) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5675332) q[1];
sx q[1];
rz(-2.1173734) q[1];
sx q[1];
rz(2.778227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9789627) q[3];
sx q[3];
rz(-1.3803326) q[3];
sx q[3];
rz(2.5951727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.268198) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(2.0568636) q[2];
rz(1.1449413) q[3];
sx q[3];
rz(-1.2866674) q[3];
sx q[3];
rz(2.1196608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5522083) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(0.4735221) q[0];
rz(-1.9206958) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(1.268505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4421688) q[0];
sx q[0];
rz(-2.0887362) q[0];
sx q[0];
rz(0.59364001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4642205) q[2];
sx q[2];
rz(-0.8607792) q[2];
sx q[2];
rz(-1.3959007) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0233588) q[1];
sx q[1];
rz(-1.9874554) q[1];
sx q[1];
rz(-0.55673846) q[1];
rz(-pi) q[2];
rz(1.1348261) q[3];
sx q[3];
rz(-1.4957168) q[3];
sx q[3];
rz(-1.0385385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(-0.24641985) q[2];
rz(-0.94270802) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(-2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3281658) q[0];
sx q[0];
rz(-1.270371) q[0];
sx q[0];
rz(0.060591977) q[0];
rz(1.6391485) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(1.5938119) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0133724) q[0];
sx q[0];
rz(-0.0039847535) q[0];
sx q[0];
rz(-0.26040034) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5752139) q[2];
sx q[2];
rz(-1.7292495) q[2];
sx q[2];
rz(2.2666933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9309153) q[1];
sx q[1];
rz(-1.1338755) q[1];
sx q[1];
rz(2.2411437) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4259858) q[3];
sx q[3];
rz(-1.2168435) q[3];
sx q[3];
rz(-0.41381124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26555201) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(1.6843686) q[2];
rz(3.0217116) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(1.1353528) q[0];
rz(0.10336939) q[1];
sx q[1];
rz(-1.4048978) q[1];
sx q[1];
rz(-0.9476544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4554868) q[0];
sx q[0];
rz(-0.64160937) q[0];
sx q[0];
rz(0.11193402) q[0];
rz(2.060076) q[2];
sx q[2];
rz(-0.7502509) q[2];
sx q[2];
rz(-0.81847755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3909437) q[1];
sx q[1];
rz(-0.90439561) q[1];
sx q[1];
rz(0.34088582) q[1];
x q[2];
rz(-0.74741628) q[3];
sx q[3];
rz(-1.0387762) q[3];
sx q[3];
rz(-0.72891392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7171628) q[2];
sx q[2];
rz(-2.5134176) q[2];
sx q[2];
rz(-0.90432811) q[2];
rz(1.2633911) q[3];
sx q[3];
rz(-1.2369913) q[3];
sx q[3];
rz(-2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(-2.1562321) q[0];
rz(-0.35161463) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(0.34686372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37736402) q[0];
sx q[0];
rz(-2.1954311) q[0];
sx q[0];
rz(-2.8714116) q[0];
rz(-3.0290481) q[2];
sx q[2];
rz(-2.0142609) q[2];
sx q[2];
rz(0.36641589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4046487) q[1];
sx q[1];
rz(-2.258806) q[1];
sx q[1];
rz(1.6860123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68315398) q[3];
sx q[3];
rz(-0.54617631) q[3];
sx q[3];
rz(2.7462296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2615307) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(3.0066709) q[2];
rz(0.12923446) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.1226817) q[0];
sx q[0];
rz(-1.8862579) q[0];
sx q[0];
rz(0.12611623) q[0];
rz(0.46161721) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(-2.2444637) q[2];
sx q[2];
rz(-1.9403602) q[2];
sx q[2];
rz(-0.42862949) q[2];
rz(0.25855385) q[3];
sx q[3];
rz(-1.7333442) q[3];
sx q[3];
rz(2.6556591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
