OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(1.9934959) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.2611715) q[0];
sx q[0];
rz(-1.5729088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2533478) q[2];
sx q[2];
rz(-1.0091072) q[2];
sx q[2];
rz(0.35866666) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4268036) q[1];
sx q[1];
rz(-1.1183294) q[1];
sx q[1];
rz(1.8718375) q[1];
x q[2];
rz(2.0066891) q[3];
sx q[3];
rz(-0.60308686) q[3];
sx q[3];
rz(3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-2.9002088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(2.7566064) q[0];
rz(-pi) q[1];
rz(-1.0788467) q[2];
sx q[2];
rz(-1.098512) q[2];
sx q[2];
rz(-0.3837331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.452255) q[1];
sx q[1];
rz(-2.472795) q[1];
sx q[1];
rz(-1.4237088) q[1];
rz(-1.4003795) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(-1.1982329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-1.1068809) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2967865) q[0];
sx q[0];
rz(-2.1150981) q[0];
sx q[0];
rz(1.9980206) q[0];
x q[1];
rz(1.9387248) q[2];
sx q[2];
rz(-0.70764467) q[2];
sx q[2];
rz(0.55756535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9991418) q[1];
sx q[1];
rz(-1.72292) q[1];
sx q[1];
rz(-1.692148) q[1];
x q[2];
rz(3.1316109) q[3];
sx q[3];
rz(-1.1407033) q[3];
sx q[3];
rz(-0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1716487) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(1.3522211) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(-0.23637493) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073233152) q[0];
sx q[0];
rz(-1.7422361) q[0];
sx q[0];
rz(-1.839848) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0868446) q[2];
sx q[2];
rz(-1.2630672) q[2];
sx q[2];
rz(-0.2211406) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(3.1349036) q[1];
x q[2];
rz(-2.866719) q[3];
sx q[3];
rz(-0.50516869) q[3];
sx q[3];
rz(1.8531909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248557) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(3.0481824) q[0];
x q[1];
rz(-0.17275177) q[2];
sx q[2];
rz(-0.68095945) q[2];
sx q[2];
rz(1.3576042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1016741) q[1];
sx q[1];
rz(-2.6933751) q[1];
sx q[1];
rz(0.43804534) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44645198) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(-2.8706467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.15774396) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(2.8463083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(0.15051145) q[0];
rz(-1.4163383) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(2.6442106) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8413137) q[1];
sx q[1];
rz(-2.0804555) q[1];
sx q[1];
rz(-1.0213486) q[1];
rz(-1.0601677) q[3];
sx q[3];
rz(-1.5139297) q[3];
sx q[3];
rz(-2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.6823403) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(2.738293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4664073) q[0];
sx q[0];
rz(-1.2012321) q[0];
sx q[0];
rz(-2.7778366) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34801872) q[2];
sx q[2];
rz(-2.4325779) q[2];
sx q[2];
rz(1.6955171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94782103) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(0.0019046849) q[1];
rz(-pi) q[2];
rz(1.4401011) q[3];
sx q[3];
rz(-0.19154597) q[3];
sx q[3];
rz(-2.6768315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531882) q[0];
sx q[0];
rz(-2.0300403) q[0];
sx q[0];
rz(-1.123239) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1152994) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(-0.80034791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.820103) q[1];
sx q[1];
rz(-0.21807018) q[1];
sx q[1];
rz(-0.58184187) q[1];
rz(-pi) q[2];
rz(-0.30386691) q[3];
sx q[3];
rz(-1.3715203) q[3];
sx q[3];
rz(2.5602333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0042469) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(2.7159193) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-0.41752648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672887) q[0];
sx q[0];
rz(-1.0397362) q[0];
sx q[0];
rz(0.051961016) q[0];
rz(-pi) q[1];
rz(-1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(0.49809581) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.233404) q[1];
sx q[1];
rz(-1.3797803) q[1];
sx q[1];
rz(-1.7407655) q[1];
rz(-pi) q[2];
rz(-2.0017654) q[3];
sx q[3];
rz(-1.3943854) q[3];
sx q[3];
rz(-2.924502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2162227) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(-0.57168369) q[0];
rz(2.4360043) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(-0.74598344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28045052) q[1];
sx q[1];
rz(-2.6965953) q[1];
sx q[1];
rz(0.54235561) q[1];
x q[2];
rz(-2.5351742) q[3];
sx q[3];
rz(-1.9740205) q[3];
sx q[3];
rz(-1.2244146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512882) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-2.5384197) q[2];
sx q[2];
rz(-1.0839673) q[2];
sx q[2];
rz(2.0187335) q[2];
rz(-0.91602305) q[3];
sx q[3];
rz(-2.3229204) q[3];
sx q[3];
rz(-1.7182072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
