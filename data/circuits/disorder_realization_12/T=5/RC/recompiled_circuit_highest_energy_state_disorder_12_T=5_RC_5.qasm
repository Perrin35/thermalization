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
rz(-0.15859088) q[0];
sx q[0];
rz(-0.51704419) q[0];
sx q[0];
rz(1.3102732) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(1.9603536) q[1];
sx q[1];
rz(11.601762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63435455) q[0];
sx q[0];
rz(-2.1956081) q[0];
sx q[0];
rz(-1.7376729) q[0];
rz(-pi) q[1];
rz(-0.50268051) q[2];
sx q[2];
rz(-2.7516904) q[2];
sx q[2];
rz(2.1250059) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.010299461) q[1];
sx q[1];
rz(-0.84761606) q[1];
sx q[1];
rz(2.7406969) q[1];
rz(-0.89972605) q[3];
sx q[3];
rz(-1.4802855) q[3];
sx q[3];
rz(0.29202005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0059119314) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(-1.4646336) q[2];
rz(1.3095193) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(-2.8215698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784073) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(1.7742668) q[0];
rz(1.6525846) q[1];
sx q[1];
rz(-0.79843489) q[1];
sx q[1];
rz(-0.29261011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45918834) q[0];
sx q[0];
rz(-1.4727797) q[0];
sx q[0];
rz(-0.3353663) q[0];
rz(-pi) q[1];
rz(0.16186773) q[2];
sx q[2];
rz(-1.4449931) q[2];
sx q[2];
rz(0.51728546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1024178) q[1];
sx q[1];
rz(-2.8264649) q[1];
sx q[1];
rz(0.71690317) q[1];
rz(1.5047362) q[3];
sx q[3];
rz(-1.8365897) q[3];
sx q[3];
rz(0.059826033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1613529) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(0.18388595) q[2];
rz(2.6889177) q[3];
sx q[3];
rz(-1.6190745) q[3];
sx q[3];
rz(-3.0116853) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9420796) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(2.0623656) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(-1.315518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53213648) q[0];
sx q[0];
rz(-2.2972496) q[0];
sx q[0];
rz(1.491602) q[0];
rz(-pi) q[1];
rz(-2.6805515) q[2];
sx q[2];
rz(-0.88141841) q[2];
sx q[2];
rz(-2.346476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1208468) q[1];
sx q[1];
rz(-1.6625424) q[1];
sx q[1];
rz(1.8534907) q[1];
rz(-1.6571548) q[3];
sx q[3];
rz(-1.2691853) q[3];
sx q[3];
rz(-0.82179606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9341854) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(-0.54971131) q[2];
rz(0.011119757) q[3];
sx q[3];
rz(-0.79922262) q[3];
sx q[3];
rz(1.8196677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(-2.8380561) q[0];
rz(0.15829463) q[1];
sx q[1];
rz(-2.0469432) q[1];
sx q[1];
rz(-2.1319481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5925668) q[0];
sx q[0];
rz(-2.0011304) q[0];
sx q[0];
rz(-2.8293912) q[0];
rz(2.539297) q[2];
sx q[2];
rz(-1.3118852) q[2];
sx q[2];
rz(-2.7342755) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6827439) q[1];
sx q[1];
rz(-1.413133) q[1];
sx q[1];
rz(-0.47994061) q[1];
rz(-1.5499653) q[3];
sx q[3];
rz(-2.0685141) q[3];
sx q[3];
rz(-1.3286984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9369072) q[2];
sx q[2];
rz(-0.82131177) q[2];
sx q[2];
rz(-2.8724907) q[2];
rz(-1.6875632) q[3];
sx q[3];
rz(-1.5917835) q[3];
sx q[3];
rz(-1.5788186) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8708385) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(2.7794072) q[0];
rz(2.4902792) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(1.6168894) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10241881) q[0];
sx q[0];
rz(-0.96511594) q[0];
sx q[0];
rz(1.4132947) q[0];
rz(-0.94793041) q[2];
sx q[2];
rz(-1.6375223) q[2];
sx q[2];
rz(-2.3579576) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.406639) q[1];
sx q[1];
rz(-1.9291153) q[1];
sx q[1];
rz(2.1291158) q[1];
rz(-pi) q[2];
rz(-0.42585856) q[3];
sx q[3];
rz(-2.3506864) q[3];
sx q[3];
rz(0.87825852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1093971) q[2];
sx q[2];
rz(-2.1869662) q[2];
sx q[2];
rz(1.3738013) q[2];
rz(-2.7023756) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(-2.5325328) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0405025) q[0];
sx q[0];
rz(-0.73223615) q[0];
sx q[0];
rz(-2.8771583) q[0];
rz(2.4624372) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(-0.98495475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0822345) q[0];
sx q[0];
rz(-1.2684039) q[0];
sx q[0];
rz(1.9134269) q[0];
x q[1];
rz(2.1169871) q[2];
sx q[2];
rz(-2.6272079) q[2];
sx q[2];
rz(0.75203224) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9912274) q[1];
sx q[1];
rz(-2.5179407) q[1];
sx q[1];
rz(-0.6536478) q[1];
rz(-2.9365083) q[3];
sx q[3];
rz(-2.0049378) q[3];
sx q[3];
rz(-1.6380426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.039845) q[2];
sx q[2];
rz(-2.9797649) q[2];
sx q[2];
rz(-2.7092773) q[2];
rz(2.1281706) q[3];
sx q[3];
rz(-1.5347967) q[3];
sx q[3];
rz(-0.54949808) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48563114) q[0];
sx q[0];
rz(-2.7483342) q[0];
sx q[0];
rz(2.0320758) q[0];
rz(0.79311496) q[1];
sx q[1];
rz(-1.3536645) q[1];
sx q[1];
rz(1.5464334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6731747) q[0];
sx q[0];
rz(-2.1270848) q[0];
sx q[0];
rz(2.2493717) q[0];
rz(-pi) q[1];
rz(-1.3605453) q[2];
sx q[2];
rz(-2.2814007) q[2];
sx q[2];
rz(2.9637314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82779694) q[1];
sx q[1];
rz(-0.89875752) q[1];
sx q[1];
rz(0.25700493) q[1];
rz(-1.6768544) q[3];
sx q[3];
rz(-2.2844188) q[3];
sx q[3];
rz(0.80435637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3194797) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(-1.7721843) q[2];
rz(-2.2111514) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9419788) q[0];
sx q[0];
rz(-1.9604585) q[0];
sx q[0];
rz(0.30366316) q[0];
rz(-2.1619201) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(2.7365275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9054149) q[0];
sx q[0];
rz(-1.9472031) q[0];
sx q[0];
rz(1.3628528) q[0];
x q[1];
rz(-1.151565) q[2];
sx q[2];
rz(-1.4128601) q[2];
sx q[2];
rz(1.4106223) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0822337) q[1];
sx q[1];
rz(-0.59525604) q[1];
sx q[1];
rz(-2.8221376) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9396254) q[3];
sx q[3];
rz(-1.151556) q[3];
sx q[3];
rz(0.48101411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8651809) q[2];
sx q[2];
rz(-0.96651912) q[2];
sx q[2];
rz(-0.69918862) q[2];
rz(-2.3326323) q[3];
sx q[3];
rz(-1.304108) q[3];
sx q[3];
rz(2.8185524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.1572105) q[0];
sx q[0];
rz(0.13741563) q[0];
rz(-3.0925062) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79530966) q[0];
sx q[0];
rz(-2.6457743) q[0];
sx q[0];
rz(0.056255503) q[0];
x q[1];
rz(-1.5586583) q[2];
sx q[2];
rz(-1.3985139) q[2];
sx q[2];
rz(-2.7207295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5180099) q[1];
sx q[1];
rz(-1.2482263) q[1];
sx q[1];
rz(-2.0913893) q[1];
rz(-pi) q[2];
rz(1.4490332) q[3];
sx q[3];
rz(-1.6596438) q[3];
sx q[3];
rz(-0.73032398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55234838) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(-2.153896) q[2];
rz(0.47932953) q[3];
sx q[3];
rz(-1.5440116) q[3];
sx q[3];
rz(2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.669303) q[0];
sx q[0];
rz(-2.9030114) q[0];
sx q[0];
rz(-0.14078374) q[0];
rz(-0.36551481) q[1];
sx q[1];
rz(-2.3194158) q[1];
sx q[1];
rz(2.4929094) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1295531) q[0];
sx q[0];
rz(-0.7551935) q[0];
sx q[0];
rz(-2.6213403) q[0];
x q[1];
rz(-1.4541523) q[2];
sx q[2];
rz(-2.3320227) q[2];
sx q[2];
rz(0.67942373) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5488447) q[1];
sx q[1];
rz(-2.5429568) q[1];
sx q[1];
rz(-0.099912481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8126103) q[3];
sx q[3];
rz(-1.8886654) q[3];
sx q[3];
rz(2.1235025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4623798) q[2];
sx q[2];
rz(-2.1361735) q[2];
sx q[2];
rz(3.0688378) q[2];
rz(-0.66271979) q[3];
sx q[3];
rz(-1.3136256) q[3];
sx q[3];
rz(2.976118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543906) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(1.6064593) q[1];
sx q[1];
rz(-2.2538593) q[1];
sx q[1];
rz(-1.5182553) q[1];
rz(-2.9130878) q[2];
sx q[2];
rz(-1.6863556) q[2];
sx q[2];
rz(-0.40876331) q[2];
rz(1.0853416) q[3];
sx q[3];
rz(-0.74182908) q[3];
sx q[3];
rz(-0.49334455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
