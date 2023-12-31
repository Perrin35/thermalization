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
rz(-1.1480968) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05605927) q[0];
sx q[0];
rz(-2.8319608) q[0];
sx q[0];
rz(3.1349896) q[0];
rz(-0.78656466) q[2];
sx q[2];
rz(-0.85430356) q[2];
sx q[2];
rz(2.5094945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4268036) q[1];
sx q[1];
rz(-2.0232632) q[1];
sx q[1];
rz(-1.2697551) q[1];
rz(2.8586219) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(0.24138385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5854908) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(2.7566064) q[0];
x q[1];
rz(-2.6163289) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(2.1936552) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0028210359) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-pi) q[2];
rz(1.7216191) q[3];
sx q[3];
rz(-1.5966468) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-1.0528475) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171779) q[0];
sx q[0];
rz(-0.6783692) q[0];
sx q[0];
rz(2.5413187) q[0];
x q[1];
rz(1.9387248) q[2];
sx q[2];
rz(-2.433948) q[2];
sx q[2];
rz(-0.55756535) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(-0.66837515) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1316109) q[3];
sx q[3];
rz(-1.1407033) q[3];
sx q[3];
rz(-0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(-0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(-2.9052177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.198092) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(0.99347465) q[0];
x q[1];
rz(-0.3503372) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(-1.9621153) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35509767) q[1];
sx q[1];
rz(-1.4031283) q[1];
sx q[1];
rz(3.1349036) q[1];
x q[2];
rz(0.48913042) q[3];
sx q[3];
rz(-1.4390579) q[3];
sx q[3];
rz(0.040442332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0376315) q[0];
sx q[0];
rz(-0.2346633) q[0];
sx q[0];
rz(-1.9734567) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4679568) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(0.078439586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039918598) q[1];
sx q[1];
rz(-0.44821757) q[1];
sx q[1];
rz(-0.43804534) q[1];
rz(-pi) q[2];
rz(-0.44645198) q[3];
sx q[3];
rz(-0.40735301) q[3];
sx q[3];
rz(-2.8706467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(-0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(2.8463083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(2.9910812) q[0];
rz(-pi) q[1];
rz(2.8304808) q[2];
sx q[2];
rz(-1.4236441) q[2];
sx q[2];
rz(1.0263024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8413137) q[1];
sx q[1];
rz(-1.0611371) q[1];
sx q[1];
rz(1.0213486) q[1];
rz(2.0814249) q[3];
sx q[3];
rz(-1.5139297) q[3];
sx q[3];
rz(0.97842723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.6823403) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(0.40329969) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.4629668) q[2];
sx q[2];
rz(-1.7947065) q[2];
sx q[2];
rz(-2.9976171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(-1.7809479) q[1];
rz(-pi) q[2];
rz(-1.7607479) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(2.1638889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91281259) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.5085295) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143215) q[0];
sx q[0];
rz(-1.1724171) q[0];
sx q[0];
rz(2.6398753) q[0];
rz(-pi) q[1];
rz(-1.5254283) q[2];
sx q[2];
rz(-0.52549911) q[2];
sx q[2];
rz(0.85277992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9143876) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(-1.449613) q[1];
rz(2.5478701) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(-0.42632521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(-2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0644181) q[0];
sx q[0];
rz(-1.525997) q[0];
sx q[0];
rz(-1.0391462) q[0];
rz(-pi) q[1];
rz(1.19403) q[2];
sx q[2];
rz(-0.28290877) q[2];
sx q[2];
rz(2.6434968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.233404) q[1];
sx q[1];
rz(-1.7618124) q[1];
sx q[1];
rz(1.4008271) q[1];
rz(2.0017654) q[3];
sx q[3];
rz(-1.3943854) q[3];
sx q[3];
rz(-0.21709066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(-2.8441692) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81998527) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(2.5489775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2106902) q[2];
sx q[2];
rz(-0.97133342) q[2];
sx q[2];
rz(0.42720366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8611421) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(0.54235561) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5351742) q[3];
sx q[3];
rz(-1.1675721) q[3];
sx q[3];
rz(-1.917178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-0.21223016) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(0.75081667) q[2];
sx q[2];
rz(-0.75559794) q[2];
sx q[2];
rz(-2.0970367) q[2];
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
