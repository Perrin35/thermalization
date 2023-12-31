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
rz(4.0822786) q[1];
sx q[1];
rz(10.818426) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5210261) q[0];
sx q[0];
rz(-1.5728083) q[0];
sx q[0];
rz(2.8319671) q[0];
rz(2.355028) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(-2.5094945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4268036) q[1];
sx q[1];
rz(-2.0232632) q[1];
sx q[1];
rz(1.2697551) q[1];
rz(2.8586219) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(-2.9623048) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0497465) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(-0.24138385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2780212) q[0];
sx q[0];
rz(-1.9262505) q[0];
sx q[0];
rz(1.1583207) q[0];
rz(-pi) q[1];
rz(2.6163289) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(0.94793749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0028210359) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(-2.2343193) q[1];
x q[2];
rz(-1.4199735) q[3];
sx q[3];
rz(-1.5966468) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.8227791) q[2];
sx q[2];
rz(-1.1068809) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(2.0887451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360639) q[0];
sx q[0];
rz(-1.2084506) q[0];
sx q[0];
rz(0.58689582) q[0];
rz(-0.89715965) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(-1.8434075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(0.66837515) q[1];
x q[2];
rz(1.5925519) q[3];
sx q[3];
rz(-0.43020159) q[3];
sx q[3];
rz(2.3019626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(-0.28918239) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(-0.23637493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94350067) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(-0.99347465) q[0];
x q[1];
rz(-2.7912555) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(1.9621153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35509767) q[1];
sx q[1];
rz(-1.4031283) q[1];
sx q[1];
rz(-0.0066890072) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27487367) q[3];
sx q[3];
rz(-0.50516869) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(-0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(-0.59590894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0376315) q[0];
sx q[0];
rz(-2.9069293) q[0];
sx q[0];
rz(-1.9734567) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17275177) q[2];
sx q[2];
rz(-0.68095945) q[2];
sx q[2];
rz(-1.7839884) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1314428) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(-2.7308982) q[1];
x q[2];
rz(2.6951407) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(2.8706467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(2.5904783) q[2];
rz(-2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.6023887) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(-0.47479182) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8897032) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(-2.9910812) q[0];
rz(-0.45093243) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(2.1692838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59776781) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(2.3901229) q[1];
rz(1.0601677) q[3];
sx q[3];
rz(-1.5139297) q[3];
sx q[3];
rz(2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(0.93377101) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-2.899535) q[0];
rz(-2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(-0.40329969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86352578) q[0];
sx q[0];
rz(-2.6289872) q[0];
sx q[0];
rz(0.82786064) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4629668) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(2.9976171) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93869081) q[1];
sx q[1];
rz(-0.21015973) q[1];
sx q[1];
rz(1.5797257) q[1];
x q[2];
rz(1.7014916) q[3];
sx q[3];
rz(-2.9500467) q[3];
sx q[3];
rz(-2.6768315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8884044) q[0];
sx q[0];
rz(-2.0300403) q[0];
sx q[0];
rz(2.0183536) q[0];
rz(-3.1152994) q[2];
sx q[2];
rz(-1.0458938) q[2];
sx q[2];
rz(-0.80034791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8200092) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(0.18305852) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30386691) q[3];
sx q[3];
rz(-1.3715203) q[3];
sx q[3];
rz(0.58135939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-11*pi/12) q[3];
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
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(-2.7240662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5696213) q[0];
sx q[0];
rz(-0.53335359) q[0];
sx q[0];
rz(1.4825975) q[0];
rz(-1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(-2.6434968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36996499) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(-0.19373993) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1398272) q[3];
sx q[3];
rz(-1.7472072) q[3];
sx q[3];
rz(-0.21709066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(2.647906) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(2.8216968) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17393728) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2162227) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(0.57168369) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70558833) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(-0.74598344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3494306) q[1];
sx q[1];
rz(-1.3467448) q[1];
sx q[1];
rz(0.38778023) q[1];
x q[2];
rz(-2.4990436) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(-2.97314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82548213) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(2.9293625) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(2.5384197) q[2];
sx q[2];
rz(-2.0576253) q[2];
sx q[2];
rz(-1.1228592) q[2];
rz(-0.57701941) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
