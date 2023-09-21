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
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049126547) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(-1.5729088) q[0];
rz(-2.2533478) q[2];
sx q[2];
rz(-2.1324854) q[2];
sx q[2];
rz(-2.782926) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3323101) q[1];
sx q[1];
rz(-2.6039632) q[1];
sx q[1];
rz(0.54772954) q[1];
rz(1.1349036) q[3];
sx q[3];
rz(-0.60308686) q[3];
sx q[3];
rz(-3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(2.9623048) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-2.9002088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86357147) q[0];
sx q[0];
rz(-1.2153421) q[0];
sx q[0];
rz(1.9832719) q[0];
x q[1];
rz(-2.6163289) q[2];
sx q[2];
rz(-2.004945) q[2];
sx q[2];
rz(-2.1936552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5026958) q[1];
sx q[1];
rz(-0.91050373) q[1];
sx q[1];
rz(3.026282) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4199735) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(0.54102708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(1.9096411) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(1.0528475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360639) q[0];
sx q[0];
rz(-1.9331421) q[0];
sx q[0];
rz(0.58689582) q[0];
rz(-2.8430976) q[2];
sx q[2];
rz(-2.2224991) q[2];
sx q[2];
rz(-3.0534844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7317253) q[1];
sx q[1];
rz(-1.6907398) q[1];
sx q[1];
rz(2.9883595) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5925519) q[3];
sx q[3];
rz(-0.43020159) q[3];
sx q[3];
rz(0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(0.012399013) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5445697) q[0];
sx q[0];
rz(-1.835808) q[0];
sx q[0];
rz(-0.17770627) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.054748) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(0.2211406) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.826556) q[1];
sx q[1];
rz(-0.16780014) q[1];
sx q[1];
rz(1.6102953) q[1];
rz(1.719791) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(-2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-0.59590894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673698) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(0.09341021) q[0];
rz(-pi) q[1];
rz(2.4679568) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(-3.0631531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0101498) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(-2.7308982) q[1];
rz(2.6951407) q[3];
sx q[3];
rz(-0.40735301) q[3];
sx q[3];
rz(0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(2.5904783) q[2];
rz(-0.20714949) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8897032) q[0];
sx q[0];
rz(-1.9648874) q[0];
sx q[0];
rz(-0.15051145) q[0];
rz(-1.7252543) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(2.6442106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8413137) q[1];
sx q[1];
rz(-2.0804555) q[1];
sx q[1];
rz(1.0213486) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0764334) q[3];
sx q[3];
rz(-1.061073) q[3];
sx q[3];
rz(-2.5173957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39020145) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(-2.899535) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(0.40329969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6751854) q[0];
sx q[0];
rz(-1.2012321) q[0];
sx q[0];
rz(0.36375605) q[0];
x q[1];
rz(2.4629668) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(-2.9976171) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1937716) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(-3.139688) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7607479) q[3];
sx q[3];
rz(-1.5956094) q[3];
sx q[3];
rz(-2.1638889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(2.1408391) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.5085295) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(0.60044926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8884044) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(1.123239) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1152994) q[2];
sx q[2];
rz(-1.0458938) q[2];
sx q[2];
rz(-0.80034791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.820103) q[1];
sx q[1];
rz(-2.9235225) q[1];
sx q[1];
rz(-0.58184187) q[1];
rz(2.5478701) q[3];
sx q[3];
rz(-0.36168081) q[3];
sx q[3];
rz(0.42632521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(3.0252769) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(0.41752648) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672887) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(-0.051961016) q[0];
rz(0.1065503) q[2];
sx q[2];
rz(-1.8333734) q[2];
sx q[2];
rz(-3.0343461) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(-0.71868371) q[1];
rz(2.9478491) q[3];
sx q[3];
rz(-1.1469517) q[3];
sx q[3];
rz(1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4801165) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(-0.68558145) q[0];
rz(-2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(1.1313653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.92537) q[0];
sx q[0];
rz(-0.67514456) q[0];
sx q[0];
rz(-2.569909) q[0];
x q[1];
rz(-2.4360043) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(0.74598344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79216204) q[1];
sx q[1];
rz(-1.7948479) q[1];
sx q[1];
rz(-0.38778023) q[1];
rz(-pi) q[2];
rz(-2.0496619) q[3];
sx q[3];
rz(-1.0189971) q[3];
sx q[3];
rz(-0.61210504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
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
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-0.60317294) q[2];
sx q[2];
rz(-2.0576253) q[2];
sx q[2];
rz(-1.1228592) q[2];
rz(-0.86757579) q[3];
sx q[3];
rz(-2.0316364) q[3];
sx q[3];
rz(0.33566725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
