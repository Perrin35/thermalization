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
rz(1.2020943) q[0];
sx q[0];
rz(10.572875) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(-1.3936477) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(1.5729088) q[0];
x q[1];
rz(-2.4601828) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(-1.5208706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3323101) q[1];
sx q[1];
rz(-0.53762943) q[1];
sx q[1];
rz(-0.54772954) q[1];
x q[2];
rz(2.8586219) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(-2.9623048) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(0.24138385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(-0.38498621) q[0];
rz(-pi) q[1];
rz(-2.6163289) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(-0.94793749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6893377) q[1];
sx q[1];
rz(-0.66879767) q[1];
sx q[1];
rz(-1.4237088) q[1];
rz(-pi) q[2];
rz(0.026147141) q[3];
sx q[3];
rz(-1.7215683) q[3];
sx q[3];
rz(-2.1157516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0138578) q[2];
sx q[2];
rz(-1.8227791) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(1.7539304) q[3];
sx q[3];
rz(-2.6219086) q[3];
sx q[3];
rz(1.2319516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(1.0528475) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171779) q[0];
sx q[0];
rz(-0.6783692) q[0];
sx q[0];
rz(-2.5413187) q[0];
rz(-pi) q[1];
rz(-2.244433) q[2];
sx q[2];
rz(-1.3348012) q[2];
sx q[2];
rz(-1.8434075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3215072) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(-0.66837515) q[1];
rz(-pi) q[2];
rz(1.1406844) q[3];
sx q[3];
rz(-1.579869) q[3];
sx q[3];
rz(0.7113925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(0.23637493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94350067) q[0];
sx q[0];
rz(-2.8236832) q[0];
sx q[0];
rz(0.99347465) q[0];
x q[1];
rz(0.3503372) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(-1.1794773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9247777) q[1];
sx q[1];
rz(-1.5642011) q[1];
sx q[1];
rz(1.4031246) q[1];
rz(2.6524622) q[3];
sx q[3];
rz(-1.7025347) q[3];
sx q[3];
rz(0.040442332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699162) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0674954) q[0];
sx q[0];
rz(-1.4795545) q[0];
sx q[0];
rz(-1.3542961) q[0];
rz(-0.67363588) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(-3.0631531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1314428) q[1];
sx q[1];
rz(-1.7556659) q[1];
sx q[1];
rz(-2.7308982) q[1];
x q[2];
rz(-0.44645198) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(-0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(0.20714949) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2518894) q[0];
sx q[0];
rz(-1.9648874) q[0];
sx q[0];
rz(-0.15051145) q[0];
rz(-pi) q[1];
rz(1.7252543) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(-0.49738202) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5807647) q[1];
sx q[1];
rz(-2.0441214) q[1];
sx q[1];
rz(-2.5614489) q[1];
x q[2];
rz(0.065159273) q[3];
sx q[3];
rz(-1.061073) q[3];
sx q[3];
rz(0.62419696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(2.2078216) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-2.738293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1093275) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(-1.9637252) q[0];
x q[1];
rz(-2.4629668) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(-0.14397552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1937716) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(-0.0019046849) q[1];
x q[2];
rz(-1.3808448) q[3];
sx q[3];
rz(-1.5956094) q[3];
sx q[3];
rz(-0.97770377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91281259) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-2.6312857) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1338761) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4279815) q[0];
sx q[0];
rz(-2.5117154) q[0];
sx q[0];
rz(2.4226818) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1152994) q[2];
sx q[2];
rz(-1.0458938) q[2];
sx q[2];
rz(2.3412447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3214896) q[1];
sx q[1];
rz(-2.9235225) q[1];
sx q[1];
rz(-2.5597508) q[1];
rz(-1.3622215) q[3];
sx q[3];
rz(-1.2731291) q[3];
sx q[3];
rz(-2.0901594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-0.41752648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.674304) q[0];
sx q[0];
rz(-1.0397362) q[0];
sx q[0];
rz(0.051961016) q[0];
rz(-1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(0.49809581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7716277) q[1];
sx q[1];
rz(-1.7376448) q[1];
sx q[1];
rz(2.9478527) q[1];
x q[2];
rz(0.19374356) q[3];
sx q[3];
rz(-1.1469517) q[3];
sx q[3];
rz(1.8684052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(-2.8441692) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(2.0102274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3216074) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(-0.59261514) q[0];
rz(-pi) q[1];
rz(0.70558833) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(-2.3956092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2724534) q[1];
sx q[1];
rz(-1.9483882) q[1];
sx q[1];
rz(-1.3294405) q[1];
x q[2];
rz(-2.0496619) q[3];
sx q[3];
rz(-1.0189971) q[3];
sx q[3];
rz(2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(0.70739174) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
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