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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0855334) q[0];
sx q[0];
rz(-0.30963184) q[0];
sx q[0];
rz(0.0066030533) q[0];
rz(-0.88824484) q[2];
sx q[2];
rz(-2.1324854) q[2];
sx q[2];
rz(-0.35866666) q[2];
sx q[3];
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
rz(-1.8718375) q[1];
rz(-pi) q[2];
rz(-2.1288793) q[3];
sx q[3];
rz(-1.8126243) q[3];
sx q[3];
rz(-1.946132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(0.17928784) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-0.15549913) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.9002088) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3791083) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(-2.3178029) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74626211) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(1.2506739) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0028210359) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(-2.2343193) q[1];
x q[2];
rz(0.026147141) q[3];
sx q[3];
rz(-1.4200243) q[3];
sx q[3];
rz(2.1157516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-1.1068809) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(1.2319516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36104193) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-2.0887451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360639) q[0];
sx q[0];
rz(-1.9331421) q[0];
sx q[0];
rz(-0.58689582) q[0];
x q[1];
rz(-0.29849507) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(0.088108206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14245089) q[1];
sx q[1];
rz(-1.4186727) q[1];
sx q[1];
rz(1.692148) q[1];
rz(-pi) q[2];
rz(1.1406844) q[3];
sx q[3];
rz(-1.579869) q[3];
sx q[3];
rz(-2.4302002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(2.5946674) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.3522211) q[0];
rz(-2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94350067) q[0];
sx q[0];
rz(-2.8236832) q[0];
sx q[0];
rz(0.99347465) q[0];
rz(0.3503372) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(1.9621153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(3.1349036) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6524622) q[3];
sx q[3];
rz(-1.7025347) q[3];
sx q[3];
rz(3.1011503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(-2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6248557) q[0];
sx q[0];
rz(-1.3552109) q[0];
sx q[0];
rz(-0.09341021) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4679568) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(0.078439586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0101498) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(0.4106945) q[1];
x q[2];
rz(2.6951407) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(-0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(-2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-0.95170784) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5135358) q[0];
sx q[0];
rz(-0.4204458) q[0];
sx q[0];
rz(-1.224731) q[0];
rz(-pi) q[1];
rz(0.45093243) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(0.97230881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59776781) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(0.75146971) q[1];
rz(-pi) q[2];
rz(-1.4548364) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(0.93377101) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(2.899535) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(-2.738293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1093275) q[0];
sx q[0];
rz(-1.2326213) q[0];
sx q[0];
rz(-1.1778675) q[0];
x q[1];
rz(-1.8554011) q[2];
sx q[2];
rz(-2.2294514) q[2];
sx q[2];
rz(-1.2496442) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1937716) q[1];
sx q[1];
rz(-1.3606451) q[1];
sx q[1];
rz(-0.0019046849) q[1];
rz(0.025267406) q[3];
sx q[3];
rz(-1.380904) q[3];
sx q[3];
rz(-2.5437298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(2.6469321) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884044) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(-1.123239) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1152994) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(0.80034791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2272051) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(-1.6919796) q[1];
rz(0.30386691) q[3];
sx q[3];
rz(-1.7700723) q[3];
sx q[3];
rz(2.5602333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-3.0252769) q[2];
rz(2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(11*pi/12) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.6631888) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5696213) q[0];
sx q[0];
rz(-0.53335359) q[0];
sx q[0];
rz(-1.6589952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.19403) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(2.6434968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49839834) q[1];
sx q[1];
rz(-2.8865951) q[1];
sx q[1];
rz(-0.71868371) q[1];
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
rz(1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(-2.647906) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.71815021) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(1.7921599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28045052) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(-0.54235561) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4990436) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(0.16845265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3161105) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(0.18856089) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512882) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.1420508) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
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