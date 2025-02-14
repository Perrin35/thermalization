OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.759999) q[0];
sx q[0];
rz(-0.17188369) q[0];
sx q[0];
rz(2.5266393) q[0];
rz(1.3568658) q[1];
sx q[1];
rz(-1.5306127) q[1];
sx q[1];
rz(0.15712486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2113094) q[0];
sx q[0];
rz(-2.0727949) q[0];
sx q[0];
rz(-2.1949578) q[0];
rz(-1.6446227) q[2];
sx q[2];
rz(-1.7968555) q[2];
sx q[2];
rz(1.3238012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8871044) q[1];
sx q[1];
rz(-0.65782065) q[1];
sx q[1];
rz(1.595934) q[1];
rz(-2.4218153) q[3];
sx q[3];
rz(-1.3941289) q[3];
sx q[3];
rz(-0.16181419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60363808) q[2];
sx q[2];
rz(-0.31542888) q[2];
sx q[2];
rz(-1.8559378) q[2];
rz(1.8042608) q[3];
sx q[3];
rz(-2.1888816) q[3];
sx q[3];
rz(-3.0281236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.86785698) q[0];
sx q[0];
rz(-1.8524167) q[0];
sx q[0];
rz(-0.56647545) q[0];
rz(1.6784338) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(2.4991551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89680371) q[0];
sx q[0];
rz(-1.519108) q[0];
sx q[0];
rz(-2.0336853) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9011097) q[2];
sx q[2];
rz(-0.95446247) q[2];
sx q[2];
rz(1.7834477) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21684619) q[1];
sx q[1];
rz(-1.5202151) q[1];
sx q[1];
rz(-0.84212007) q[1];
x q[2];
rz(2.1945454) q[3];
sx q[3];
rz(-2.1998341) q[3];
sx q[3];
rz(1.3082023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88175875) q[2];
sx q[2];
rz(-0.79462785) q[2];
sx q[2];
rz(0.83414042) q[2];
rz(-3.0547018) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(2.0104008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46383816) q[0];
sx q[0];
rz(-1.0967655) q[0];
sx q[0];
rz(1.0595529) q[0];
rz(0.73486596) q[1];
sx q[1];
rz(-0.015345416) q[1];
sx q[1];
rz(-2.9954092) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8627166) q[0];
sx q[0];
rz(-1.4116316) q[0];
sx q[0];
rz(2.0508641) q[0];
rz(1.6950111) q[2];
sx q[2];
rz(-1.5746279) q[2];
sx q[2];
rz(1.4909378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66524341) q[1];
sx q[1];
rz(-0.2201597) q[1];
sx q[1];
rz(2.317628) q[1];
rz(-2.2495117) q[3];
sx q[3];
rz(-2.7872501) q[3];
sx q[3];
rz(-2.645849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77201468) q[2];
sx q[2];
rz(-0.39083189) q[2];
sx q[2];
rz(-2.3005627) q[2];
rz(-3.0817025) q[3];
sx q[3];
rz(-1.6570897) q[3];
sx q[3];
rz(-2.4976775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.197072) q[0];
sx q[0];
rz(-2.6332025) q[0];
sx q[0];
rz(-3.0149241) q[0];
rz(1.4177167) q[1];
sx q[1];
rz(-3.1210777) q[1];
sx q[1];
rz(-2.1518478) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0084466) q[0];
sx q[0];
rz(-1.5494816) q[0];
sx q[0];
rz(-0.7288148) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1737972) q[2];
sx q[2];
rz(-1.6493031) q[2];
sx q[2];
rz(2.5369413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77074487) q[1];
sx q[1];
rz(-2.641692) q[1];
sx q[1];
rz(1.1558394) q[1];
rz(-pi) q[2];
rz(0.56405385) q[3];
sx q[3];
rz(-2.8587935) q[3];
sx q[3];
rz(-2.0334311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7030455) q[2];
sx q[2];
rz(-0.35312167) q[2];
sx q[2];
rz(-0.14567308) q[2];
rz(0.4438256) q[3];
sx q[3];
rz(-1.5747109) q[3];
sx q[3];
rz(2.023196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038427453) q[0];
sx q[0];
rz(-1.225809) q[0];
sx q[0];
rz(-0.78381729) q[0];
rz(1.2603643) q[1];
sx q[1];
rz(-0.0030219373) q[1];
sx q[1];
rz(0.22617117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.579826) q[0];
sx q[0];
rz(-1.34403) q[0];
sx q[0];
rz(1.0933601) q[0];
rz(2.1820081) q[2];
sx q[2];
rz(-0.49641616) q[2];
sx q[2];
rz(0.69545262) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1395654) q[1];
sx q[1];
rz(-2.0420549) q[1];
sx q[1];
rz(0.084080701) q[1];
rz(-pi) q[2];
x q[2];
rz(2.466684) q[3];
sx q[3];
rz(-2.2332472) q[3];
sx q[3];
rz(-1.1250594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0164531) q[2];
sx q[2];
rz(-0.21802248) q[2];
sx q[2];
rz(-1.3583292) q[2];
rz(2.6839117) q[3];
sx q[3];
rz(-1.7031368) q[3];
sx q[3];
rz(-0.55041909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0871171) q[0];
sx q[0];
rz(-0.0019019141) q[0];
sx q[0];
rz(-3.0892293) q[0];
rz(-0.26854435) q[1];
sx q[1];
rz(-1.9895376) q[1];
sx q[1];
rz(2.9103738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044088783) q[0];
sx q[0];
rz(-1.8262926) q[0];
sx q[0];
rz(-1.3895504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8753239) q[2];
sx q[2];
rz(-0.84919676) q[2];
sx q[2];
rz(0.90291427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0760804) q[1];
sx q[1];
rz(-1.8022746) q[1];
sx q[1];
rz(2.8006018) q[1];
rz(-2.1778132) q[3];
sx q[3];
rz(-1.7208817) q[3];
sx q[3];
rz(-3.1021189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3117567) q[2];
sx q[2];
rz(-2.629403) q[2];
sx q[2];
rz(-2.3066985) q[2];
rz(0.087415047) q[3];
sx q[3];
rz(-1.089047) q[3];
sx q[3];
rz(-2.1277229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.3111303) q[0];
sx q[0];
rz(-0.99263793) q[0];
sx q[0];
rz(-2.2845238) q[0];
rz(-0.79062194) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(-2.0649921) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6031736) q[0];
sx q[0];
rz(-0.40791528) q[0];
sx q[0];
rz(0.84533219) q[0];
rz(-2.1250379) q[2];
sx q[2];
rz(-2.2735201) q[2];
sx q[2];
rz(1.6464485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5416661) q[1];
sx q[1];
rz(-1.2866719) q[1];
sx q[1];
rz(-2.3521573) q[1];
rz(1.7067327) q[3];
sx q[3];
rz(-1.9099209) q[3];
sx q[3];
rz(1.1992169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2662346) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(0.88179624) q[2];
rz(1.734123) q[3];
sx q[3];
rz(-0.93972814) q[3];
sx q[3];
rz(0.56434977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8565916) q[0];
sx q[0];
rz(-0.0017310062) q[0];
sx q[0];
rz(-2.8518387) q[0];
rz(-1.9081135) q[1];
sx q[1];
rz(-1.2954442) q[1];
sx q[1];
rz(-0.27898702) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10418071) q[0];
sx q[0];
rz(-1.54561) q[0];
sx q[0];
rz(-1.7392147) q[0];
rz(-1.2867497) q[2];
sx q[2];
rz(-1.1835754) q[2];
sx q[2];
rz(-0.40867886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.853914) q[1];
sx q[1];
rz(-2.8077237) q[1];
sx q[1];
rz(-1.7823269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2905228) q[3];
sx q[3];
rz(-1.2382231) q[3];
sx q[3];
rz(0.50902396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6799307) q[2];
sx q[2];
rz(-1.2744023) q[2];
sx q[2];
rz(0.28110853) q[2];
rz(-0.57349652) q[3];
sx q[3];
rz(-0.95279136) q[3];
sx q[3];
rz(0.41962418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.24116521) q[0];
sx q[0];
rz(-2.219491) q[0];
sx q[0];
rz(-1.7036555) q[0];
rz(0.18962139) q[1];
sx q[1];
rz(-0.01556839) q[1];
sx q[1];
rz(1.9059034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8164809) q[0];
sx q[0];
rz(-0.50704038) q[0];
sx q[0];
rz(-1.1757686) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3310555) q[2];
sx q[2];
rz(-0.97522465) q[2];
sx q[2];
rz(-3.0308804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37285629) q[1];
sx q[1];
rz(-2.8017178) q[1];
sx q[1];
rz(2.0716785) q[1];
x q[2];
rz(-1.5639717) q[3];
sx q[3];
rz(-1.5532074) q[3];
sx q[3];
rz(2.7138055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8427061) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(2.5617808) q[2];
rz(-1.3696356) q[3];
sx q[3];
rz(-2.7176791) q[3];
sx q[3];
rz(2.6778609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.46691) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(1.7036194) q[0];
rz(1.2314318) q[1];
sx q[1];
rz(-0.91130251) q[1];
sx q[1];
rz(1.4968754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8498488) q[0];
sx q[0];
rz(-0.8546517) q[0];
sx q[0];
rz(-0.35479389) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6193498) q[2];
sx q[2];
rz(-1.2127611) q[2];
sx q[2];
rz(0.71292646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0383427) q[1];
sx q[1];
rz(-1.553086) q[1];
sx q[1];
rz(-1.5679541) q[1];
rz(-pi) q[2];
rz(1.0874463) q[3];
sx q[3];
rz(-1.8385244) q[3];
sx q[3];
rz(-0.30743956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15184312) q[2];
sx q[2];
rz(-1.5785549) q[2];
sx q[2];
rz(2.6452276) q[2];
rz(2.1967891) q[3];
sx q[3];
rz(-0.0050408575) q[3];
sx q[3];
rz(2.4387824) q[3];
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
rz(-1.6481358) q[0];
sx q[0];
rz(-1.428816) q[0];
sx q[0];
rz(-1.1270123) q[0];
rz(-1.5827178) q[1];
sx q[1];
rz(-0.3180779) q[1];
sx q[1];
rz(0.19837468) q[1];
rz(-1.4737829) q[2];
sx q[2];
rz(-2.8837268) q[2];
sx q[2];
rz(-2.9261598) q[2];
rz(0.55806969) q[3];
sx q[3];
rz(-1.8841828) q[3];
sx q[3];
rz(1.7876865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
