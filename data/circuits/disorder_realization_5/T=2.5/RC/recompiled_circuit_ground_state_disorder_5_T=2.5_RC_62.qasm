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
rz(-0.61495334) q[0];
rz(1.3568658) q[1];
sx q[1];
rz(-1.5306127) q[1];
sx q[1];
rz(0.15712486) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8346773) q[0];
sx q[0];
rz(-2.1086725) q[0];
sx q[0];
rz(-2.5468537) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22665568) q[2];
sx q[2];
rz(-1.6427411) q[2];
sx q[2];
rz(-2.8780216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8871044) q[1];
sx q[1];
rz(-0.65782065) q[1];
sx q[1];
rz(-1.5456587) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71977736) q[3];
sx q[3];
rz(-1.3941289) q[3];
sx q[3];
rz(-0.16181419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60363808) q[2];
sx q[2];
rz(-2.8261638) q[2];
sx q[2];
rz(1.8559378) q[2];
rz(1.8042608) q[3];
sx q[3];
rz(-0.95271102) q[3];
sx q[3];
rz(3.0281236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2737357) q[0];
sx q[0];
rz(-1.289176) q[0];
sx q[0];
rz(2.5751172) q[0];
rz(-1.6784338) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(0.64243752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3644407) q[0];
sx q[0];
rz(-0.46555799) q[0];
sx q[0];
rz(-1.4554532) q[0];
x q[1];
rz(-0.64277896) q[2];
sx q[2];
rz(-1.3029429) q[2];
sx q[2];
rz(0.016989022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3088544) q[1];
sx q[1];
rz(-2.2983317) q[1];
sx q[1];
rz(-0.067752167) q[1];
rz(-pi) q[2];
rz(2.1945454) q[3];
sx q[3];
rz(-2.1998341) q[3];
sx q[3];
rz(-1.8333904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2598339) q[2];
sx q[2];
rz(-0.79462785) q[2];
sx q[2];
rz(-0.83414042) q[2];
rz(0.086890876) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(2.0104008) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6777545) q[0];
sx q[0];
rz(-1.0967655) q[0];
sx q[0];
rz(-2.0820397) q[0];
rz(-2.4067267) q[1];
sx q[1];
rz(-0.015345416) q[1];
sx q[1];
rz(0.14618348) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37426078) q[0];
sx q[0];
rz(-1.097299) q[0];
sx q[0];
rz(2.9625508) q[0];
rz(-pi) q[1];
rz(-1.4465815) q[2];
sx q[2];
rz(-1.5669647) q[2];
sx q[2];
rz(1.6506548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9707036) q[1];
sx q[1];
rz(-1.4218937) q[1];
sx q[1];
rz(-1.7335684) q[1];
rz(-0.89208093) q[3];
sx q[3];
rz(-0.35434252) q[3];
sx q[3];
rz(0.49574363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77201468) q[2];
sx q[2];
rz(-2.7507608) q[2];
sx q[2];
rz(0.84103) q[2];
rz(0.05989017) q[3];
sx q[3];
rz(-1.4845029) q[3];
sx q[3];
rz(2.4976775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445207) q[0];
sx q[0];
rz(-2.6332025) q[0];
sx q[0];
rz(-3.0149241) q[0];
rz(-1.4177167) q[1];
sx q[1];
rz(-0.020514943) q[1];
sx q[1];
rz(-2.1518478) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6031109) q[0];
sx q[0];
rz(-0.72906919) q[0];
sx q[0];
rz(3.109594) q[0];
x q[1];
rz(-0.085096882) q[2];
sx q[2];
rz(-1.1750882) q[2];
sx q[2];
rz(2.1425785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7104678) q[1];
sx q[1];
rz(-1.3763274) q[1];
sx q[1];
rz(1.1072913) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2408086) q[3];
sx q[3];
rz(-1.4210555) q[3];
sx q[3];
rz(0.083268585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4385472) q[2];
sx q[2];
rz(-2.788471) q[2];
sx q[2];
rz(0.14567308) q[2];
rz(0.4438256) q[3];
sx q[3];
rz(-1.5668818) q[3];
sx q[3];
rz(1.1183967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038427453) q[0];
sx q[0];
rz(-1.9157836) q[0];
sx q[0];
rz(0.78381729) q[0];
rz(-1.8812284) q[1];
sx q[1];
rz(-0.0030219373) q[1];
sx q[1];
rz(-2.9154215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5617666) q[0];
sx q[0];
rz(-1.34403) q[0];
sx q[0];
rz(-1.0933601) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30136828) q[2];
sx q[2];
rz(-1.971481) q[2];
sx q[2];
rz(-3.1189975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0020272) q[1];
sx q[1];
rz(-2.0420549) q[1];
sx q[1];
rz(-0.084080701) q[1];
rz(-pi) q[2];
rz(-2.2461579) q[3];
sx q[3];
rz(-2.2339719) q[3];
sx q[3];
rz(-1.1010045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1251395) q[2];
sx q[2];
rz(-0.21802248) q[2];
sx q[2];
rz(1.7832635) q[2];
rz(-2.6839117) q[3];
sx q[3];
rz(-1.4384559) q[3];
sx q[3];
rz(2.5911736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0871171) q[0];
sx q[0];
rz(-0.0019019141) q[0];
sx q[0];
rz(0.052363366) q[0];
rz(-2.8730483) q[1];
sx q[1];
rz(-1.152055) q[1];
sx q[1];
rz(-0.23121887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0975039) q[0];
sx q[0];
rz(-1.8262926) q[0];
sx q[0];
rz(1.7520422) q[0];
rz(-pi) q[1];
rz(2.3965712) q[2];
sx q[2];
rz(-1.797849) q[2];
sx q[2];
rz(-0.46316564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0760804) q[1];
sx q[1];
rz(-1.8022746) q[1];
sx q[1];
rz(0.3409909) q[1];
x q[2];
rz(2.1778132) q[3];
sx q[3];
rz(-1.7208817) q[3];
sx q[3];
rz(3.1021189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3117567) q[2];
sx q[2];
rz(-0.51218963) q[2];
sx q[2];
rz(-0.83489418) q[2];
rz(3.0541776) q[3];
sx q[3];
rz(-1.089047) q[3];
sx q[3];
rz(2.1277229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8304623) q[0];
sx q[0];
rz(-2.1489547) q[0];
sx q[0];
rz(-2.2845238) q[0];
rz(2.3509707) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(1.0766006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6031736) q[0];
sx q[0];
rz(-2.7336774) q[0];
sx q[0];
rz(-0.84533219) q[0];
rz(-pi) q[1];
rz(-2.1250379) q[2];
sx q[2];
rz(-0.86807251) q[2];
sx q[2];
rz(-1.6464485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5416661) q[1];
sx q[1];
rz(-1.2866719) q[1];
sx q[1];
rz(-0.78943531) q[1];
x q[2];
rz(2.7995501) q[3];
sx q[3];
rz(-1.6989467) q[3];
sx q[3];
rz(2.8154833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2662346) q[2];
sx q[2];
rz(-2.715832) q[2];
sx q[2];
rz(2.2597964) q[2];
rz(-1.4074696) q[3];
sx q[3];
rz(-0.93972814) q[3];
sx q[3];
rz(0.56434977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.28500104) q[0];
sx q[0];
rz(-0.0017310062) q[0];
sx q[0];
rz(-2.8518387) q[0];
rz(1.2334791) q[1];
sx q[1];
rz(-1.2954442) q[1];
sx q[1];
rz(2.8626056) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.679259) q[0];
sx q[0];
rz(-1.4024319) q[0];
sx q[0];
rz(0.025547693) q[0];
rz(-0.40172949) q[2];
sx q[2];
rz(-1.8332943) q[2];
sx q[2];
rz(-1.0523192) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.853914) q[1];
sx q[1];
rz(-2.8077237) q[1];
sx q[1];
rz(-1.3592657) q[1];
rz(-pi) q[2];
rz(2.7110149) q[3];
sx q[3];
rz(-0.89814808) q[3];
sx q[3];
rz(-2.3585542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6799307) q[2];
sx q[2];
rz(-1.8671904) q[2];
sx q[2];
rz(-2.8604841) q[2];
rz(0.57349652) q[3];
sx q[3];
rz(-0.95279136) q[3];
sx q[3];
rz(2.7219685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24116521) q[0];
sx q[0];
rz(-2.219491) q[0];
sx q[0];
rz(-1.7036555) q[0];
rz(-0.18962139) q[1];
sx q[1];
rz(-0.01556839) q[1];
sx q[1];
rz(1.2356893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3251117) q[0];
sx q[0];
rz(-0.50704038) q[0];
sx q[0];
rz(1.1757686) q[0];
x q[1];
rz(2.3896994) q[2];
sx q[2];
rz(-2.1778669) q[2];
sx q[2];
rz(2.1715669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7687364) q[1];
sx q[1];
rz(-0.33987486) q[1];
sx q[1];
rz(2.0716785) q[1];
x q[2];
rz(-0.37008567) q[3];
sx q[3];
rz(-3.1227263) q[3];
sx q[3];
rz(-2.3436598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8427061) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(-2.5617808) q[2];
rz(1.3696356) q[3];
sx q[3];
rz(-0.42391351) q[3];
sx q[3];
rz(-0.4637318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6746826) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(1.4379733) q[0];
rz(1.9101608) q[1];
sx q[1];
rz(-2.2302901) q[1];
sx q[1];
rz(1.4968754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.623965) q[0];
sx q[0];
rz(-1.8359487) q[0];
sx q[0];
rz(2.318906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52224285) q[2];
sx q[2];
rz(-1.9288315) q[2];
sx q[2];
rz(-0.71292646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.673996) q[1];
sx q[1];
rz(-1.5679545) q[1];
sx q[1];
rz(0.01771042) q[1];
rz(2.8411658) q[3];
sx q[3];
rz(-2.0355371) q[3];
sx q[3];
rz(2.0162051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9897495) q[2];
sx q[2];
rz(-1.5630378) q[2];
sx q[2];
rz(-2.6452276) q[2];
rz(-2.1967891) q[3];
sx q[3];
rz(-0.0050408575) q[3];
sx q[3];
rz(0.70281023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.7127767) q[0];
sx q[0];
rz(2.0145804) q[0];
rz(-1.5827178) q[1];
sx q[1];
rz(-0.3180779) q[1];
sx q[1];
rz(0.19837468) q[1];
rz(1.6678098) q[2];
sx q[2];
rz(-2.8837268) q[2];
sx q[2];
rz(-2.9261598) q[2];
rz(-1.9357135) q[3];
sx q[3];
rz(-2.09874) q[3];
sx q[3];
rz(-2.7345777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
