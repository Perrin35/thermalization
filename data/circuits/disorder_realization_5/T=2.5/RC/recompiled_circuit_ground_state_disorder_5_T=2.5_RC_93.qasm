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
rz(2.969709) q[0];
sx q[0];
rz(10.039731) q[0];
rz(1.3568658) q[1];
sx q[1];
rz(4.7525726) q[1];
sx q[1];
rz(9.5819028) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93028322) q[0];
sx q[0];
rz(-1.0687978) q[0];
sx q[0];
rz(0.94663488) q[0];
x q[1];
rz(0.31034361) q[2];
sx q[2];
rz(-0.23761286) q[2];
sx q[2];
rz(-1.0050736) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3362017) q[1];
sx q[1];
rz(-1.5861643) q[1];
sx q[1];
rz(0.91312855) q[1];
x q[2];
rz(-2.8771175) q[3];
sx q[3];
rz(-2.4042359) q[3];
sx q[3];
rz(1.5347598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5379546) q[2];
sx q[2];
rz(-2.8261638) q[2];
sx q[2];
rz(1.2856548) q[2];
rz(-1.8042608) q[3];
sx q[3];
rz(-0.95271102) q[3];
sx q[3];
rz(0.11346909) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86785698) q[0];
sx q[0];
rz(-1.289176) q[0];
sx q[0];
rz(-2.5751172) q[0];
rz(-1.4631588) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(2.4991551) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77715194) q[0];
sx q[0];
rz(-0.46555799) q[0];
sx q[0];
rz(-1.4554532) q[0];
rz(-pi) q[1];
rz(-2.4988137) q[2];
sx q[2];
rz(-1.3029429) q[2];
sx q[2];
rz(3.1246036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21684619) q[1];
sx q[1];
rz(-1.6213776) q[1];
sx q[1];
rz(-2.2994726) q[1];
rz(-pi) q[2];
rz(-2.1945454) q[3];
sx q[3];
rz(-2.1998341) q[3];
sx q[3];
rz(1.8333904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88175875) q[2];
sx q[2];
rz(-2.3469648) q[2];
sx q[2];
rz(-2.3074522) q[2];
rz(-0.086890876) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(-2.0104008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6777545) q[0];
sx q[0];
rz(-1.0967655) q[0];
sx q[0];
rz(1.0595529) q[0];
rz(-0.73486596) q[1];
sx q[1];
rz(-0.015345416) q[1];
sx q[1];
rz(2.9954092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8627166) q[0];
sx q[0];
rz(-1.729961) q[0];
sx q[0];
rz(-1.0907286) q[0];
x q[1];
rz(1.4465815) q[2];
sx q[2];
rz(-1.5746279) q[2];
sx q[2];
rz(-1.4909378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7173269) q[1];
sx q[1];
rz(-1.4098412) q[1];
sx q[1];
rz(-2.9907254) q[1];
x q[2];
rz(1.2904148) q[3];
sx q[3];
rz(-1.7903847) q[3];
sx q[3];
rz(-1.7226579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77201468) q[2];
sx q[2];
rz(-0.39083189) q[2];
sx q[2];
rz(-2.3005627) q[2];
rz(0.05989017) q[3];
sx q[3];
rz(-1.6570897) q[3];
sx q[3];
rz(-2.4976775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445207) q[0];
sx q[0];
rz(-2.6332025) q[0];
sx q[0];
rz(3.0149241) q[0];
rz(1.7238759) q[1];
sx q[1];
rz(-0.020514943) q[1];
sx q[1];
rz(-2.1518478) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58137547) q[0];
sx q[0];
rz(-0.84218431) q[0];
sx q[0];
rz(-1.5993662) q[0];
rz(-pi) q[1];
x q[1];
rz(0.085096882) q[2];
sx q[2];
rz(-1.1750882) q[2];
sx q[2];
rz(-2.1425785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2359611) q[1];
sx q[1];
rz(-1.1166908) q[1];
sx q[1];
rz(-0.21673165) q[1];
rz(-1.7249171) q[3];
sx q[3];
rz(-1.8088578) q[3];
sx q[3];
rz(-1.4509089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1031652) q[0];
sx q[0];
rz(-1.9157836) q[0];
sx q[0];
rz(-0.78381729) q[0];
rz(-1.8812284) q[1];
sx q[1];
rz(-3.1385707) q[1];
sx q[1];
rz(-0.22617117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5617666) q[0];
sx q[0];
rz(-1.7975627) q[0];
sx q[0];
rz(1.0933601) q[0];
rz(-0.95958454) q[2];
sx q[2];
rz(-0.49641616) q[2];
sx q[2];
rz(-2.44614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0020272) q[1];
sx q[1];
rz(-1.0995378) q[1];
sx q[1];
rz(-3.057512) q[1];
rz(-pi) q[2];
rz(-0.67490863) q[3];
sx q[3];
rz(-0.90834544) q[3];
sx q[3];
rz(1.1250594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1251395) q[2];
sx q[2];
rz(-2.9235702) q[2];
sx q[2];
rz(-1.7832635) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054475527) q[0];
sx q[0];
rz(-3.1396907) q[0];
sx q[0];
rz(-0.052363366) q[0];
rz(2.8730483) q[1];
sx q[1];
rz(-1.9895376) q[1];
sx q[1];
rz(-0.23121887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58329952) q[0];
sx q[0];
rz(-0.31210408) q[0];
sx q[0];
rz(2.5375478) q[0];
rz(2.3965712) q[2];
sx q[2];
rz(-1.797849) q[2];
sx q[2];
rz(-0.46316564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42405427) q[1];
sx q[1];
rz(-1.9023384) q[1];
sx q[1];
rz(1.3257205) q[1];
rz(-0.96377947) q[3];
sx q[3];
rz(-1.4207109) q[3];
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
rz(-0.51218963) q[2];
sx q[2];
rz(2.3066985) q[2];
rz(-3.0541776) q[3];
sx q[3];
rz(-2.0525457) q[3];
sx q[3];
rz(2.1277229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8304623) q[0];
sx q[0];
rz(-0.99263793) q[0];
sx q[0];
rz(0.8570689) q[0];
rz(-0.79062194) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(1.0766006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83504377) q[0];
sx q[0];
rz(-1.2694512) q[0];
sx q[0];
rz(0.27923601) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55600682) q[2];
sx q[2];
rz(-0.86472325) q[2];
sx q[2];
rz(0.88265693) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4461313) q[1];
sx q[1];
rz(-2.3205955) q[1];
sx q[1];
rz(-1.9638792) q[1];
rz(-pi) q[2];
rz(2.7748037) q[3];
sx q[3];
rz(-2.7772172) q[3];
sx q[3];
rz(-1.5522788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2662346) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(-2.2597964) q[2];
rz(-1.4074696) q[3];
sx q[3];
rz(-2.2018645) q[3];
sx q[3];
rz(-0.56434977) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8565916) q[0];
sx q[0];
rz(-0.0017310062) q[0];
sx q[0];
rz(-0.28975394) q[0];
rz(1.9081135) q[1];
sx q[1];
rz(-1.2954442) q[1];
sx q[1];
rz(0.27898702) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.679259) q[0];
sx q[0];
rz(-1.4024319) q[0];
sx q[0];
rz(0.025547693) q[0];
rz(-pi) q[1];
rz(2.7398632) q[2];
sx q[2];
rz(-1.8332943) q[2];
sx q[2];
rz(2.0892734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51117544) q[1];
sx q[1];
rz(-1.8969444) q[1];
sx q[1];
rz(0.072695865) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0881422) q[3];
sx q[3];
rz(-0.78018809) q[3];
sx q[3];
rz(-1.4182981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4616619) q[2];
sx q[2];
rz(-1.2744023) q[2];
sx q[2];
rz(-0.28110853) q[2];
rz(2.5680961) q[3];
sx q[3];
rz(-2.1888013) q[3];
sx q[3];
rz(-0.41962418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.24116521) q[0];
sx q[0];
rz(-0.92210162) q[0];
sx q[0];
rz(1.4379372) q[0];
rz(-0.18962139) q[1];
sx q[1];
rz(-0.01556839) q[1];
sx q[1];
rz(1.2356893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10383725) q[0];
sx q[0];
rz(-1.7587737) q[0];
sx q[0];
rz(2.044552) q[0];
rz(-pi) q[1];
rz(2.3478144) q[2];
sx q[2];
rz(-2.214162) q[2];
sx q[2];
rz(-1.9930594) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6744102) q[1];
sx q[1];
rz(-1.7315718) q[1];
sx q[1];
rz(-1.2700446) q[1];
x q[2];
rz(-1.5639717) q[3];
sx q[3];
rz(-1.5883852) q[3];
sx q[3];
rz(0.42778711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29888657) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(0.57981181) q[2];
rz(1.3696356) q[3];
sx q[3];
rz(-2.7176791) q[3];
sx q[3];
rz(-2.6778609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.46691) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(1.7036194) q[0];
rz(-1.9101608) q[1];
sx q[1];
rz(-0.91130251) q[1];
sx q[1];
rz(1.4968754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2220331) q[0];
sx q[0];
rz(-0.78503937) q[0];
sx q[0];
rz(1.9505984) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6193498) q[2];
sx q[2];
rz(-1.9288315) q[2];
sx q[2];
rz(2.4286662) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4675967) q[1];
sx q[1];
rz(-1.5679545) q[1];
sx q[1];
rz(0.01771042) q[1];
rz(2.8411658) q[3];
sx q[3];
rz(-1.1060556) q[3];
sx q[3];
rz(1.1253875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9897495) q[2];
sx q[2];
rz(-1.5785549) q[2];
sx q[2];
rz(-2.6452276) q[2];
rz(-0.94480354) q[3];
sx q[3];
rz(-3.1365518) q[3];
sx q[3];
rz(-2.4387824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.6678098) q[2];
sx q[2];
rz(-0.2578659) q[2];
sx q[2];
rz(0.21543287) q[2];
rz(2.583523) q[3];
sx q[3];
rz(-1.2574099) q[3];
sx q[3];
rz(-1.3539061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
