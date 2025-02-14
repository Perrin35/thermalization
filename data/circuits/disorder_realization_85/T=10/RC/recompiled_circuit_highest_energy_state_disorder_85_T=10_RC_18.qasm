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
rz(0.99211168) q[0];
sx q[0];
rz(3.8881128) q[0];
sx q[0];
rz(9.9915656) q[0];
rz(-1.8911288) q[1];
sx q[1];
rz(-1.1486147) q[1];
sx q[1];
rz(-2.221938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88094372) q[0];
sx q[0];
rz(-2.3365417) q[0];
sx q[0];
rz(-0.22329231) q[0];
rz(-2.0777474) q[2];
sx q[2];
rz(-2.2301444) q[2];
sx q[2];
rz(-0.44753513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30324591) q[1];
sx q[1];
rz(-0.64323264) q[1];
sx q[1];
rz(0.077416181) q[1];
x q[2];
rz(0.53594671) q[3];
sx q[3];
rz(-1.4466373) q[3];
sx q[3];
rz(2.9298327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93996843) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(2.0882108) q[2];
rz(2.4273704) q[3];
sx q[3];
rz(-2.768399) q[3];
sx q[3];
rz(1.2936973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5193704) q[0];
sx q[0];
rz(-0.70835963) q[0];
sx q[0];
rz(-0.57193065) q[0];
rz(-1.8427303) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(0.78972185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1173218) q[0];
sx q[0];
rz(-0.3171176) q[0];
sx q[0];
rz(-1.8834524) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9339173) q[2];
sx q[2];
rz(-2.4768314) q[2];
sx q[2];
rz(1.1663811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8487723) q[1];
sx q[1];
rz(-2.1629984) q[1];
sx q[1];
rz(1.2315537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5511207) q[3];
sx q[3];
rz(-0.76667029) q[3];
sx q[3];
rz(1.1696512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.3624462) q[2];
rz(-0.29081523) q[3];
sx q[3];
rz(-1.3664061) q[3];
sx q[3];
rz(-3.0349558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0027851) q[0];
sx q[0];
rz(-2.0464351) q[0];
sx q[0];
rz(-0.7005257) q[0];
rz(0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-2.9170759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1248847) q[0];
sx q[0];
rz(-1.2205692) q[0];
sx q[0];
rz(1.4136397) q[0];
rz(-pi) q[1];
rz(-1.4204558) q[2];
sx q[2];
rz(-2.5579961) q[2];
sx q[2];
rz(-2.4395669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3923334) q[1];
sx q[1];
rz(-0.9888923) q[1];
sx q[1];
rz(-3.0514293) q[1];
rz(-0.11364524) q[3];
sx q[3];
rz(-0.84535852) q[3];
sx q[3];
rz(0.64021969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38480467) q[2];
sx q[2];
rz(-2.2497358) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(3.0745506) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(-2.4750347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(-3.1224342) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.7117585) q[1];
sx q[1];
rz(0.0044936831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2990103) q[0];
sx q[0];
rz(-1.0753462) q[0];
sx q[0];
rz(0.61758496) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.121663) q[2];
sx q[2];
rz(-2.0259078) q[2];
sx q[2];
rz(-3.1246076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42378615) q[1];
sx q[1];
rz(-2.7423334) q[1];
sx q[1];
rz(2.0171793) q[1];
rz(-pi) q[2];
rz(1.6636623) q[3];
sx q[3];
rz(-1.3680653) q[3];
sx q[3];
rz(-0.40846881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4517453) q[2];
sx q[2];
rz(-2.0514026) q[2];
sx q[2];
rz(-1.4534265) q[2];
rz(1.8870185) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(-2.5793251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8101863) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(1.1619262) q[0];
rz(-0.76238531) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(-2.7120178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3439826) q[0];
sx q[0];
rz(-1.6071885) q[0];
sx q[0];
rz(-1.9790566) q[0];
rz(-1.3061021) q[2];
sx q[2];
rz(-2.3117723) q[2];
sx q[2];
rz(3.004937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1602729) q[1];
sx q[1];
rz(-1.2339051) q[1];
sx q[1];
rz(2.261922) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3394903) q[3];
sx q[3];
rz(-0.94565839) q[3];
sx q[3];
rz(0.021566077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9261711) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(1.8298836) q[2];
rz(2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-0.1596182) q[0];
sx q[0];
rz(1.106369) q[0];
rz(3.0270992) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(0.053650275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4729484) q[0];
sx q[0];
rz(-1.3151752) q[0];
sx q[0];
rz(-3.1356407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6752376) q[2];
sx q[2];
rz(-2.375583) q[2];
sx q[2];
rz(0.05482373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70354706) q[1];
sx q[1];
rz(-1.5168191) q[1];
sx q[1];
rz(-1.0590068) q[1];
rz(0.76308454) q[3];
sx q[3];
rz(-2.3474135) q[3];
sx q[3];
rz(1.1556243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2758808) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-0.50714058) q[2];
rz(-1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5116665) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(-3.1003057) q[0];
rz(-0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(-0.98947492) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7558003) q[0];
sx q[0];
rz(-2.4210701) q[0];
sx q[0];
rz(2.3260443) q[0];
rz(-pi) q[1];
rz(-2.6648952) q[2];
sx q[2];
rz(-0.39604353) q[2];
sx q[2];
rz(0.11644289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7724645) q[1];
sx q[1];
rz(-2.1434692) q[1];
sx q[1];
rz(-0.84267183) q[1];
rz(-0.37240828) q[3];
sx q[3];
rz(-0.95015991) q[3];
sx q[3];
rz(-0.90639508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(0.011693444) q[3];
sx q[3];
rz(-1.6418567) q[3];
sx q[3];
rz(3.0288127) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3149253) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(0.23183204) q[0];
rz(-0.58206093) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(1.9225165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2363289) q[0];
sx q[0];
rz(-1.4538611) q[0];
sx q[0];
rz(0.56044436) q[0];
rz(0.35301669) q[2];
sx q[2];
rz(-1.6200049) q[2];
sx q[2];
rz(2.4491058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92182175) q[1];
sx q[1];
rz(-1.3857462) q[1];
sx q[1];
rz(-2.9875523) q[1];
rz(-2.1120511) q[3];
sx q[3];
rz(-2.2502021) q[3];
sx q[3];
rz(-0.32493362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(0.78682023) q[2];
rz(-1.6400853) q[3];
sx q[3];
rz(-2.7101176) q[3];
sx q[3];
rz(-1.7155581) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96104923) q[0];
sx q[0];
rz(-1.7683872) q[0];
sx q[0];
rz(0.94026646) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(0.14642265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347552) q[0];
sx q[0];
rz(-1.869258) q[0];
sx q[0];
rz(1.7025856) q[0];
rz(-pi) q[1];
rz(-0.27250473) q[2];
sx q[2];
rz(-1.8159869) q[2];
sx q[2];
rz(1.6134855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89288974) q[1];
sx q[1];
rz(-1.3254524) q[1];
sx q[1];
rz(2.4470083) q[1];
rz(-pi) q[2];
rz(-3.0412263) q[3];
sx q[3];
rz(-2.3477738) q[3];
sx q[3];
rz(-2.9507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.08708295) q[2];
sx q[2];
rz(-1.6528218) q[2];
sx q[2];
rz(0.83756891) q[2];
rz(-2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(-0.78105175) q[3];
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
rz(0.26302108) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(-1.7735057) q[0];
rz(-3.0655762) q[1];
sx q[1];
rz(-0.58034211) q[1];
sx q[1];
rz(1.1200294) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.66066) q[0];
sx q[0];
rz(-1.2613861) q[0];
sx q[0];
rz(-0.27746986) q[0];
rz(-0.35181184) q[2];
sx q[2];
rz(-2.4016671) q[2];
sx q[2];
rz(-2.6463395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.315584) q[1];
sx q[1];
rz(-0.30696973) q[1];
sx q[1];
rz(-1.2493285) q[1];
x q[2];
rz(0.10346966) q[3];
sx q[3];
rz(-1.7333687) q[3];
sx q[3];
rz(1.9601456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40176216) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(-2.8945727) q[2];
rz(2.5714827) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(-2.0232239) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0982672) q[0];
sx q[0];
rz(-1.387431) q[0];
sx q[0];
rz(-0.36547216) q[0];
rz(3.0430766) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-1.6016207) q[2];
sx q[2];
rz(-0.41308232) q[2];
sx q[2];
rz(-0.272258) q[2];
rz(-2.5539342) q[3];
sx q[3];
rz(-2.6937204) q[3];
sx q[3];
rz(3.1093521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
