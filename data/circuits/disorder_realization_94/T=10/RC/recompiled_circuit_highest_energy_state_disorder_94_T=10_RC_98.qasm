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
rz(-1.8985768) q[0];
sx q[0];
rz(-2.0199116) q[0];
sx q[0];
rz(-0.46410528) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(-0.90944374) q[1];
sx q[1];
rz(1.8250725) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5342193) q[0];
sx q[0];
rz(-1.6053622) q[0];
sx q[0];
rz(-0.065163356) q[0];
rz(-pi) q[1];
rz(1.6208956) q[2];
sx q[2];
rz(-1.230579) q[2];
sx q[2];
rz(0.17776793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1471828) q[1];
sx q[1];
rz(-1.3003682) q[1];
sx q[1];
rz(2.5308591) q[1];
rz(1.8620231) q[3];
sx q[3];
rz(-1.3253477) q[3];
sx q[3];
rz(-3.0994741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2509649) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(-0.29699057) q[2];
rz(-2.4885528) q[3];
sx q[3];
rz(-1.5831455) q[3];
sx q[3];
rz(-1.0562586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5121269) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(-0.28489354) q[0];
rz(-0.051305436) q[1];
sx q[1];
rz(-0.76301328) q[1];
sx q[1];
rz(0.6368534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2158686) q[0];
sx q[0];
rz(-2.0179739) q[0];
sx q[0];
rz(-0.53430064) q[0];
rz(2.7899488) q[2];
sx q[2];
rz(-1.571901) q[2];
sx q[2];
rz(-1.0467093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11733774) q[1];
sx q[1];
rz(-1.7452612) q[1];
sx q[1];
rz(-1.4937449) q[1];
x q[2];
rz(0.70592441) q[3];
sx q[3];
rz(-1.4691938) q[3];
sx q[3];
rz(-1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18179831) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(0.96724969) q[2];
rz(-2.8773384) q[3];
sx q[3];
rz(-1.5033009) q[3];
sx q[3];
rz(1.2359469) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-1.0775403) q[1];
sx q[1];
rz(-2.8376104) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39897746) q[0];
sx q[0];
rz(-1.0053867) q[0];
sx q[0];
rz(3.1088367) q[0];
rz(-pi) q[1];
rz(2.9037232) q[2];
sx q[2];
rz(-1.6697236) q[2];
sx q[2];
rz(1.400927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5604104) q[1];
sx q[1];
rz(-0.85542233) q[1];
sx q[1];
rz(0.13611273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68328698) q[3];
sx q[3];
rz(-1.7257455) q[3];
sx q[3];
rz(2.4187627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1659282) q[2];
sx q[2];
rz(-1.0378446) q[2];
sx q[2];
rz(-0.69592875) q[2];
rz(2.0078697) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(2.5886562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(-1.0473921) q[0];
rz(-1.1737431) q[1];
sx q[1];
rz(-2.4672697) q[1];
sx q[1];
rz(1.6204999) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0516617) q[0];
sx q[0];
rz(-1.0703846) q[0];
sx q[0];
rz(-1.7080281) q[0];
x q[1];
rz(-1.1999454) q[2];
sx q[2];
rz(-0.74127889) q[2];
sx q[2];
rz(1.0560869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4075267) q[1];
sx q[1];
rz(-0.95469771) q[1];
sx q[1];
rz(1.0259088) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32471809) q[3];
sx q[3];
rz(-1.3007264) q[3];
sx q[3];
rz(2.3809759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(-2.6037237) q[2];
rz(-1.6543903) q[3];
sx q[3];
rz(-0.40707773) q[3];
sx q[3];
rz(0.68638221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.2122413) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(-0.60337639) q[0];
rz(-0.76250184) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(-0.71294436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0423316) q[0];
sx q[0];
rz(-1.5047899) q[0];
sx q[0];
rz(1.1724654) q[0];
rz(-3.001956) q[2];
sx q[2];
rz(-2.0296445) q[2];
sx q[2];
rz(-0.81880867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11345574) q[1];
sx q[1];
rz(-2.0638247) q[1];
sx q[1];
rz(-2.9131469) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98865786) q[3];
sx q[3];
rz(-1.0352367) q[3];
sx q[3];
rz(-2.5861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1684299) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(-1.6707576) q[2];
rz(-2.6245608) q[3];
sx q[3];
rz(-2.0972589) q[3];
sx q[3];
rz(-1.2928591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9789199) q[0];
sx q[0];
rz(-2.6277268) q[0];
sx q[0];
rz(-0.55147076) q[0];
rz(-3.1061213) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(1.6112526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1288958) q[0];
sx q[0];
rz(-0.91980359) q[0];
sx q[0];
rz(-0.76437073) q[0];
rz(-pi) q[1];
rz(-1.6257203) q[2];
sx q[2];
rz(-2.3936205) q[2];
sx q[2];
rz(0.8665646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9236218) q[1];
sx q[1];
rz(-1.21232) q[1];
sx q[1];
rz(2.3450407) q[1];
rz(-pi) q[2];
rz(0.084613581) q[3];
sx q[3];
rz(-2.222098) q[3];
sx q[3];
rz(-0.29825975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0643206) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.2426097) q[2];
rz(-2.6284435) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(-1.3750403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.649491) q[0];
sx q[0];
rz(-2.212337) q[0];
sx q[0];
rz(0.11665601) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(1.6366417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0610675) q[0];
sx q[0];
rz(-1.8634029) q[0];
sx q[0];
rz(-0.82340296) q[0];
x q[1];
rz(1.9595615) q[2];
sx q[2];
rz(-1.6259527) q[2];
sx q[2];
rz(-2.9832632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.44918647) q[1];
sx q[1];
rz(-0.80185651) q[1];
sx q[1];
rz(2.0795601) q[1];
x q[2];
rz(-2.4890763) q[3];
sx q[3];
rz(-1.8015141) q[3];
sx q[3];
rz(1.1200865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.92631212) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(0.8872633) q[2];
rz(0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(-2.5086856) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55650869) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(-2.4853117) q[0];
rz(-2.9346924) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042699634) q[0];
sx q[0];
rz(-1.6535361) q[0];
sx q[0];
rz(1.6059884) q[0];
rz(2.0711259) q[2];
sx q[2];
rz(-2.4917951) q[2];
sx q[2];
rz(-1.4184679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43133538) q[1];
sx q[1];
rz(-1.4568304) q[1];
sx q[1];
rz(2.9821787) q[1];
rz(-pi) q[2];
rz(-2.5303305) q[3];
sx q[3];
rz(-1.3484133) q[3];
sx q[3];
rz(0.56764795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8021585) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(-2.459724) q[2];
rz(1.9631466) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(-1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449529) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(-0.21016453) q[0];
rz(-2.919803) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(-0.04235696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9440749) q[0];
sx q[0];
rz(-2.952708) q[0];
sx q[0];
rz(-0.11795363) q[0];
x q[1];
rz(-2.1788939) q[2];
sx q[2];
rz(-2.1129932) q[2];
sx q[2];
rz(1.0350943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90396229) q[1];
sx q[1];
rz(-2.1231066) q[1];
sx q[1];
rz(-1.2321074) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5333535) q[3];
sx q[3];
rz(-1.9755529) q[3];
sx q[3];
rz(-0.08431708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8611341) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(0.66656485) q[2];
rz(-2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53369451) q[0];
sx q[0];
rz(-2.7602637) q[0];
sx q[0];
rz(0.68699849) q[0];
rz(-1.8580565) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(-0.57354617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6015778) q[0];
sx q[0];
rz(-0.83160831) q[0];
sx q[0];
rz(2.4180146) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9098689) q[2];
sx q[2];
rz(-1.1292233) q[2];
sx q[2];
rz(-1.9927466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.191205) q[1];
sx q[1];
rz(-0.6362241) q[1];
sx q[1];
rz(-2.1649182) q[1];
rz(-0.73856662) q[3];
sx q[3];
rz(-0.99792483) q[3];
sx q[3];
rz(-2.3923739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0987229) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(-0.086611835) q[2];
rz(-3.0232271) q[3];
sx q[3];
rz(-1.5425073) q[3];
sx q[3];
rz(3.0983483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31556986) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(0.28620537) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(-2.8092842) q[2];
sx q[2];
rz(-1.341991) q[2];
sx q[2];
rz(1.1534635) q[2];
rz(1.2729011) q[3];
sx q[3];
rz(-2.0023228) q[3];
sx q[3];
rz(1.7154233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
