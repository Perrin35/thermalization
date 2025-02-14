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
rz(-0.22696683) q[0];
sx q[0];
rz(-2.5300955) q[0];
sx q[0];
rz(-0.55846941) q[0];
rz(-2.5418169) q[1];
sx q[1];
rz(-1.37473) q[1];
sx q[1];
rz(3.0629646) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087705165) q[0];
sx q[0];
rz(-1.0358521) q[0];
sx q[0];
rz(0.91607262) q[0];
rz(-pi) q[1];
rz(-1.6958649) q[2];
sx q[2];
rz(-2.2457221) q[2];
sx q[2];
rz(-1.1488016) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6650162) q[1];
sx q[1];
rz(-1.652583) q[1];
sx q[1];
rz(0.25340124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8923678) q[3];
sx q[3];
rz(-1.9059637) q[3];
sx q[3];
rz(0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2416396) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-1.0345577) q[2];
rz(-3.0311846) q[3];
sx q[3];
rz(-1.6610828) q[3];
sx q[3];
rz(-0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13650525) q[0];
sx q[0];
rz(-1.6965447) q[0];
sx q[0];
rz(1.1862296) q[0];
rz(1.8738481) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.3892106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1339552) q[0];
sx q[0];
rz(-1.8764155) q[0];
sx q[0];
rz(-2.1307178) q[0];
rz(1.6824426) q[2];
sx q[2];
rz(-1.0918587) q[2];
sx q[2];
rz(-3.0618596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6450883) q[1];
sx q[1];
rz(-1.8594655) q[1];
sx q[1];
rz(1.9157857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4382382) q[3];
sx q[3];
rz(-2.5131559) q[3];
sx q[3];
rz(-2.4579163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(-1.7051075) q[2];
rz(-1.2596333) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(-3.1148541) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36419511) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.2170323) q[0];
rz(2.9265535) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(-1.4395813) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3830945) q[0];
sx q[0];
rz(-1.6949843) q[0];
sx q[0];
rz(-0.66207768) q[0];
rz(-pi) q[1];
rz(-1.0595752) q[2];
sx q[2];
rz(-1.4892092) q[2];
sx q[2];
rz(2.4961584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7172723) q[1];
sx q[1];
rz(-1.3417146) q[1];
sx q[1];
rz(-3.0113264) q[1];
x q[2];
rz(-2.2398364) q[3];
sx q[3];
rz(-2.7983411) q[3];
sx q[3];
rz(-0.95860976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3415459) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(2.3254507) q[2];
rz(-0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(-1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4207882) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(-2.8724443) q[0];
rz(1.3828269) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-0.47163481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971805) q[0];
sx q[0];
rz(-2.4025859) q[0];
sx q[0];
rz(1.7790545) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6064862) q[2];
sx q[2];
rz(-1.9743414) q[2];
sx q[2];
rz(-1.923179) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9593352) q[1];
sx q[1];
rz(-1.253009) q[1];
sx q[1];
rz(0.55303581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0105693) q[3];
sx q[3];
rz(-0.71478715) q[3];
sx q[3];
rz(-0.21372114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.65583324) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(1.2539697) q[2];
rz(-2.8433825) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(-1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-0.75697672) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.6964361) q[1];
sx q[1];
rz(2.7991378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22638527) q[0];
sx q[0];
rz(-1.3801551) q[0];
sx q[0];
rz(-1.1676428) q[0];
x q[1];
rz(-0.32462281) q[2];
sx q[2];
rz(-2.2461736) q[2];
sx q[2];
rz(1.7243232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1736502) q[1];
sx q[1];
rz(-1.0999803) q[1];
sx q[1];
rz(-1.1556666) q[1];
rz(-pi) q[2];
rz(-2.7259521) q[3];
sx q[3];
rz(-1.8438253) q[3];
sx q[3];
rz(1.331783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(2.5214419) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(-1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81724375) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(-1.2212344) q[0];
rz(1.1034032) q[1];
sx q[1];
rz(-1.9788479) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4909523) q[0];
sx q[0];
rz(-1.3155259) q[0];
sx q[0];
rz(1.7351236) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76580183) q[2];
sx q[2];
rz(-0.81379902) q[2];
sx q[2];
rz(1.162093) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6891046) q[1];
sx q[1];
rz(-2.9924722) q[1];
sx q[1];
rz(1.1959056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8437988) q[3];
sx q[3];
rz(-1.2074561) q[3];
sx q[3];
rz(-0.17223528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3409884) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(3.0321339) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(3.0808595) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3133746) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(0.24318801) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(2.5004255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.084448) q[0];
sx q[0];
rz(-2.085745) q[0];
sx q[0];
rz(-0.65339974) q[0];
rz(-pi) q[1];
rz(1.2095318) q[2];
sx q[2];
rz(-2.2881977) q[2];
sx q[2];
rz(-0.016005767) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9037364) q[1];
sx q[1];
rz(-2.3163539) q[1];
sx q[1];
rz(1.115487) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9092714) q[3];
sx q[3];
rz(-1.182175) q[3];
sx q[3];
rz(2.6040833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0082561) q[2];
sx q[2];
rz(-1.1587605) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(0.64822316) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-0.42914036) q[0];
rz(0.40402135) q[1];
sx q[1];
rz(-0.41530135) q[1];
sx q[1];
rz(0.9187575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0149834) q[0];
sx q[0];
rz(-2.1169239) q[0];
sx q[0];
rz(-1.855143) q[0];
rz(-pi) q[1];
rz(0.12523052) q[2];
sx q[2];
rz(-1.8298143) q[2];
sx q[2];
rz(-0.7746402) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4705953) q[1];
sx q[1];
rz(-1.8806702) q[1];
sx q[1];
rz(-2.9122374) q[1];
rz(-pi) q[2];
rz(0.073831255) q[3];
sx q[3];
rz(-1.257535) q[3];
sx q[3];
rz(-1.8262124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1227485) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(1.7913691) q[2];
rz(0.33251479) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(1.7727218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74828446) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(2.6848324) q[0];
rz(-1.7204174) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(-3.0866887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47726705) q[0];
sx q[0];
rz(-2.3033963) q[0];
sx q[0];
rz(-0.038756702) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.730092) q[2];
sx q[2];
rz(-1.1724645) q[2];
sx q[2];
rz(-0.27137941) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0662288) q[1];
sx q[1];
rz(-1.7875331) q[1];
sx q[1];
rz(2.0677057) q[1];
x q[2];
rz(0.53590448) q[3];
sx q[3];
rz(-1.8632938) q[3];
sx q[3];
rz(-1.747594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2533337) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(2.9448275) q[2];
rz(-1.0737859) q[3];
sx q[3];
rz(-2.00627) q[3];
sx q[3];
rz(-0.73608583) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9842904) q[0];
sx q[0];
rz(-2.2881916) q[0];
sx q[0];
rz(0.89944696) q[0];
rz(-0.31117123) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.4003632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104855) q[0];
sx q[0];
rz(-2.541158) q[0];
sx q[0];
rz(-1.1945634) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.898959) q[2];
sx q[2];
rz(-2.0443235) q[2];
sx q[2];
rz(2.1514016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77620095) q[1];
sx q[1];
rz(-1.9579871) q[1];
sx q[1];
rz(-0.64761083) q[1];
rz(-pi) q[2];
rz(0.71615852) q[3];
sx q[3];
rz(-1.0993488) q[3];
sx q[3];
rz(2.6177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(-2.0628085) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-0.88206464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494736) q[0];
sx q[0];
rz(-1.619512) q[0];
sx q[0];
rz(-1.6721538) q[0];
rz(-3.0081765) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-3.1405814) q[2];
sx q[2];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5569729) q[2];
rz(-1.6422938) q[3];
sx q[3];
rz(-1.0701978) q[3];
sx q[3];
rz(1.3842907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
