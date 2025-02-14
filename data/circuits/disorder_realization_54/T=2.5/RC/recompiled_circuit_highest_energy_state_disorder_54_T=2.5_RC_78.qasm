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
rz(-0.32831353) q[0];
sx q[0];
rz(5.1483122) q[0];
sx q[0];
rz(6.8490646) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(-0.53425962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1318783) q[0];
sx q[0];
rz(-1.8410972) q[0];
sx q[0];
rz(2.2755095) q[0];
x q[1];
rz(1.8741605) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(2.0801185) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6001662) q[1];
sx q[1];
rz(-1.2703589) q[1];
sx q[1];
rz(-1.5460127) q[1];
rz(0.45716469) q[3];
sx q[3];
rz(-1.346753) q[3];
sx q[3];
rz(-0.30686298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(-2.8924083) q[2];
rz(-1.3439882) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(-3.0715166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(-2.1517854) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(2.8299423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79651901) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(-0.89181283) q[0];
rz(-0.10820724) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(2.3600887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5144044) q[1];
sx q[1];
rz(-1.3986821) q[1];
sx q[1];
rz(1.3283511) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95590638) q[3];
sx q[3];
rz(-1.2887738) q[3];
sx q[3];
rz(0.36263719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0011757294) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(2.1902093) q[2];
rz(-2.9133255) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(3.0271295) q[0];
rz(2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(0.1618298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069050633) q[0];
sx q[0];
rz(-2.4009573) q[0];
sx q[0];
rz(-2.1669751) q[0];
x q[1];
rz(-2.7527005) q[2];
sx q[2];
rz(-1.4894247) q[2];
sx q[2];
rz(0.55079764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5787631) q[1];
sx q[1];
rz(-2.17318) q[1];
sx q[1];
rz(-1.8627326) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9818241) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(2.4730554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-3.0008345) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(-0.67489135) q[0];
rz(-0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-0.45796576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77112326) q[0];
sx q[0];
rz(-1.1827977) q[0];
sx q[0];
rz(-2.6113135) q[0];
rz(-0.9651009) q[2];
sx q[2];
rz(-0.5083681) q[2];
sx q[2];
rz(1.8338211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1110036) q[1];
sx q[1];
rz(-1.9408424) q[1];
sx q[1];
rz(1.3732598) q[1];
rz(-pi) q[2];
rz(1.6776039) q[3];
sx q[3];
rz(-1.065101) q[3];
sx q[3];
rz(-0.13733521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13350479) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(-2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1595681) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(2.9087861) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(2.1585042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1637835) q[0];
sx q[0];
rz(-2.0931232) q[0];
sx q[0];
rz(-0.99951007) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11396046) q[2];
sx q[2];
rz(-2.1884005) q[2];
sx q[2];
rz(-1.8007939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59683387) q[1];
sx q[1];
rz(-1.7167517) q[1];
sx q[1];
rz(-2.9533141) q[1];
rz(-pi) q[2];
rz(-1.1484723) q[3];
sx q[3];
rz(-1.7230801) q[3];
sx q[3];
rz(-1.1642758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7588707) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(-0.37528428) q[2];
rz(-1.075607) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(-2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4009092) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(1.7042473) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(2.5934503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63614908) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(0.33238073) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0847381) q[2];
sx q[2];
rz(-1.0722188) q[2];
sx q[2];
rz(2.9996897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3972612) q[1];
sx q[1];
rz(-1.9588699) q[1];
sx q[1];
rz(0.71340386) q[1];
rz(-2.3856134) q[3];
sx q[3];
rz(-1.090831) q[3];
sx q[3];
rz(-2.2868613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(-0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-0.4959929) q[0];
rz(1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(2.0248263) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46653173) q[0];
sx q[0];
rz(-1.2197615) q[0];
sx q[0];
rz(-2.768633) q[0];
rz(-pi) q[1];
rz(-2.1082868) q[2];
sx q[2];
rz(-1.337572) q[2];
sx q[2];
rz(1.6526412) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3254814) q[1];
sx q[1];
rz(-1.2524091) q[1];
sx q[1];
rz(2.7834784) q[1];
x q[2];
rz(-3.0159166) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(-2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3108959) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(-2.7867219) q[2];
rz(2.3762459) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(-0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527721) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(-2.3759957) q[0];
rz(-0.87144026) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.3227468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2521556) q[0];
sx q[0];
rz(-1.3103292) q[0];
sx q[0];
rz(-1.5171264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23643146) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(-0.98275358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.821956) q[1];
sx q[1];
rz(-2.0791035) q[1];
sx q[1];
rz(-1.5896569) q[1];
rz(0.20669528) q[3];
sx q[3];
rz(-1.350025) q[3];
sx q[3];
rz(-1.7305935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677277) q[0];
sx q[0];
rz(-1.2465957) q[0];
sx q[0];
rz(-2.8644323) q[0];
rz(-0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.1433196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(-1.8822549) q[0];
rz(-pi) q[1];
rz(-2.9306445) q[2];
sx q[2];
rz(-0.8767414) q[2];
sx q[2];
rz(1.8927434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9710247) q[1];
sx q[1];
rz(-1.8877561) q[1];
sx q[1];
rz(-0.16522878) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36181776) q[3];
sx q[3];
rz(-1.4103977) q[3];
sx q[3];
rz(3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9144168) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(1.0660727) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(3.1220806) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(1.5628409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514141) q[0];
sx q[0];
rz(-1.7008243) q[0];
sx q[0];
rz(-3.0517654) q[0];
rz(-pi) q[1];
rz(-2.1807272) q[2];
sx q[2];
rz(-1.7210519) q[2];
sx q[2];
rz(-3.1203258) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0622618) q[1];
sx q[1];
rz(-0.76594662) q[1];
sx q[1];
rz(-2.2205611) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14871493) q[3];
sx q[3];
rz(-2.2991355) q[3];
sx q[3];
rz(-0.088012849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(-0.16555244) q[2];
rz(2.5340269) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(-0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3792569) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(-1.5070076) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(2.3661939) q[2];
sx q[2];
rz(-1.4682894) q[2];
sx q[2];
rz(-0.52649047) q[2];
rz(-1.7692825) q[3];
sx q[3];
rz(-2.1673421) q[3];
sx q[3];
rz(2.0094595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
