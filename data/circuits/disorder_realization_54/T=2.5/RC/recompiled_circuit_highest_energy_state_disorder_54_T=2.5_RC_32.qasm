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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(-0.56587926) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(2.607333) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.256828) q[0];
sx q[0];
rz(-0.74640154) q[0];
sx q[0];
rz(1.9749667) q[0];
x q[1];
rz(0.59492971) q[2];
sx q[2];
rz(-1.317136) q[2];
sx q[2];
rz(0.67753917) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45785832) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(-3.0617759) q[1];
rz(-pi) q[2];
rz(-0.45716469) q[3];
sx q[3];
rz(-1.7948397) q[3];
sx q[3];
rz(-0.30686298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(-1.3439882) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4942112) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(2.3509332) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(2.8299423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4379398) q[0];
sx q[0];
rz(-0.92267413) q[0];
sx q[0];
rz(2.4841643) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0333854) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(0.78150392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6271882) q[1];
sx q[1];
rz(-1.7429105) q[1];
sx q[1];
rz(-1.3283511) q[1];
rz(-pi) q[2];
rz(-1.1053322) q[3];
sx q[3];
rz(-2.4727945) q[3];
sx q[3];
rz(-1.5835403) q[3];
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
rz(-0.95138335) q[2];
rz(-2.9133255) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(1.2736646) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15762873) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(-0.1618298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67437088) q[0];
sx q[0];
rz(-2.1631952) q[0];
sx q[0];
rz(0.47426275) q[0];
rz(-0.38889216) q[2];
sx q[2];
rz(-1.6521679) q[2];
sx q[2];
rz(-2.590795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0664252) q[1];
sx q[1];
rz(-0.6614162) q[1];
sx q[1];
rz(2.7452046) q[1];
rz(1.9818241) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(-2.4730554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74730045) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(0.14075819) q[2];
rz(0.55177871) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(-2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(2.4667013) q[0];
rz(-2.9871125) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-2.6836269) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58139153) q[0];
sx q[0];
rz(-1.0835674) q[0];
sx q[0];
rz(2.013252) q[0];
rz(-pi) q[1];
rz(0.9651009) q[2];
sx q[2];
rz(-2.6332246) q[2];
sx q[2];
rz(1.8338211) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6667216) q[1];
sx q[1];
rz(-0.41731605) q[1];
sx q[1];
rz(0.46837193) q[1];
rz(0.19019698) q[3];
sx q[3];
rz(-2.6256997) q[3];
sx q[3];
rz(2.7864393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0080879) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(-2.9087861) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(0.98308841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9778091) q[0];
sx q[0];
rz(-2.0931232) q[0];
sx q[0];
rz(0.99951007) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4120502) q[2];
sx q[2];
rz(-0.6266784) q[2];
sx q[2];
rz(1.9959297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8195962) q[1];
sx q[1];
rz(-2.9038972) q[1];
sx q[1];
rz(-2.4758215) q[1];
rz(1.2125254) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(3.0606396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(-2.7663084) q[2];
rz(1.075607) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63614908) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(-2.8092119) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4324722) q[2];
sx q[2];
rz(-2.4599584) q[2];
sx q[2];
rz(2.1644985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3972612) q[1];
sx q[1];
rz(-1.9588699) q[1];
sx q[1];
rz(-0.71340386) q[1];
rz(-2.4924566) q[3];
sx q[3];
rz(-2.2723291) q[3];
sx q[3];
rz(-2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(1.5843377) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976582) q[0];
sx q[0];
rz(-1.1827396) q[0];
sx q[0];
rz(2.6455998) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(2.0248263) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46653173) q[0];
sx q[0];
rz(-1.2197615) q[0];
sx q[0];
rz(-2.768633) q[0];
x q[1];
rz(2.871795) q[2];
sx q[2];
rz(-1.0493663) q[2];
sx q[2];
rz(-0.21873378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36193725) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(1.9091868) q[1];
x q[2];
rz(0.12567606) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(0.3885551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3108959) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(-2.7867219) q[2];
rz(2.3762459) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(-2.2701524) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(-1.3227468) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2521556) q[0];
sx q[0];
rz(-1.8312635) q[0];
sx q[0];
rz(-1.5171264) q[0];
x q[1];
rz(-2.2501588) q[2];
sx q[2];
rz(-1.3855033) q[2];
sx q[2];
rz(0.44033716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28089954) q[1];
sx q[1];
rz(-2.6329663) q[1];
sx q[1];
rz(0.033837576) q[1];
rz(-pi) q[2];
rz(2.3115308) q[3];
sx q[3];
rz(-0.30127159) q[3];
sx q[3];
rz(2.4943182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(2.4634585) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073864989) q[0];
sx q[0];
rz(-1.2465957) q[0];
sx q[0];
rz(0.27716032) q[0];
rz(-0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(-1.998273) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80625929) q[0];
sx q[0];
rz(-1.8758095) q[0];
sx q[0];
rz(0.2103896) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2758946) q[2];
sx q[2];
rz(-1.4091461) q[2];
sx q[2];
rz(-0.18582144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4797693) q[1];
sx q[1];
rz(-0.35614888) q[1];
sx q[1];
rz(-2.0355862) q[1];
x q[2];
rz(2.7797749) q[3];
sx q[3];
rz(-1.4103977) q[3];
sx q[3];
rz(3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2271759) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(-0.048967036) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-2.5132096) q[3];
sx q[3];
rz(-0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0186036) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(1.5628409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21694788) q[0];
sx q[0];
rz(-0.15789444) q[0];
sx q[0];
rz(-2.1720706) q[0];
rz(-pi) q[1];
rz(2.9589505) q[2];
sx q[2];
rz(-2.172875) q[2];
sx q[2];
rz(1.6537567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73247468) q[1];
sx q[1];
rz(-2.1555087) q[1];
sx q[1];
rz(0.52701108) q[1];
x q[2];
rz(1.4061798) q[3];
sx q[3];
rz(-0.74062982) q[3];
sx q[3];
rz(3.0082138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(2.9760402) q[2];
rz(-2.5340269) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.3792569) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
rz(1.634585) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(-0.14590895) q[2];
sx q[2];
rz(-0.78073954) q[2];
sx q[2];
rz(-2.2013046) q[2];
rz(-1.3723102) q[3];
sx q[3];
rz(-0.9742506) q[3];
sx q[3];
rz(-1.1321332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
