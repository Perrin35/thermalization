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
rz(2.5757134) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(-1.3574418) q[1];
sx q[1];
rz(0.53425962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847647) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(1.9749667) q[0];
x q[1];
rz(-1.2674321) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(-1.0614741) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6001662) q[1];
sx q[1];
rz(-1.2703589) q[1];
sx q[1];
rz(-1.59558) q[1];
rz(-pi) q[2];
rz(-2.665042) q[3];
sx q[3];
rz(-2.635987) q[3];
sx q[3];
rz(1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1285706) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(-2.8924083) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-3.0715166) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473815) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(2.1517854) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(2.8299423) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383331) q[0];
sx q[0];
rz(-2.0796392) q[0];
sx q[0];
rz(2.3340387) q[0];
rz(-pi) q[1];
rz(1.7575532) q[2];
sx q[2];
rz(-1.6771363) q[2];
sx q[2];
rz(0.76923907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6271882) q[1];
sx q[1];
rz(-1.3986821) q[1];
sx q[1];
rz(1.3283511) q[1];
rz(-pi) q[2];
rz(-0.95590638) q[3];
sx q[3];
rz(-1.2887738) q[3];
sx q[3];
rz(2.7789555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1404169) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(-1.2736646) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15762873) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(-0.1618298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.072542) q[0];
sx q[0];
rz(-0.74063535) q[0];
sx q[0];
rz(-0.97461755) q[0];
rz(-2.7527005) q[2];
sx q[2];
rz(-1.4894247) q[2];
sx q[2];
rz(-2.590795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.075167478) q[1];
sx q[1];
rz(-0.6614162) q[1];
sx q[1];
rz(2.7452046) q[1];
x q[2];
rz(2.3409445) q[3];
sx q[3];
rz(-1.8654239) q[3];
sx q[3];
rz(1.1931452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(-2.4667013) q[0];
rz(0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(0.45796576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7688528) q[0];
sx q[0];
rz(-0.64592664) q[0];
sx q[0];
rz(0.67966184) q[0];
x q[1];
rz(0.30720023) q[2];
sx q[2];
rz(-1.9824902) q[2];
sx q[2];
rz(1.1634942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.612466) q[1];
sx q[1];
rz(-1.3867897) q[1];
sx q[1];
rz(-0.37671582) q[1];
rz(1.6776039) q[3];
sx q[3];
rz(-2.0764917) q[3];
sx q[3];
rz(0.13733521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(0.31944719) q[2];
rz(1.3301298) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(-2.5813848) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.2215479) q[0];
rz(-0.23280652) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(-2.1585042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670416) q[0];
sx q[0];
rz(-2.3877151) q[0];
sx q[0];
rz(-0.75410944) q[0];
rz(-pi) q[1];
rz(1.7295425) q[2];
sx q[2];
rz(-2.5149143) q[2];
sx q[2];
rz(-1.9959297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1399263) q[1];
sx q[1];
rz(-1.3845433) q[1];
sx q[1];
rz(-1.4222533) q[1];
rz(-pi) q[2];
rz(-1.9290673) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(3.0606396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(-0.37528428) q[2];
rz(-2.0659857) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(-2.5934503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9288612) q[0];
sx q[0];
rz(-2.2240046) q[0];
sx q[0];
rz(1.3034588) q[0];
rz(-pi) q[1];
rz(-0.55193837) q[2];
sx q[2];
rz(-1.1480398) q[2];
sx q[2];
rz(-1.960159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50997616) q[1];
sx q[1];
rz(-2.2215034) q[1];
sx q[1];
rz(1.0751616) q[1];
rz(2.4924566) q[3];
sx q[3];
rz(-2.2723291) q[3];
sx q[3];
rz(-0.26065206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(0.6753298) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(-0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(2.6455998) q[0];
rz(1.4542106) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(-2.0248263) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9705212) q[0];
sx q[0];
rz(-1.9200033) q[0];
sx q[0];
rz(-1.9454576) q[0];
x q[1];
rz(0.26979763) q[2];
sx q[2];
rz(-2.0922263) q[2];
sx q[2];
rz(-0.21873378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3254814) q[1];
sx q[1];
rz(-1.8891836) q[1];
sx q[1];
rz(2.7834784) q[1];
rz(-pi) q[2];
rz(2.6855822) q[3];
sx q[3];
rz(-1.6263762) q[3];
sx q[3];
rz(-2.0721276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(0.76534671) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4527721) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(1.8188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046503456) q[0];
sx q[0];
rz(-2.8757767) q[0];
sx q[0];
rz(2.9429716) q[0];
rz(1.2808676) q[2];
sx q[2];
rz(-2.4413042) q[2];
sx q[2];
rz(1.354745) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28089954) q[1];
sx q[1];
rz(-0.50862632) q[1];
sx q[1];
rz(-0.033837576) q[1];
rz(0.20669528) q[3];
sx q[3];
rz(-1.7915676) q[3];
sx q[3];
rz(1.7305935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(-0.22171177) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-0.67813412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0677277) q[0];
sx q[0];
rz(-1.2465957) q[0];
sx q[0];
rz(-0.27716032) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.998273) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4244014) q[0];
sx q[0];
rz(-2.7729308) q[0];
sx q[0];
rz(-2.1564846) q[0];
rz(-0.86569806) q[2];
sx q[2];
rz(-1.4091461) q[2];
sx q[2];
rz(2.9557712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34830293) q[1];
sx q[1];
rz(-1.4138681) q[1];
sx q[1];
rz(1.2497529) q[1];
rz(-pi) q[2];
rz(0.42871126) q[3];
sx q[3];
rz(-2.747251) q[3];
sx q[3];
rz(-2.0286146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2271759) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(-3.0926256) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(-2.9858885) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(2.3628269) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5787517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9492968) q[0];
sx q[0];
rz(-1.4817294) q[0];
sx q[0];
rz(1.7013447) q[0];
rz(-pi) q[1];
rz(-1.3124002) q[2];
sx q[2];
rz(-0.62587291) q[2];
sx q[2];
rz(-1.3385119) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.07933087) q[1];
sx q[1];
rz(-0.76594662) q[1];
sx q[1];
rz(2.2205611) q[1];
rz(1.4061798) q[3];
sx q[3];
rz(-0.74062982) q[3];
sx q[3];
rz(3.0082138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(2.9760402) q[2];
rz(-2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.634585) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(0.77539879) q[2];
sx q[2];
rz(-1.6733032) q[2];
sx q[2];
rz(2.6151022) q[2];
rz(2.8589917) q[3];
sx q[3];
rz(-0.62487124) q[3];
sx q[3];
rz(1.665653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
