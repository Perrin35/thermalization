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
rz(0.24428754) q[1];
sx q[1];
rz(4.9257435) q[1];
sx q[1];
rz(9.9590376) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0097144011) q[0];
sx q[0];
rz(-1.8410972) q[0];
sx q[0];
rz(-2.2755095) q[0];
x q[1];
rz(-0.59492971) q[2];
sx q[2];
rz(-1.8244566) q[2];
sx q[2];
rz(-2.4640535) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.104887) q[1];
sx q[1];
rz(-1.547123) q[1];
sx q[1];
rz(0.30052431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45716469) q[3];
sx q[3];
rz(-1.346753) q[3];
sx q[3];
rz(-2.8347297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-3.0715166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(-2.1517854) q[0];
rz(-0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-2.8299423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4379398) q[0];
sx q[0];
rz(-2.2189185) q[0];
sx q[0];
rz(-0.65742832) q[0];
rz(-0.10820724) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(2.3600887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5923002) q[1];
sx q[1];
rz(-0.29634297) q[1];
sx q[1];
rz(-0.94409385) q[1];
rz(-pi) q[2];
rz(2.0362605) q[3];
sx q[3];
rz(-0.66879818) q[3];
sx q[3];
rz(-1.5580524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1404169) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(2.9133255) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(-1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15762873) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(-0.11446318) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(-2.9797629) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9660258) q[0];
sx q[0];
rz(-1.9593666) q[0];
sx q[0];
rz(-2.2184664) q[0];
x q[1];
rz(2.7527005) q[2];
sx q[2];
rz(-1.4894247) q[2];
sx q[2];
rz(2.590795) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.964965) q[1];
sx q[1];
rz(-1.3313659) q[1];
sx q[1];
rz(0.62271948) q[1];
rz(-pi) q[2];
rz(-1.9818241) q[3];
sx q[3];
rz(-2.3279802) q[3];
sx q[3];
rz(0.66853722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74730045) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-0.14075819) q[2];
rz(0.55177871) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(-0.67489135) q[0];
rz(2.9871125) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(0.45796576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602011) q[0];
sx q[0];
rz(-2.0580252) q[0];
sx q[0];
rz(-2.013252) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1764917) q[2];
sx q[2];
rz(-0.5083681) q[2];
sx q[2];
rz(-1.8338211) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47487101) q[1];
sx q[1];
rz(-2.7242766) q[1];
sx q[1];
rz(-0.46837193) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9513957) q[3];
sx q[3];
rz(-0.51589291) q[3];
sx q[3];
rz(0.35515337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(-1.8114629) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-2.5722952) q[1];
sx q[1];
rz(2.1585042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670416) q[0];
sx q[0];
rz(-2.3877151) q[0];
sx q[0];
rz(-0.75410944) q[0];
rz(3.0276322) q[2];
sx q[2];
rz(-2.1884005) q[2];
sx q[2];
rz(1.3407988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1399263) q[1];
sx q[1];
rz(-1.7570494) q[1];
sx q[1];
rz(-1.4222533) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9931204) q[3];
sx q[3];
rz(-1.7230801) q[3];
sx q[3];
rz(-1.1642758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7588707) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(0.37528428) q[2];
rz(2.0659857) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4009092) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5054436) q[0];
sx q[0];
rz(-0.6983122) q[0];
sx q[0];
rz(-0.33238073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4324722) q[2];
sx q[2];
rz(-2.4599584) q[2];
sx q[2];
rz(-2.1644985) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74433148) q[1];
sx q[1];
rz(-1.9588699) q[1];
sx q[1];
rz(0.71340386) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3856134) q[3];
sx q[3];
rz(-2.0507617) q[3];
sx q[3];
rz(-2.2868613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(0.6753298) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(0.61711446) q[3];
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
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
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
rz(0.46653173) q[0];
sx q[0];
rz(-1.9218311) q[0];
sx q[0];
rz(0.37295966) q[0];
rz(-pi) q[1];
rz(2.1082868) q[2];
sx q[2];
rz(-1.8040207) q[2];
sx q[2];
rz(-1.4889515) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6903444) q[1];
sx q[1];
rz(-0.4745698) q[1];
sx q[1];
rz(-0.75466538) q[1];
rz(-0.12567606) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(-2.3762459) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(-0.76559693) q[0];
rz(0.87144026) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(-1.3227468) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8894371) q[0];
sx q[0];
rz(-1.3103292) q[0];
sx q[0];
rz(-1.6244662) q[0];
x q[1];
rz(0.23643146) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(-0.98275358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8812528) q[1];
sx q[1];
rz(-1.5872721) q[1];
sx q[1];
rz(0.50838281) q[1];
x q[2];
rz(-1.345383) q[3];
sx q[3];
rz(-1.7724049) q[3];
sx q[3];
rz(2.9359093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(-0.22171177) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(0.67813412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-2.6942418) q[1];
sx q[1];
rz(-1.998273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4244014) q[0];
sx q[0];
rz(-0.36866185) q[0];
sx q[0];
rz(-2.1564846) q[0];
x q[1];
rz(1.3242993) q[2];
sx q[2];
rz(-2.4213104) q[2];
sx q[2];
rz(1.5718854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9710247) q[1];
sx q[1];
rz(-1.2538365) q[1];
sx q[1];
rz(0.16522878) q[1];
x q[2];
rz(-2.7797749) q[3];
sx q[3];
rz(-1.4103977) q[3];
sx q[3];
rz(-3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9144168) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(-3.0926256) q[2];
rz(-1.0660727) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.5787517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9492968) q[0];
sx q[0];
rz(-1.6598633) q[0];
sx q[0];
rz(-1.7013447) q[0];
x q[1];
rz(1.8291924) q[2];
sx q[2];
rz(-0.62587291) q[2];
sx q[2];
rz(-1.3385119) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1491039) q[1];
sx q[1];
rz(-2.0035775) q[1];
sx q[1];
rz(-2.2243568) q[1];
rz(-0.14871493) q[3];
sx q[3];
rz(-2.2991355) q[3];
sx q[3];
rz(0.088012849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.94893) q[2];
sx q[2];
rz(-2.2781585) q[2];
sx q[2];
rz(-2.9760402) q[2];
rz(0.60756573) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7623357) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(1.634585) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(2.9956837) q[2];
sx q[2];
rz(-0.78073954) q[2];
sx q[2];
rz(-2.2013046) q[2];
rz(2.535798) q[3];
sx q[3];
rz(-1.4069362) q[3];
sx q[3];
rz(-2.8154472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
