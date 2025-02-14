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
rz(0.24428754) q[1];
sx q[1];
rz(4.9257435) q[1];
sx q[1];
rz(9.9590376) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256828) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(-1.166626) q[0];
x q[1];
rz(-1.2674321) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(-1.0614741) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54142648) q[1];
sx q[1];
rz(-1.8712338) q[1];
sx q[1];
rz(1.59558) q[1];
rz(-2.684428) q[3];
sx q[3];
rz(-1.7948397) q[3];
sx q[3];
rz(-2.8347297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(-1.7976044) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-3.0715166) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(-2.3509332) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(-2.8299423) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79651901) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(-2.2497798) q[0];
rz(-1.0490408) q[2];
sx q[2];
rz(-2.9269896) q[2];
sx q[2];
rz(-1.3133446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0428704) q[1];
sx q[1];
rz(-1.8095892) q[1];
sx q[1];
rz(-2.9643994) q[1];
rz(-pi) q[2];
rz(-1.1053322) q[3];
sx q[3];
rz(-2.4727945) q[3];
sx q[3];
rz(1.5580524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(2.1902093) q[2];
rz(-2.9133255) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(1.2736646) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(2.9797629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4672218) q[0];
sx q[0];
rz(-0.97839744) q[0];
sx q[0];
rz(2.6673299) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38889216) q[2];
sx q[2];
rz(-1.6521679) q[2];
sx q[2];
rz(-2.590795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5787631) q[1];
sx q[1];
rz(-0.96841267) q[1];
sx q[1];
rz(-1.8627326) q[1];
rz(2.7416218) q[3];
sx q[3];
rz(-2.2999527) q[3];
sx q[3];
rz(3.0384881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(3.0008345) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(-0.67489135) q[0];
rz(-2.9871125) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(-0.45796576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727399) q[0];
sx q[0];
rz(-2.495666) q[0];
sx q[0];
rz(-0.67966184) q[0];
rz(-pi) q[1];
rz(-2.8343924) q[2];
sx q[2];
rz(-1.9824902) q[2];
sx q[2];
rz(1.1634942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47487101) q[1];
sx q[1];
rz(-0.41731605) q[1];
sx q[1];
rz(2.6732207) q[1];
x q[2];
rz(-2.6334718) q[3];
sx q[3];
rz(-1.6641938) q[3];
sx q[3];
rz(-1.7600218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.9820246) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(-1.2215479) q[0];
rz(0.23280652) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(0.98308841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0449033) q[0];
sx q[0];
rz(-2.0585198) q[0];
sx q[0];
rz(-2.5414667) q[0];
rz(-1.4120502) q[2];
sx q[2];
rz(-0.6266784) q[2];
sx q[2];
rz(1.9959297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59683387) q[1];
sx q[1];
rz(-1.4248409) q[1];
sx q[1];
rz(-2.9533141) q[1];
rz(1.9931204) q[3];
sx q[3];
rz(-1.4185126) q[3];
sx q[3];
rz(1.1642758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(-2.7663084) q[2];
rz(1.075607) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(-2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(-0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(-0.54814235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054436) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(2.8092119) q[0];
rz(-pi) q[1];
rz(2.4324722) q[2];
sx q[2];
rz(-2.4599584) q[2];
sx q[2];
rz(-2.1644985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2386475) q[1];
sx q[1];
rz(-0.79552197) q[1];
sx q[1];
rz(2.5832157) q[1];
rz(-pi) q[2];
rz(-0.75597922) q[3];
sx q[3];
rz(-2.0507617) q[3];
sx q[3];
rz(-2.2868613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(0.6753298) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(-1.687382) q[1];
sx q[1];
rz(-1.0554689) q[1];
sx q[1];
rz(1.1167663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46653173) q[0];
sx q[0];
rz(-1.9218311) q[0];
sx q[0];
rz(-0.37295966) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1082868) q[2];
sx q[2];
rz(-1.8040207) q[2];
sx q[2];
rz(-1.6526412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7796554) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(1.2324059) q[1];
rz(-pi) q[2];
rz(-0.45601042) q[3];
sx q[3];
rz(-1.5152165) q[3];
sx q[3];
rz(2.0721276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3108959) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(2.3762459) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(-2.3759957) q[0];
rz(-2.2701524) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(1.8188459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046503456) q[0];
sx q[0];
rz(-0.26581598) q[0];
sx q[0];
rz(-0.19862108) q[0];
rz(1.8607251) q[2];
sx q[2];
rz(-2.4413042) q[2];
sx q[2];
rz(1.7868477) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8606931) q[1];
sx q[1];
rz(-2.6329663) q[1];
sx q[1];
rz(-0.033837576) q[1];
rz(1.345383) q[3];
sx q[3];
rz(-1.3691878) q[3];
sx q[3];
rz(2.9359093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4837997) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(-1.0929557) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073864989) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(0.27716032) q[0];
rz(-0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.1433196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(1.2593377) q[0];
rz(0.86569806) q[2];
sx q[2];
rz(-1.4091461) q[2];
sx q[2];
rz(0.18582144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7932897) q[1];
sx q[1];
rz(-1.7277246) q[1];
sx q[1];
rz(1.8918397) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3995029) q[3];
sx q[3];
rz(-1.2138324) q[3];
sx q[3];
rz(1.5726643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(3.0926256) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-2.5132096) q[3];
sx q[3];
rz(2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12298909) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5628409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9246448) q[0];
sx q[0];
rz(-0.15789444) q[0];
sx q[0];
rz(0.96952207) q[0];
rz(-pi) q[1];
rz(-2.9589505) q[2];
sx q[2];
rz(-0.96871766) q[2];
sx q[2];
rz(-1.4878359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9924888) q[1];
sx q[1];
rz(-1.1380151) q[1];
sx q[1];
rz(-0.91723587) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83694038) q[3];
sx q[3];
rz(-1.4599953) q[3];
sx q[3];
rz(-1.5821804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19266263) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(-0.16555244) q[2];
rz(0.60756573) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.5070076) q[1];
sx q[1];
rz(-1.6804809) q[1];
sx q[1];
rz(1.5483821) q[1];
rz(0.14590895) q[2];
sx q[2];
rz(-2.3608531) q[2];
sx q[2];
rz(0.94028801) q[2];
rz(-0.60579469) q[3];
sx q[3];
rz(-1.4069362) q[3];
sx q[3];
rz(-2.8154472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
