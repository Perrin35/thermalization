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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3572278) q[0];
sx q[0];
rz(-0.89656943) q[0];
sx q[0];
rz(-2.7927464) q[0];
rz(-pi) q[1];
rz(-0.43325328) q[2];
sx q[2];
rz(-0.64068551) q[2];
sx q[2];
rz(-2.6034466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54142648) q[1];
sx q[1];
rz(-1.2703589) q[1];
sx q[1];
rz(-1.59558) q[1];
rz(-pi) q[2];
x q[2];
rz(2.665042) q[3];
sx q[3];
rz(-0.50560564) q[3];
sx q[3];
rz(-1.6881534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1285706) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(-0.24918431) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4942112) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(-0.98980728) q[0];
rz(-0.79065943) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(-0.31165037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30325952) q[0];
sx q[0];
rz(-1.0619535) q[0];
sx q[0];
rz(-2.3340387) q[0];
rz(-pi) q[1];
rz(3.0333854) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(-0.78150392) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5492925) q[1];
sx q[1];
rz(-0.29634297) q[1];
sx q[1];
rz(-2.1974988) q[1];
x q[2];
rz(0.95590638) q[3];
sx q[3];
rz(-1.8528189) q[3];
sx q[3];
rz(2.7789555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(2.9133255) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(-1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(3.0271295) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(0.1618298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4672218) q[0];
sx q[0];
rz(-2.1631952) q[0];
sx q[0];
rz(0.47426275) q[0];
rz(1.6587018) q[2];
sx q[2];
rz(-1.958333) q[2];
sx q[2];
rz(2.0883002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0664252) q[1];
sx q[1];
rz(-0.6614162) q[1];
sx q[1];
rz(-2.7452046) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7416218) q[3];
sx q[3];
rz(-2.2999527) q[3];
sx q[3];
rz(-0.10310452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74730045) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-3.0008345) q[2];
rz(-0.55177871) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9363339) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(2.4667013) q[0];
rz(2.9871125) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-0.45796576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5602011) q[0];
sx q[0];
rz(-1.0835674) q[0];
sx q[0];
rz(1.1283406) q[0];
x q[1];
rz(2.0003597) q[2];
sx q[2];
rz(-1.8515966) q[2];
sx q[2];
rz(-2.8605638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.612466) q[1];
sx q[1];
rz(-1.3867897) q[1];
sx q[1];
rz(2.7648768) q[1];
x q[2];
rz(0.50812085) q[3];
sx q[3];
rz(-1.6641938) q[3];
sx q[3];
rz(-1.7600218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13350479) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(-1.3301298) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(-2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1595681) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.2215479) q[0];
rz(-0.23280652) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(-0.98308841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0449033) q[0];
sx q[0];
rz(-1.0830729) q[0];
sx q[0];
rz(-0.60012598) q[0];
rz(-1.4120502) q[2];
sx q[2];
rz(-2.5149143) q[2];
sx q[2];
rz(-1.9959297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1399263) q[1];
sx q[1];
rz(-1.7570494) q[1];
sx q[1];
rz(-1.4222533) q[1];
x q[2];
rz(-2.9748989) q[3];
sx q[3];
rz(-1.1536666) q[3];
sx q[3];
rz(0.47458173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38272196) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(-2.7663084) q[2];
rz(2.0659857) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4009092) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(2.5934503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1931216) q[0];
sx q[0];
rz(-1.7821494) q[0];
sx q[0];
rz(2.4708492) q[0];
rz(-2.4324722) q[2];
sx q[2];
rz(-0.68163423) q[2];
sx q[2];
rz(0.97709419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50997616) q[1];
sx q[1];
rz(-0.92008925) q[1];
sx q[1];
rz(2.0664311) q[1];
rz(0.75597922) q[3];
sx q[3];
rz(-2.0507617) q[3];
sx q[3];
rz(2.2868613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(-1.5572549) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5439344) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(1.1167663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46653173) q[0];
sx q[0];
rz(-1.9218311) q[0];
sx q[0];
rz(2.768633) q[0];
rz(-pi) q[1];
rz(2.1082868) q[2];
sx q[2];
rz(-1.337572) q[2];
sx q[2];
rz(-1.6526412) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7796554) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(1.9091868) q[1];
x q[2];
rz(-0.12567606) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(-0.3885551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(2.7867219) q[2];
rz(0.76534671) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(-3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(0.76559693) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(1.8188459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0950892) q[0];
sx q[0];
rz(-2.8757767) q[0];
sx q[0];
rz(-2.9429716) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23643146) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(-0.98275358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28089954) q[1];
sx q[1];
rz(-2.6329663) q[1];
sx q[1];
rz(-3.1077551) q[1];
rz(-1.345383) q[3];
sx q[3];
rz(-1.3691878) q[3];
sx q[3];
rz(0.20568337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65779296) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(1.0929557) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(-2.4634585) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-0.27716032) q[0];
rz(-0.74835888) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(1.998273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4410981) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(1.2593377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21094811) q[2];
sx q[2];
rz(-2.2648513) q[2];
sx q[2];
rz(-1.8927434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9710247) q[1];
sx q[1];
rz(-1.2538365) q[1];
sx q[1];
rz(2.9763639) q[1];
x q[2];
rz(-2.7797749) q[3];
sx q[3];
rz(-1.731195) q[3];
sx q[3];
rz(-0.058506207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2271759) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(-3.0926256) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-0.019512026) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5628409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21694788) q[0];
sx q[0];
rz(-0.15789444) q[0];
sx q[0];
rz(2.1720706) q[0];
x q[1];
rz(2.1807272) q[2];
sx q[2];
rz(-1.7210519) q[2];
sx q[2];
rz(-0.021266887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1491039) q[1];
sx q[1];
rz(-2.0035775) q[1];
sx q[1];
rz(2.2243568) q[1];
rz(-1.7354129) q[3];
sx q[3];
rz(-2.4009628) q[3];
sx q[3];
rz(0.1333789) q[3];
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
rz(2.5340269) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.7623357) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(-1.634585) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(-1.4277369) q[2];
sx q[2];
rz(-0.80052994) q[2];
sx q[2];
rz(1.1442727) q[2];
rz(0.28260091) q[3];
sx q[3];
rz(-2.5167214) q[3];
sx q[3];
rz(-1.4759397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
