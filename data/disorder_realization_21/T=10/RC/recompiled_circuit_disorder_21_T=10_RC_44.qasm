OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48243385) q[0];
sx q[0];
rz(-1.4855488) q[0];
sx q[0];
rz(-1.8929338) q[0];
rz(2.5311004) q[2];
sx q[2];
rz(-0.98449003) q[2];
sx q[2];
rz(2.5551978) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1134125) q[1];
sx q[1];
rz(-0.98787381) q[1];
sx q[1];
rz(-2.7229573) q[1];
rz(-pi) q[2];
rz(-0.46408848) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(2.1634963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(-2.5855529) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125688) q[0];
sx q[0];
rz(-1.7984263) q[0];
sx q[0];
rz(0.17507041) q[0];
rz(-1.3729587) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(2.3193662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.038139548) q[1];
sx q[1];
rz(-0.98587576) q[1];
sx q[1];
rz(2.1409047) q[1];
rz(-1.7440967) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(-1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.8781352) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(-0.24761565) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945064) q[0];
sx q[0];
rz(-1.9425834) q[0];
sx q[0];
rz(1.0466172) q[0];
rz(-pi) q[1];
rz(-2.1842314) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(-0.36188175) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-1.0291568) q[1];
sx q[1];
rz(0.72802131) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2948202) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(1.5848974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(-0.49989191) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.8968556) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8053999) q[0];
sx q[0];
rz(-1.5691783) q[0];
sx q[0];
rz(2.8021115) q[0];
rz(-pi) q[1];
rz(-1.3696026) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(0.97845562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9089531) q[1];
sx q[1];
rz(-2.5563572) q[1];
sx q[1];
rz(1.1802243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3454868) q[3];
sx q[3];
rz(-1.011519) q[3];
sx q[3];
rz(-2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.9909031) q[2];
rz(1.4771279) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-0.081469014) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1254079) q[0];
sx q[0];
rz(-0.59690079) q[0];
sx q[0];
rz(-1.3422658) q[0];
rz(-pi) q[1];
rz(0.91707768) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(1.9249141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4245783) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(-1.8841519) q[1];
rz(-1.8893858) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(-0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-2.6200263) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(0.7985324) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0621588) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(-2.9955203) q[0];
x q[1];
rz(-1.2217667) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(0.40086056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3244464) q[1];
sx q[1];
rz(-1.6231316) q[1];
sx q[1];
rz(-0.99919341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(-0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(0.46404776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-0.95373017) q[0];
sx q[0];
rz(-2.6130136) q[0];
rz(-pi) q[1];
rz(1.6527962) q[2];
sx q[2];
rz(-0.24157761) q[2];
sx q[2];
rz(-1.8919924) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0276427) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(0.34160683) q[1];
x q[2];
rz(0.23630996) q[3];
sx q[3];
rz(-1.4910306) q[3];
sx q[3];
rz(-2.7304756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(-3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.7000748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38948108) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.585929) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61261119) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(-2.999246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.085233363) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(-1.4443881) q[1];
rz(-2.8076595) q[3];
sx q[3];
rz(-1.5302368) q[3];
sx q[3];
rz(-0.12727748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(0.98659602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6331659) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(-1.3915865) q[0];
rz(-pi) q[1];
rz(2.6096452) q[2];
sx q[2];
rz(-2.1956458) q[2];
sx q[2];
rz(0.22235409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36531891) q[1];
sx q[1];
rz(-0.96818189) q[1];
sx q[1];
rz(1.7425294) q[1];
rz(0.17296965) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(0.017661974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(2.731936) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.4155037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9916519) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.3209692) q[0];
rz(-0.98244046) q[2];
sx q[2];
rz(-2.7124565) q[2];
sx q[2];
rz(2.2619801) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0701329) q[1];
sx q[1];
rz(-1.7839285) q[1];
sx q[1];
rz(0.15503426) q[1];
rz(-pi) q[2];
rz(-0.16898445) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(0.60097354) q[2];
sx q[2];
rz(-0.31243159) q[2];
sx q[2];
rz(2.5675788) q[2];
rz(0.23871213) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];