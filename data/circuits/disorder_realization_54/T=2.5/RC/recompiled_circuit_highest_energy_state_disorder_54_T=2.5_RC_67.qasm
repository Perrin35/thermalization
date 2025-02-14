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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7843648) q[0];
sx q[0];
rz(-0.89656943) q[0];
sx q[0];
rz(-2.7927464) q[0];
rz(-pi) q[1];
rz(-2.7083394) q[2];
sx q[2];
rz(-2.5009071) q[2];
sx q[2];
rz(-2.6034466) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45785832) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(0.079816743) q[1];
rz(-pi) q[2];
rz(2.665042) q[3];
sx q[3];
rz(-2.635987) q[3];
sx q[3];
rz(-1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(1.7976044) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(-0.98980728) q[0];
rz(2.3509332) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-2.8299423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7036529) q[0];
sx q[0];
rz(-2.2189185) q[0];
sx q[0];
rz(0.65742832) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10820724) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(0.78150392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.098722301) q[1];
sx q[1];
rz(-1.8095892) q[1];
sx q[1];
rz(0.17719321) q[1];
rz(-pi) q[2];
rz(-2.8007224) q[3];
sx q[3];
rz(-0.98353993) q[3];
sx q[3];
rz(-1.0140918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0011757294) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15762873) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(3.0271295) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(-2.9797629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1755668) q[0];
sx q[0];
rz(-1.9593666) q[0];
sx q[0];
rz(0.92312621) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6587018) q[2];
sx q[2];
rz(-1.1832596) q[2];
sx q[2];
rz(-1.0532925) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56282957) q[1];
sx q[1];
rz(-0.96841267) q[1];
sx q[1];
rz(1.8627326) q[1];
rz(2.3409445) q[3];
sx q[3];
rz(-1.8654239) q[3];
sx q[3];
rz(-1.9484474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(3.0008345) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9363339) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(2.4667013) q[0];
rz(0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-2.6836269) q[1];
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
rz(1.1283406) q[0];
x q[1];
rz(-1.141233) q[2];
sx q[2];
rz(-1.2899961) q[2];
sx q[2];
rz(-0.28102885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47487101) q[1];
sx q[1];
rz(-2.7242766) q[1];
sx q[1];
rz(-2.6732207) q[1];
rz(-pi) q[2];
rz(-0.50812085) q[3];
sx q[3];
rz(-1.6641938) q[3];
sx q[3];
rz(-1.3815708) q[3];
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
rz(-2.1250686) q[3];
sx q[3];
rz(2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1595681) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.9200448) q[0];
rz(-2.9087861) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(2.1585042) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096689396) q[0];
sx q[0];
rz(-1.0830729) q[0];
sx q[0];
rz(0.60012598) q[0];
x q[1];
rz(-1.7295425) q[2];
sx q[2];
rz(-2.5149143) q[2];
sx q[2];
rz(-1.145663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1399263) q[1];
sx q[1];
rz(-1.7570494) q[1];
sx q[1];
rz(-1.4222533) q[1];
rz(1.9290673) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(0.080953065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38272196) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(0.37528428) q[2];
rz(-1.075607) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(1.7042473) q[0];
rz(-3.0628693) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054436) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(-0.33238073) q[0];
x q[1];
rz(2.5896543) q[2];
sx q[2];
rz(-1.1480398) q[2];
sx q[2];
rz(-1.960159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9029451) q[1];
sx q[1];
rz(-0.79552197) q[1];
sx q[1];
rz(-2.5832157) q[1];
rz(-pi) q[2];
rz(0.6491361) q[3];
sx q[3];
rz(-2.2723291) q[3];
sx q[3];
rz(-2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(-1.5572549) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(-1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(-2.0248263) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6750609) q[0];
sx q[0];
rz(-1.2197615) q[0];
sx q[0];
rz(0.37295966) q[0];
rz(-pi) q[1];
rz(1.1363813) q[2];
sx q[2];
rz(-0.58131733) q[2];
sx q[2];
rz(2.8535247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45124829) q[1];
sx q[1];
rz(-0.4745698) q[1];
sx q[1];
rz(0.75466538) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0159166) q[3];
sx q[3];
rz(-2.6824441) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(-0.35487077) q[2];
rz(0.76534671) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(1.3227468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3324748) q[0];
sx q[0];
rz(-1.6226543) q[0];
sx q[0];
rz(-2.8807667) q[0];
x q[1];
rz(2.9051612) q[2];
sx q[2];
rz(-2.2364103) q[2];
sx q[2];
rz(-0.98275358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.1077551) q[1];
rz(-pi) q[2];
rz(2.3115308) q[3];
sx q[3];
rz(-2.8403211) q[3];
sx q[3];
rz(0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(-1.0929557) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(-0.67813412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(2.8644323) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(1.998273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(1.2593377) q[0];
rz(-2.2758946) q[2];
sx q[2];
rz(-1.7324466) q[2];
sx q[2];
rz(-0.18582144) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7932897) q[1];
sx q[1];
rz(-1.7277246) q[1];
sx q[1];
rz(1.8918397) q[1];
rz(-pi) q[2];
rz(1.7420898) q[3];
sx q[3];
rz(-1.2138324) q[3];
sx q[3];
rz(-1.5689284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0186036) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(-2.3628269) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5628409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7514141) q[0];
sx q[0];
rz(-1.4407684) q[0];
sx q[0];
rz(-3.0517654) q[0];
rz(-0.96086545) q[2];
sx q[2];
rz(-1.4205407) q[2];
sx q[2];
rz(0.021266887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.409118) q[1];
sx q[1];
rz(-2.1555087) q[1];
sx q[1];
rz(2.6145816) q[1];
rz(-pi) q[2];
rz(-1.4061798) q[3];
sx q[3];
rz(-0.74062982) q[3];
sx q[3];
rz(0.1333789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(0.16555244) q[2];
rz(-0.60756573) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3792569) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(1.5070076) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(0.14590895) q[2];
sx q[2];
rz(-2.3608531) q[2];
sx q[2];
rz(0.94028801) q[2];
rz(-2.535798) q[3];
sx q[3];
rz(-1.7346564) q[3];
sx q[3];
rz(0.32614542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
