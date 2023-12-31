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
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6591588) q[0];
sx q[0];
rz(-1.4855488) q[0];
sx q[0];
rz(-1.2486588) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88959496) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(1.7878469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-2.4383713) q[1];
sx q[1];
rz(-1.0183079) q[1];
x q[2];
rz(-2.6775042) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(-2.1634963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(2.1957943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6933724) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5599247) q[0];
sx q[0];
rz(-1.4002869) q[0];
sx q[0];
rz(-1.8018363) q[0];
rz(-pi) q[1];
rz(-0.63983812) q[2];
sx q[2];
rz(-0.32391641) q[2];
sx q[2];
rz(-1.4878291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8728767) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(-0.66653911) q[1];
x q[2];
rz(1.7440967) q[3];
sx q[3];
rz(-2.0850075) q[3];
sx q[3];
rz(1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2382425) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.7984614) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1789078) q[0];
sx q[0];
rz(-2.509153) q[0];
sx q[0];
rz(2.2326367) q[0];
rz(-1.1142251) q[2];
sx q[2];
rz(-0.66514981) q[2];
sx q[2];
rz(2.3014724) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3232705) q[1];
sx q[1];
rz(-0.96410492) q[1];
sx q[1];
rz(2.2491749) q[1];
rz(2.2504911) q[3];
sx q[3];
rz(-2.292233) q[3];
sx q[3];
rz(2.6654055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.2447371) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361928) q[0];
sx q[0];
rz(-1.5691783) q[0];
sx q[0];
rz(2.8021115) q[0];
rz(-0.130529) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(-0.56632698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-2.1069063) q[1];
sx q[1];
rz(-2.8944573) q[1];
rz(-pi) q[2];
rz(1.3454868) q[3];
sx q[3];
rz(-1.011519) q[3];
sx q[3];
rz(2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.9909031) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.5030456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7446639) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(0.98608195) q[0];
rz(-pi) q[1];
rz(-1.6790752) q[2];
sx q[2];
rz(-1.4881954) q[2];
sx q[2];
rz(0.29512197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(1.8841519) q[1];
x q[2];
rz(2.6322707) q[3];
sx q[3];
rz(-1.29042) q[3];
sx q[3];
rz(-2.544739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6333255) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-0.26671985) q[0];
rz(0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-0.7985324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.417056) q[0];
sx q[0];
rz(-0.31053156) q[0];
sx q[0];
rz(2.0470371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49662923) q[2];
sx q[2];
rz(-1.8804669) q[2];
sx q[2];
rz(1.8079828) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3244464) q[1];
sx q[1];
rz(-1.6231316) q[1];
sx q[1];
rz(0.99919341) q[1];
rz(2.3409307) q[3];
sx q[3];
rz(-1.1122397) q[3];
sx q[3];
rz(2.3199905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.8926228) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2460829) q[0];
sx q[0];
rz(-1.9946788) q[0];
sx q[0];
rz(0.8830107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3299994) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(0.400825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5210919) q[1];
sx q[1];
rz(-1.2695841) q[1];
sx q[1];
rz(2.0786933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(-1.1788648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15526074) q[0];
sx q[0];
rz(-3.0763456) q[0];
sx q[0];
rz(-2.9078527) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61261119) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(-2.999246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4691094) q[1];
sx q[1];
rz(-1.4454578) q[1];
sx q[1];
rz(3.010716) q[1];
rz(-2.8076595) q[3];
sx q[3];
rz(-1.5302368) q[3];
sx q[3];
rz(-0.12727748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-0.98659602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9771377) q[0];
sx q[0];
rz(-1.7286321) q[0];
sx q[0];
rz(-0.49789238) q[0];
x q[1];
rz(0.95790205) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(0.56607841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(2.5320401) q[1];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.7289366) q[3];
sx q[3];
rz(-1.4827673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9916519) q[0];
sx q[0];
rz(-1.9016001) q[0];
sx q[0];
rz(1.8206235) q[0];
rz(2.1591522) q[2];
sx q[2];
rz(-2.7124565) q[2];
sx q[2];
rz(-0.8796126) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.532383) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(1.7864368) q[1];
rz(-pi) q[2];
rz(-1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(0.48081276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(3.0174875) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-2.4805099) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(1.7514501) q[2];
sx q[2];
rz(-1.8271108) q[2];
sx q[2];
rz(1.9432632) q[2];
rz(-2.5849871) q[3];
sx q[3];
rz(-1.6986596) q[3];
sx q[3];
rz(-1.2996775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
