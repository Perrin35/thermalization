OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4386908) q[0];
sx q[0];
rz(-0.54649788) q[0];
sx q[0];
rz(1.7706857) q[0];
rz(-pi) q[1];
rz(-1.8194524) q[2];
sx q[2];
rz(-2.6373632) q[2];
sx q[2];
rz(-2.3766975) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6403113) q[1];
sx q[1];
rz(-2.8456563) q[1];
sx q[1];
rz(-2.179115) q[1];
rz(-pi) q[2];
rz(-1.3931307) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(-0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(0.5509848) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-0.93634161) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124356) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-2.9108414) q[0];
rz(0.78511946) q[2];
sx q[2];
rz(-1.8765419) q[2];
sx q[2];
rz(2.608125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7649819) q[1];
sx q[1];
rz(-2.2602343) q[1];
sx q[1];
rz(-1.7060243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.068088) q[3];
sx q[3];
rz(-1.5503251) q[3];
sx q[3];
rz(2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-0.42638391) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8957829) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(0.93908969) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1547326) q[0];
sx q[0];
rz(-1.2504471) q[0];
sx q[0];
rz(-0.31517584) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5005641) q[2];
sx q[2];
rz(-0.68968455) q[2];
sx q[2];
rz(-0.94674142) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8507104) q[1];
sx q[1];
rz(-1.1516654) q[1];
sx q[1];
rz(0.23371975) q[1];
x q[2];
rz(-2.5251758) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(-2.4544231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(0.38763186) q[3];
sx q[3];
rz(-1.5165611) q[3];
sx q[3];
rz(0.38813996) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664292) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(2.373383) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216953) q[0];
sx q[0];
rz(-1.2384402) q[0];
sx q[0];
rz(2.7601348) q[0];
rz(-1.5966162) q[2];
sx q[2];
rz(-0.46586793) q[2];
sx q[2];
rz(-2.4398838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0087352) q[1];
sx q[1];
rz(-2.1454304) q[1];
sx q[1];
rz(0.43032129) q[1];
rz(-pi) q[2];
rz(-2.4079011) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(1.2804077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(-2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30848635) q[0];
sx q[0];
rz(-1.3163438) q[0];
sx q[0];
rz(1.893505) q[0];
rz(-2.167335) q[2];
sx q[2];
rz(-1.1968808) q[2];
sx q[2];
rz(-2.8982382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1041303) q[1];
sx q[1];
rz(-2.3298652) q[1];
sx q[1];
rz(-2.3292259) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1220783) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(0.36908484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(-1.5340012) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896478) q[0];
sx q[0];
rz(-1.6501353) q[0];
sx q[0];
rz(0.21110714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8842823) q[2];
sx q[2];
rz(-0.67337155) q[2];
sx q[2];
rz(1.1416669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51965442) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-2.7154891) q[1];
rz(-1.9147647) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(-0.96836585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.9645875) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(2.0578407) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18400684) q[0];
sx q[0];
rz(-1.1054966) q[0];
sx q[0];
rz(0.9853978) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8132642) q[2];
sx q[2];
rz(-2.7331181) q[2];
sx q[2];
rz(2.7834053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5092897) q[1];
sx q[1];
rz(-1.0273233) q[1];
sx q[1];
rz(1.6093045) q[1];
x q[2];
rz(-0.13109644) q[3];
sx q[3];
rz(-0.56250611) q[3];
sx q[3];
rz(1.8545811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5291418) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(1.051349) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(2.8578551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34127125) q[0];
sx q[0];
rz(-0.30297908) q[0];
sx q[0];
rz(-3.0269701) q[0];
x q[1];
rz(2.9494638) q[2];
sx q[2];
rz(-1.2461975) q[2];
sx q[2];
rz(-2.0786376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6730496) q[1];
sx q[1];
rz(-2.0646411) q[1];
sx q[1];
rz(1.3480575) q[1];
x q[2];
rz(0.92054263) q[3];
sx q[3];
rz(-1.0245819) q[3];
sx q[3];
rz(1.2342681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.4452176) q[2];
rz(1.5971659) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504869) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42800316) q[0];
sx q[0];
rz(-2.028095) q[0];
sx q[0];
rz(2.2767115) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32585085) q[2];
sx q[2];
rz(-1.4928276) q[2];
sx q[2];
rz(0.98380145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20977565) q[1];
sx q[1];
rz(-1.6546384) q[1];
sx q[1];
rz(-0.2340338) q[1];
x q[2];
rz(2.2194355) q[3];
sx q[3];
rz(-1.4033917) q[3];
sx q[3];
rz(0.64099121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(-1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-2.8737601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7627015) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(2.1428109) q[0];
rz(-pi) q[1];
rz(-1.7540625) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(-2.0545517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2512868) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(-0.73015405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1166499) q[3];
sx q[3];
rz(-2.5513259) q[3];
sx q[3];
rz(2.2772307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(-0.5464856) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(0.71074625) q[2];
sx q[2];
rz(-1.1087316) q[2];
sx q[2];
rz(1.0089594) q[2];
rz(0.6775425) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
