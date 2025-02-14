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
rz(2.607333) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7843648) q[0];
sx q[0];
rz(-0.89656943) q[0];
sx q[0];
rz(-2.7927464) q[0];
x q[1];
rz(0.59492971) q[2];
sx q[2];
rz(-1.317136) q[2];
sx q[2];
rz(-2.4640535) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45785832) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(-3.0617759) q[1];
rz(-pi) q[2];
rz(2.684428) q[3];
sx q[3];
rz(-1.7948397) q[3];
sx q[3];
rz(-0.30686298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(1.7976044) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-0.070076076) q[3];
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
rz(2.6473815) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-0.31165037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3450736) q[0];
sx q[0];
rz(-2.2537044) q[0];
sx q[0];
rz(0.89181283) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0925518) q[2];
sx q[2];
rz(-0.21460303) q[2];
sx q[2];
rz(1.8282481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6271882) q[1];
sx q[1];
rz(-1.3986821) q[1];
sx q[1];
rz(1.8132416) q[1];
rz(-pi) q[2];
rz(-2.8007224) q[3];
sx q[3];
rz(-2.1580527) q[3];
sx q[3];
rz(-2.1275008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-2.1902093) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15762873) q[0];
sx q[0];
rz(-1.5887337) q[0];
sx q[0];
rz(0.11446318) q[0];
rz(2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(-2.9797629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4672218) q[0];
sx q[0];
rz(-2.1631952) q[0];
sx q[0];
rz(2.6673299) q[0];
rz(0.21185565) q[2];
sx q[2];
rz(-0.39688928) q[2];
sx q[2];
rz(-0.8241764) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.964965) q[1];
sx q[1];
rz(-1.3313659) q[1];
sx q[1];
rz(0.62271948) q[1];
x q[2];
rz(1.9818241) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(0.66853722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3942922) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(-3.0008345) q[2];
rz(-2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(2.4667013) q[0];
rz(2.9871125) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(-2.6836269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7688528) q[0];
sx q[0];
rz(-2.495666) q[0];
sx q[0];
rz(0.67966184) q[0];
x q[1];
rz(2.0003597) q[2];
sx q[2];
rz(-1.2899961) q[2];
sx q[2];
rz(-0.28102885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1110036) q[1];
sx q[1];
rz(-1.9408424) q[1];
sx q[1];
rz(1.3732598) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4639888) q[3];
sx q[3];
rz(-1.065101) q[3];
sx q[3];
rz(3.0042574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(2.8221455) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.2215479) q[0];
rz(0.23280652) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(-2.1585042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096689396) q[0];
sx q[0];
rz(-1.0830729) q[0];
sx q[0];
rz(-2.5414667) q[0];
rz(-pi) q[1];
rz(2.1914761) q[2];
sx q[2];
rz(-1.663637) q[2];
sx q[2];
rz(-0.29618057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0016664) q[1];
sx q[1];
rz(-1.3845433) q[1];
sx q[1];
rz(-1.4222533) q[1];
x q[2];
rz(1.9290673) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(0.080953065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38272196) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(2.0659857) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7406834) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(-1.4373454) q[0];
rz(-3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(2.5934503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9288612) q[0];
sx q[0];
rz(-0.91758801) q[0];
sx q[0];
rz(-1.3034588) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4324722) q[2];
sx q[2];
rz(-0.68163423) q[2];
sx q[2];
rz(-2.1644985) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74433148) q[1];
sx q[1];
rz(-1.1827227) q[1];
sx q[1];
rz(0.71340386) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75597922) q[3];
sx q[3];
rz(-2.0507617) q[3];
sx q[3];
rz(0.85473138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4569725) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-0.4959929) q[0];
rz(-1.687382) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(-1.1167663) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1710715) q[0];
sx q[0];
rz(-1.2215893) q[0];
sx q[0];
rz(1.1961351) q[0];
x q[1];
rz(-1.1363813) q[2];
sx q[2];
rz(-0.58131733) q[2];
sx q[2];
rz(0.28806799) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8161113) q[1];
sx q[1];
rz(-1.8891836) q[1];
sx q[1];
rz(0.35811425) q[1];
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
rz(-0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(0.35487077) q[2];
rz(-0.76534671) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68882051) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.3227468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2521556) q[0];
sx q[0];
rz(-1.8312635) q[0];
sx q[0];
rz(-1.5171264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2501588) q[2];
sx q[2];
rz(-1.3855033) q[2];
sx q[2];
rz(-2.7012555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2603399) q[1];
sx q[1];
rz(-1.5543206) q[1];
sx q[1];
rz(2.6332098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9348974) q[3];
sx q[3];
rz(-1.7915676) q[3];
sx q[3];
rz(-1.7305935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(-1.0929557) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073864989) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(0.27716032) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(1.1433196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(1.2593377) q[0];
x q[1];
rz(1.3242993) q[2];
sx q[2];
rz(-2.4213104) q[2];
sx q[2];
rz(-1.5697073) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7932897) q[1];
sx q[1];
rz(-1.4138681) q[1];
sx q[1];
rz(1.2497529) q[1];
rz(-2.7797749) q[3];
sx q[3];
rz(-1.731195) q[3];
sx q[3];
rz(3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2271759) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(-0.048967036) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.12298909) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(1.5787517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922958) q[0];
sx q[0];
rz(-1.6598633) q[0];
sx q[0];
rz(1.440248) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8291924) q[2];
sx q[2];
rz(-0.62587291) q[2];
sx q[2];
rz(-1.3385119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0622618) q[1];
sx q[1];
rz(-0.76594662) q[1];
sx q[1];
rz(0.92103157) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83694038) q[3];
sx q[3];
rz(-1.4599953) q[3];
sx q[3];
rz(1.5821804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(2.9760402) q[2];
rz(-0.60756573) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(1.4277369) q[2];
sx q[2];
rz(-2.3410627) q[2];
sx q[2];
rz(-1.9973199) q[2];
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
