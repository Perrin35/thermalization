OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(0.62358207) q[0];
rz(3.4317598) q[1];
sx q[1];
rz(5.5640339) q[1];
sx q[1];
rz(13.066864) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5252285) q[0];
sx q[0];
rz(-2.541757) q[0];
sx q[0];
rz(-3.1347549) q[0];
rz(-1.0153158) q[2];
sx q[2];
rz(-2.9657288) q[2];
sx q[2];
rz(-1.6989087) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1653633) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(-0.43433365) q[1];
rz(-pi) q[2];
rz(-1.8960564) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(-1.1577275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(-1.0789385) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21214813) q[0];
sx q[0];
rz(-2.0245027) q[0];
sx q[0];
rz(-1.8871904) q[0];
rz(-pi) q[1];
rz(-2.6639054) q[2];
sx q[2];
rz(-0.3590695) q[2];
sx q[2];
rz(-1.8600841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7484819) q[1];
sx q[1];
rz(-1.0903653) q[1];
sx q[1];
rz(1.837681) q[1];
x q[2];
rz(-1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(2.1982847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(-0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(-1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.741515) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6436359) q[0];
sx q[0];
rz(-2.6640577) q[0];
sx q[0];
rz(0.60995539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4351575) q[2];
sx q[2];
rz(-1.7559768) q[2];
sx q[2];
rz(-0.26091012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53576614) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(-2.1893326) q[1];
rz(-pi) q[2];
rz(2.9754144) q[3];
sx q[3];
rz(-1.8718534) q[3];
sx q[3];
rz(-2.500246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.489958) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79384971) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(-0.34805527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30663438) q[0];
sx q[0];
rz(-1.8638896) q[0];
sx q[0];
rz(-0.65960633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9449171) q[2];
sx q[2];
rz(-1.7110363) q[2];
sx q[2];
rz(1.8082878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.066612331) q[1];
sx q[1];
rz(-2.1790824) q[1];
sx q[1];
rz(3.0327256) q[1];
rz(-pi) q[2];
rz(-0.16163687) q[3];
sx q[3];
rz(-1.703754) q[3];
sx q[3];
rz(1.1383575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-2.7187738) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(3.056934) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87930644) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(2.6111531) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.2984498) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12621524) q[0];
sx q[0];
rz(-1.7324289) q[0];
sx q[0];
rz(-0.14590185) q[0];
rz(-pi) q[1];
rz(-1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(-2.1538018) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1227222) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(0.36032569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1634444) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(-3.0654207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-2.6679664) q[2];
rz(-0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6376003) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(-0.87431327) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186819) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(-2.8055311) q[0];
rz(-pi) q[1];
rz(0.45852197) q[2];
sx q[2];
rz(-2.6715171) q[2];
sx q[2];
rz(2.2058861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14086831) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(2.9904757) q[1];
x q[2];
rz(-1.7344597) q[3];
sx q[3];
rz(-1.6657077) q[3];
sx q[3];
rz(-1.6223736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(-2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(-2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2433462) q[0];
sx q[0];
rz(-1.4833437) q[0];
sx q[0];
rz(1.5623708) q[0];
rz(1.2862455) q[2];
sx q[2];
rz(-1.3853405) q[2];
sx q[2];
rz(2.4539349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0737887) q[1];
sx q[1];
rz(-1.5859335) q[1];
sx q[1];
rz(-3.036036) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94957955) q[3];
sx q[3];
rz(-2.5887244) q[3];
sx q[3];
rz(-2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.9246624) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970012) q[0];
sx q[0];
rz(-2.3750711) q[0];
sx q[0];
rz(-2.9826829) q[0];
rz(-1.2411225) q[2];
sx q[2];
rz(-0.85061073) q[2];
sx q[2];
rz(2.5230797) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5319547) q[1];
sx q[1];
rz(-1.6913809) q[1];
sx q[1];
rz(3.0871255) q[1];
x q[2];
rz(-1.8955599) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(1.8079545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(2.111964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(-1.2458941) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-0.54661173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194191) q[0];
sx q[0];
rz(-1.1328567) q[0];
sx q[0];
rz(2.7287448) q[0];
rz(-pi) q[1];
rz(-1.1999646) q[2];
sx q[2];
rz(-1.2002581) q[2];
sx q[2];
rz(-1.3818936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6555772) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(-1.1737215) q[1];
x q[2];
rz(-2.6528477) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(-2.3084156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0299915) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.4477504) q[2];
rz(1.1374121) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428064) q[0];
sx q[0];
rz(-1.2278623) q[0];
sx q[0];
rz(2.2862611) q[0];
rz(-pi) q[1];
rz(-2.131358) q[2];
sx q[2];
rz(-2.8786504) q[2];
sx q[2];
rz(1.5514785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.087079436) q[1];
sx q[1];
rz(-0.88216773) q[1];
sx q[1];
rz(0.41333945) q[1];
x q[2];
rz(-2.7865949) q[3];
sx q[3];
rz(-2.1600351) q[3];
sx q[3];
rz(2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(-2.2383402) q[2];
rz(1.5385657) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(0.8159591) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(-2.3564561) q[2];
sx q[2];
rz(-2.5265836) q[2];
sx q[2];
rz(2.1267736) q[2];
rz(-1.611931) q[3];
sx q[3];
rz(-1.6251246) q[3];
sx q[3];
rz(2.3390935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];