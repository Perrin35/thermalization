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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7843648) q[0];
sx q[0];
rz(-0.89656943) q[0];
sx q[0];
rz(-2.7927464) q[0];
rz(1.8741605) q[2];
sx q[2];
rz(-0.99736664) q[2];
sx q[2];
rz(1.0614741) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6837343) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(-3.0617759) q[1];
x q[2];
rz(-0.47655063) q[3];
sx q[3];
rz(-0.50560564) q[3];
sx q[3];
rz(1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(-1.7976044) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3450736) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(-2.2497798) q[0];
rz(2.0925518) q[2];
sx q[2];
rz(-0.21460303) q[2];
sx q[2];
rz(1.3133446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5492925) q[1];
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
x q[1];
rz(0.0011757294) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-0.95138335) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(-1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.552859) q[0];
sx q[0];
rz(3.0271295) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(2.9797629) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4672218) q[0];
sx q[0];
rz(-0.97839744) q[0];
sx q[0];
rz(0.47426275) q[0];
x q[1];
rz(1.4828909) q[2];
sx q[2];
rz(-1.1832596) q[2];
sx q[2];
rz(-1.0532925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1766277) q[1];
sx q[1];
rz(-1.8102268) q[1];
sx q[1];
rz(-0.62271948) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1597685) q[3];
sx q[3];
rz(-2.3279802) q[3];
sx q[3];
rz(-0.66853722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-3.0008345) q[2];
rz(-0.55177871) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(0.75132918) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(0.67489135) q[0];
rz(-0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-0.45796576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602011) q[0];
sx q[0];
rz(-1.0835674) q[0];
sx q[0];
rz(-1.1283406) q[0];
x q[1];
rz(2.0003597) q[2];
sx q[2];
rz(-1.8515966) q[2];
sx q[2];
rz(-2.8605638) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47487101) q[1];
sx q[1];
rz(-0.41731605) q[1];
sx q[1];
rz(-0.46837193) q[1];
x q[2];
rz(-2.6334718) q[3];
sx q[3];
rz(-1.4773989) q[3];
sx q[3];
rz(1.7600218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-2.8221455) q[2];
rz(-1.8114629) q[3];
sx q[3];
rz(-1.0165241) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.2215479) q[0];
rz(2.9087861) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(0.98308841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9778091) q[0];
sx q[0];
rz(-1.0484694) q[0];
sx q[0];
rz(0.99951007) q[0];
x q[1];
rz(-1.4120502) q[2];
sx q[2];
rz(-0.6266784) q[2];
sx q[2];
rz(1.9959297) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5447588) q[1];
sx q[1];
rz(-1.4248409) q[1];
sx q[1];
rz(2.9533141) q[1];
rz(-pi) q[2];
rz(0.1666937) q[3];
sx q[3];
rz(-1.1536666) q[3];
sx q[3];
rz(0.47458173) q[3];
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
rz(2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.7406834) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(0.54814235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948471) q[0];
sx q[0];
rz(-1.7821494) q[0];
sx q[0];
rz(-2.4708492) q[0];
x q[1];
rz(0.55193837) q[2];
sx q[2];
rz(-1.1480398) q[2];
sx q[2];
rz(-1.1814337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9029451) q[1];
sx q[1];
rz(-2.3460707) q[1];
sx q[1];
rz(-2.5832157) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75597922) q[3];
sx q[3];
rz(-1.090831) q[3];
sx q[3];
rz(0.85473138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6846201) q[2];
sx q[2];
rz(-0.62070864) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.5439344) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-0.4959929) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(1.1167663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3163213) q[0];
sx q[0];
rz(-0.50646913) q[0];
sx q[0];
rz(-0.78790085) q[0];
x q[1];
rz(-1.1363813) q[2];
sx q[2];
rz(-2.5602753) q[2];
sx q[2];
rz(-0.28806799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6903444) q[1];
sx q[1];
rz(-2.6670229) q[1];
sx q[1];
rz(-2.3869273) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6326866) q[3];
sx q[3];
rz(-1.1155432) q[3];
sx q[3];
rz(0.52857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3108959) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(-0.35487077) q[2];
rz(-0.76534671) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(-3.0520458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(-0.76559693) q[0];
rz(0.87144026) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.8188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8894371) q[0];
sx q[0];
rz(-1.3103292) q[0];
sx q[0];
rz(1.6244662) q[0];
rz(-pi) q[1];
rz(-2.2501588) q[2];
sx q[2];
rz(-1.3855033) q[2];
sx q[2];
rz(-2.7012555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8812528) q[1];
sx q[1];
rz(-1.5543206) q[1];
sx q[1];
rz(-0.50838281) q[1];
rz(-1.345383) q[3];
sx q[3];
rz(-1.3691878) q[3];
sx q[3];
rz(0.20568337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(-2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-0.27716032) q[0];
rz(2.3932338) q[1];
sx q[1];
rz(-0.44735083) q[1];
sx q[1];
rz(-1.998273) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80625929) q[0];
sx q[0];
rz(-1.2657832) q[0];
sx q[0];
rz(0.2103896) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2758946) q[2];
sx q[2];
rz(-1.7324466) q[2];
sx q[2];
rz(2.9557712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9710247) q[1];
sx q[1];
rz(-1.8877561) q[1];
sx q[1];
rz(-0.16522878) q[1];
rz(-0.36181776) q[3];
sx q[3];
rz(-1.731195) q[3];
sx q[3];
rz(-3.0830864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12298909) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(-0.019512026) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(1.5787517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39017856) q[0];
sx q[0];
rz(-1.7008243) q[0];
sx q[0];
rz(3.0517654) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8291924) q[2];
sx q[2];
rz(-0.62587291) q[2];
sx q[2];
rz(-1.3385119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.409118) q[1];
sx q[1];
rz(-2.1555087) q[1];
sx q[1];
rz(-0.52701108) q[1];
rz(-2.9928777) q[3];
sx q[3];
rz(-0.84245719) q[3];
sx q[3];
rz(-3.0535798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7623357) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
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
rz(1.3723102) q[3];
sx q[3];
rz(-2.1673421) q[3];
sx q[3];
rz(2.0094595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
