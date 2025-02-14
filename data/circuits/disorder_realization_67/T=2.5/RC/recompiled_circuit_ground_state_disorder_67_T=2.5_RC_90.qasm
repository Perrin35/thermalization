OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(-1.1986873) q[0];
rz(1.1341473) q[1];
sx q[1];
rz(-2.3448047) q[1];
sx q[1];
rz(-0.30153433) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59936675) q[0];
sx q[0];
rz(-2.8950373) q[0];
sx q[0];
rz(0.85041468) q[0];
rz(-pi) q[1];
rz(2.1211336) q[2];
sx q[2];
rz(-2.7937903) q[2];
sx q[2];
rz(2.2344799) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49820013) q[1];
sx q[1];
rz(-2.3215356) q[1];
sx q[1];
rz(-0.092685164) q[1];
rz(-1.1490378) q[3];
sx q[3];
rz(-0.43614498) q[3];
sx q[3];
rz(2.8835322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4528759) q[2];
sx q[2];
rz(-1.0998925) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(0.12198837) q[3];
sx q[3];
rz(-0.9361836) q[3];
sx q[3];
rz(-0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695456) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(-0.0023512996) q[0];
rz(-0.40070847) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(-2.2854663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99994217) q[0];
sx q[0];
rz(-2.460833) q[0];
sx q[0];
rz(1.9826243) q[0];
rz(-pi) q[1];
rz(1.0547423) q[2];
sx q[2];
rz(-1.3203838) q[2];
sx q[2];
rz(2.5977787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5634671) q[1];
sx q[1];
rz(-1.620786) q[1];
sx q[1];
rz(-1.0432788) q[1];
rz(-pi) q[2];
rz(-1.2049742) q[3];
sx q[3];
rz(-1.0351887) q[3];
sx q[3];
rz(-2.7887695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1043642) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(2.6715703) q[2];
rz(-0.59703279) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(-1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715348) q[0];
sx q[0];
rz(-2.6486588) q[0];
sx q[0];
rz(2.9135627) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(2.1953348) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957009) q[0];
sx q[0];
rz(-1.8270703) q[0];
sx q[0];
rz(-1.9836224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1990132) q[2];
sx q[2];
rz(-1.0120856) q[2];
sx q[2];
rz(-3.1329346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98639224) q[1];
sx q[1];
rz(-2.7554535) q[1];
sx q[1];
rz(2.6394352) q[1];
x q[2];
rz(-0.15839307) q[3];
sx q[3];
rz(-1.8509281) q[3];
sx q[3];
rz(-1.9060077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2366644) q[2];
sx q[2];
rz(-1.4651352) q[2];
sx q[2];
rz(-1.6620592) q[2];
rz(2.9605401) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(2.2553867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7774696) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(-2.2430578) q[0];
rz(-2.6354375) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047193479) q[0];
sx q[0];
rz(-1.5051821) q[0];
sx q[0];
rz(2.608247) q[0];
rz(-pi) q[1];
rz(-2.6674472) q[2];
sx q[2];
rz(-0.80171004) q[2];
sx q[2];
rz(-2.2849883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0713363) q[1];
sx q[1];
rz(-1.500529) q[1];
sx q[1];
rz(1.6855168) q[1];
x q[2];
rz(0.91588275) q[3];
sx q[3];
rz(-1.5228027) q[3];
sx q[3];
rz(-2.1700883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(-1.6391222) q[2];
rz(1.255704) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(-0.046253117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2332377) q[0];
sx q[0];
rz(-0.73035705) q[0];
sx q[0];
rz(0.18049845) q[0];
rz(-1.3795229) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(-1.0951805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9126606) q[0];
sx q[0];
rz(-1.4353509) q[0];
sx q[0];
rz(1.6333461) q[0];
rz(-pi) q[1];
rz(1.6542475) q[2];
sx q[2];
rz(-2.3660283) q[2];
sx q[2];
rz(2.6706539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68183652) q[1];
sx q[1];
rz(-2.3533258) q[1];
sx q[1];
rz(-2.9440109) q[1];
rz(-pi) q[2];
rz(0.88222701) q[3];
sx q[3];
rz(-1.6739356) q[3];
sx q[3];
rz(-0.28447275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77202648) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(2.0728716) q[2];
rz(-0.42701834) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0291979) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(0.29119626) q[0];
rz(1.4578106) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(2.2529032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4584258) q[0];
sx q[0];
rz(-0.99686503) q[0];
sx q[0];
rz(-2.5794585) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3843147) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(0.64964408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1322358) q[1];
sx q[1];
rz(-1.0601383) q[1];
sx q[1];
rz(-2.4978994) q[1];
x q[2];
rz(0.8055775) q[3];
sx q[3];
rz(-2.71851) q[3];
sx q[3];
rz(2.7132636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39151829) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(-2.1938426) q[2];
rz(-0.1968955) q[3];
sx q[3];
rz(-2.9010549) q[3];
sx q[3];
rz(-0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8295558) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(-0.46992508) q[0];
rz(-1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(-1.7376815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.9187154) q[0];
sx q[0];
rz(2.1676284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.092488) q[2];
sx q[2];
rz(-1.7564536) q[2];
sx q[2];
rz(-2.3379679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95937377) q[1];
sx q[1];
rz(-1.2101908) q[1];
sx q[1];
rz(-2.1350767) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9150195) q[3];
sx q[3];
rz(-1.0148409) q[3];
sx q[3];
rz(-1.8575625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64138428) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(-3.0339962) q[2];
rz(2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-2.0893607) q[0];
sx q[0];
rz(-0.25303823) q[0];
rz(0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(2.1926682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.388577) q[0];
sx q[0];
rz(-2.8070794) q[0];
sx q[0];
rz(-2.0434414) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19472943) q[2];
sx q[2];
rz(-1.9720594) q[2];
sx q[2];
rz(-1.8068318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1250985) q[1];
sx q[1];
rz(-2.0501425) q[1];
sx q[1];
rz(0.060430364) q[1];
rz(-pi) q[2];
rz(0.83794318) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(1.4670682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-0.88699114) q[2];
sx q[2];
rz(2.8308947) q[2];
rz(0.24142309) q[3];
sx q[3];
rz(-1.2025236) q[3];
sx q[3];
rz(2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0017515) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(1.9061506) q[0];
rz(2.6605117) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.7135886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9426711) q[0];
sx q[0];
rz(-2.3301972) q[0];
sx q[0];
rz(1.5553724) q[0];
rz(-pi) q[1];
rz(2.253654) q[2];
sx q[2];
rz(-2.7882034) q[2];
sx q[2];
rz(-1.2438174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.075632) q[1];
sx q[1];
rz(-1.8143225) q[1];
sx q[1];
rz(-1.3058409) q[1];
x q[2];
rz(1.2441745) q[3];
sx q[3];
rz(-0.87509586) q[3];
sx q[3];
rz(2.7040834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9674176) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(3.1136801) q[2];
rz(2.2715955) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1160527) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(1.0593587) q[0];
rz(1.7268044) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(-2.6639604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896509) q[0];
sx q[0];
rz(-2.383259) q[0];
sx q[0];
rz(0.21679057) q[0];
rz(-1.2349434) q[2];
sx q[2];
rz(-1.5726461) q[2];
sx q[2];
rz(-2.2950878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2207909) q[1];
sx q[1];
rz(-3.0169636) q[1];
sx q[1];
rz(1.1225379) q[1];
rz(0.82179339) q[3];
sx q[3];
rz(-2.0374277) q[3];
sx q[3];
rz(2.1958283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9719505) q[2];
sx q[2];
rz(-0.98196882) q[2];
sx q[2];
rz(0.36894813) q[2];
rz(1.87489) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(-2.1784311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0030768) q[0];
sx q[0];
rz(-1.2110447) q[0];
sx q[0];
rz(-1.3883653) q[0];
rz(-0.44454642) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(-2.5196471) q[2];
sx q[2];
rz(-2.6126886) q[2];
sx q[2];
rz(-2.9168203) q[2];
rz(-0.46753771) q[3];
sx q[3];
rz(-1.5880006) q[3];
sx q[3];
rz(1.4505462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
