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
rz(-2.0074453) q[1];
sx q[1];
rz(2.3448047) q[1];
sx q[1];
rz(12.264836) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8753531) q[0];
sx q[0];
rz(-1.7325028) q[0];
sx q[0];
rz(-1.7577359) q[0];
rz(1.8704988) q[2];
sx q[2];
rz(-1.7500008) q[2];
sx q[2];
rz(-1.9546997) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63362279) q[1];
sx q[1];
rz(-0.75532856) q[1];
sx q[1];
rz(1.4719187) q[1];
x q[2];
rz(1.9728679) q[3];
sx q[3];
rz(-1.744606) q[3];
sx q[3];
rz(-0.92648348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6887168) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-2.0067046) q[2];
rz(-0.12198837) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(3.0854935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720471) q[0];
sx q[0];
rz(-0.24709728) q[0];
sx q[0];
rz(0.0023512996) q[0];
rz(0.40070847) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(2.2854663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99994217) q[0];
sx q[0];
rz(-2.460833) q[0];
sx q[0];
rz(-1.1589684) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8555774) q[2];
sx q[2];
rz(-2.0692424) q[2];
sx q[2];
rz(1.1666544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90710956) q[1];
sx q[1];
rz(-0.52965783) q[1];
sx q[1];
rz(-1.6698599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56613381) q[3];
sx q[3];
rz(-1.8834891) q[3];
sx q[3];
rz(-1.4110402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1043642) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(-2.6715703) q[2];
rz(-0.59703279) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(1.443583) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-2.6486588) q[0];
sx q[0];
rz(0.22802995) q[0];
rz(-2.6920964) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(-0.94625783) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0272574) q[0];
sx q[0];
rz(-1.9693621) q[0];
sx q[0];
rz(-0.27862676) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65785711) q[2];
sx q[2];
rz(-1.0491174) q[2];
sx q[2];
rz(1.1945832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5214649) q[1];
sx q[1];
rz(-1.9072272) q[1];
sx q[1];
rz(-1.7640339) q[1];
x q[2];
rz(-1.069293) q[3];
sx q[3];
rz(-0.32078241) q[3];
sx q[3];
rz(-0.71170413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9049282) q[2];
sx q[2];
rz(-1.4651352) q[2];
sx q[2];
rz(-1.6620592) q[2];
rz(-2.9605401) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(-2.2553867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(-0.89853483) q[0];
rz(0.50615519) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(-1.4599962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.579297) q[0];
sx q[0];
rz(-1.0387207) q[0];
sx q[0];
rz(1.4946372) q[0];
rz(-2.6674472) q[2];
sx q[2];
rz(-2.3398826) q[2];
sx q[2];
rz(2.2849883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0702563) q[1];
sx q[1];
rz(-1.500529) q[1];
sx q[1];
rz(-1.6855168) q[1];
x q[2];
rz(2.2257099) q[3];
sx q[3];
rz(-1.5228027) q[3];
sx q[3];
rz(-0.97150436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(1.8858887) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(-3.0953395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
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
rz(1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(2.0464121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9126606) q[0];
sx q[0];
rz(-1.4353509) q[0];
sx q[0];
rz(1.6333461) q[0];
rz(-pi) q[1];
rz(0.081549598) q[2];
sx q[2];
rz(-2.3429541) q[2];
sx q[2];
rz(0.35432409) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4597561) q[1];
sx q[1];
rz(-2.3533258) q[1];
sx q[1];
rz(2.9440109) q[1];
rz(-1.7322695) q[3];
sx q[3];
rz(-0.69500178) q[3];
sx q[3];
rz(-1.410786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3695662) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(-2.0728716) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0291979) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(-0.29119626) q[0];
rz(-1.683782) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(-0.88868946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55808961) q[0];
sx q[0];
rz(-2.0348685) q[0];
sx q[0];
rz(0.91832535) q[0];
rz(-pi) q[1];
rz(-3.1142919) q[2];
sx q[2];
rz(-0.75746398) q[2];
sx q[2];
rz(-2.2006019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1322358) q[1];
sx q[1];
rz(-2.0814544) q[1];
sx q[1];
rz(2.4978994) q[1];
x q[2];
rz(2.3360152) q[3];
sx q[3];
rz(-2.71851) q[3];
sx q[3];
rz(0.42832908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7500744) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(-2.1938426) q[2];
rz(0.1968955) q[3];
sx q[3];
rz(-2.9010549) q[3];
sx q[3];
rz(-2.8314364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8295558) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(0.46992508) q[0];
rz(1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(-1.4039111) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023177308) q[0];
sx q[0];
rz(-2.1274533) q[0];
sx q[0];
rz(-0.413229) q[0];
x q[1];
rz(1.092488) q[2];
sx q[2];
rz(-1.7564536) q[2];
sx q[2];
rz(-0.80362475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0216071) q[1];
sx q[1];
rz(-2.48263) q[1];
sx q[1];
rz(-0.95665415) q[1];
x q[2];
rz(-2.5581854) q[3];
sx q[3];
rz(-1.2800763) q[3];
sx q[3];
rz(-0.099778508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5002084) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(-3.0339962) q[2];
rz(-0.61947668) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(-1.7776325) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(0.25303823) q[0];
rz(-0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(0.94892445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8924361) q[0];
sx q[0];
rz(-1.2741486) q[0];
sx q[0];
rz(-0.15693024) q[0];
x q[1];
rz(-1.1429092) q[2];
sx q[2];
rz(-2.6978931) q[2];
sx q[2];
rz(-2.2744409) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6151877) q[1];
sx q[1];
rz(-1.5171836) q[1];
sx q[1];
rz(-2.0508906) q[1];
rz(-pi) q[2];
rz(2.6213245) q[3];
sx q[3];
rz(-2.0754077) q[3];
sx q[3];
rz(2.5443973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(-0.31069791) q[2];
rz(-0.24142309) q[3];
sx q[3];
rz(-1.2025236) q[3];
sx q[3];
rz(-2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(-1.9061506) q[0];
rz(0.48108092) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(1.7135886) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1989216) q[0];
sx q[0];
rz(-2.3301972) q[0];
sx q[0];
rz(-1.5862203) q[0];
rz(-0.22869208) q[2];
sx q[2];
rz(-1.8426101) q[2];
sx q[2];
rz(2.6120409) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9199099) q[1];
sx q[1];
rz(-0.35791985) q[1];
sx q[1];
rz(-0.81166761) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2441745) q[3];
sx q[3];
rz(-2.2664968) q[3];
sx q[3];
rz(0.43750924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9674176) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(-0.027912557) q[2];
rz(-2.2715955) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(1.1869259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025539909) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(1.0593587) q[0];
rz(-1.4147883) q[1];
sx q[1];
rz(-1.4780412) q[1];
sx q[1];
rz(-0.47763225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896509) q[0];
sx q[0];
rz(-0.7583337) q[0];
sx q[0];
rz(0.21679057) q[0];
rz(-pi) q[1];
rz(-1.9066493) q[2];
sx q[2];
rz(-1.5726461) q[2];
sx q[2];
rz(2.2950878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76948386) q[1];
sx q[1];
rz(-1.4585351) q[1];
sx q[1];
rz(3.0873507) q[1];
rz(-pi) q[2];
rz(0.60250256) q[3];
sx q[3];
rz(-0.91711603) q[3];
sx q[3];
rz(2.9127171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16964218) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(2.7726445) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(2.1784311) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030768) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(0.44454642) q[1];
sx q[1];
rz(-2.3311756) q[1];
sx q[1];
rz(-0.12745007) q[1];
rz(-2.5196471) q[2];
sx q[2];
rz(-2.6126886) q[2];
sx q[2];
rz(-2.9168203) q[2];
rz(-1.5515242) q[3];
sx q[3];
rz(-2.0382593) q[3];
sx q[3];
rz(-0.12893547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
