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
rz(-2.969279) q[0];
sx q[0];
rz(-0.20792374) q[0];
sx q[0];
rz(-1.2907668) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(-0.42833498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9865532) q[0];
sx q[0];
rz(-1.4615834) q[0];
sx q[0];
rz(0.10515736) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7784155) q[2];
sx q[2];
rz(-0.70549772) q[2];
sx q[2];
rz(0.35340912) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2178035) q[1];
sx q[1];
rz(-2.3084967) q[1];
sx q[1];
rz(-1.6290725) q[1];
rz(-0.74083565) q[3];
sx q[3];
rz(-0.29262283) q[3];
sx q[3];
rz(0.13340852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20994818) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(1.2098562) q[2];
rz(3.0620388) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2570268) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(-2.7575745) q[0];
rz(2.2513023) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(0.30337897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7913032) q[0];
sx q[0];
rz(-1.3801738) q[0];
sx q[0];
rz(1.0269357) q[0];
x q[1];
rz(-2.882134) q[2];
sx q[2];
rz(-2.0739809) q[2];
sx q[2];
rz(0.34215701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0133923) q[1];
sx q[1];
rz(-1.6329995) q[1];
sx q[1];
rz(1.661646) q[1];
rz(-pi) q[2];
rz(-0.65360131) q[3];
sx q[3];
rz(-1.2051688) q[3];
sx q[3];
rz(-1.0753683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7552135) q[2];
sx q[2];
rz(-1.1073802) q[2];
sx q[2];
rz(2.4066822) q[2];
rz(-2.2928061) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-2.7809704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81247771) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(-3.062881) q[0];
rz(-2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(1.8471921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98662739) q[0];
sx q[0];
rz(-1.5091578) q[0];
sx q[0];
rz(0.060527965) q[0];
rz(-pi) q[1];
rz(0.90581494) q[2];
sx q[2];
rz(-0.73630263) q[2];
sx q[2];
rz(-1.6396963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6523167) q[1];
sx q[1];
rz(-1.4536263) q[1];
sx q[1];
rz(0.28276171) q[1];
x q[2];
rz(-2.2749316) q[3];
sx q[3];
rz(-1.3876283) q[3];
sx q[3];
rz(-2.4742692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2417629) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(-1.5166616) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4888332) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(-2.5732727) q[0];
rz(-1.9991416) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.4521339) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947493) q[0];
sx q[0];
rz(-2.1921625) q[0];
sx q[0];
rz(-1.2970379) q[0];
rz(-1.4771039) q[2];
sx q[2];
rz(-1.7498651) q[2];
sx q[2];
rz(2.8833517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4751079) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(0.38352769) q[1];
x q[2];
rz(-2.9630757) q[3];
sx q[3];
rz(-1.9048637) q[3];
sx q[3];
rz(2.8802383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(-3.0647035) q[2];
rz(-2.7422089) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15453108) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(-2.8416908) q[0];
rz(-1.1497644) q[1];
sx q[1];
rz(-2.1297784) q[1];
sx q[1];
rz(-0.48847517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91575275) q[0];
sx q[0];
rz(-1.7560648) q[0];
sx q[0];
rz(-0.033323296) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68814338) q[2];
sx q[2];
rz(-0.71758413) q[2];
sx q[2];
rz(0.46844278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3255182) q[1];
sx q[1];
rz(-1.5518922) q[1];
sx q[1];
rz(-2.5579648) q[1];
rz(2.4998922) q[3];
sx q[3];
rz(-1.3535168) q[3];
sx q[3];
rz(0.02441306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(-1.4786973) q[2];
rz(1.1100769) q[3];
sx q[3];
rz(-0.9980945) q[3];
sx q[3];
rz(-0.64322513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.4019302) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(-0.31841835) q[0];
rz(-3.1069801) q[1];
sx q[1];
rz(-0.56762677) q[1];
sx q[1];
rz(-0.45077032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0552669) q[0];
sx q[0];
rz(-1.0604211) q[0];
sx q[0];
rz(0.016252131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8209967) q[2];
sx q[2];
rz(-2.2911173) q[2];
sx q[2];
rz(2.3483417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83799441) q[1];
sx q[1];
rz(-0.7571836) q[1];
sx q[1];
rz(2.1794469) q[1];
rz(-pi) q[2];
rz(-0.81574635) q[3];
sx q[3];
rz(-1.8459847) q[3];
sx q[3];
rz(-1.311917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0938809) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(-0.06165687) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(2.4390167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75673574) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(0.039948832) q[0];
rz(2.3710251) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-2.9133453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4994482) q[0];
sx q[0];
rz(-0.11760437) q[0];
sx q[0];
rz(0.85310491) q[0];
rz(-pi) q[1];
rz(3.1399916) q[2];
sx q[2];
rz(-1.8769771) q[2];
sx q[2];
rz(1.9322576) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3747663) q[1];
sx q[1];
rz(-3.0438381) q[1];
sx q[1];
rz(-0.575211) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93668117) q[3];
sx q[3];
rz(-1.2290658) q[3];
sx q[3];
rz(3.0724728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.063529) q[2];
sx q[2];
rz(-1.0588131) q[2];
sx q[2];
rz(-0.25827363) q[2];
rz(-2.9410948) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(0.18331461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(2.8343416) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(-2.5603851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191588) q[0];
sx q[0];
rz(-0.40589505) q[0];
sx q[0];
rz(-1.4916363) q[0];
rz(0.70971428) q[2];
sx q[2];
rz(-1.6953371) q[2];
sx q[2];
rz(2.1309851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4631066) q[1];
sx q[1];
rz(-2.3358279) q[1];
sx q[1];
rz(-1.8496978) q[1];
x q[2];
rz(0.56116207) q[3];
sx q[3];
rz(-0.61120874) q[3];
sx q[3];
rz(-1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4622978) q[2];
sx q[2];
rz(-1.1335979) q[2];
sx q[2];
rz(-3.1104258) q[2];
rz(-0.22552414) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(2.1112736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088260055) q[0];
sx q[0];
rz(-0.044476155) q[0];
sx q[0];
rz(-2.4326676) q[0];
rz(2.9290579) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(-2.2841891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5851201) q[0];
sx q[0];
rz(-1.5657288) q[0];
sx q[0];
rz(1.4410254) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8959664) q[2];
sx q[2];
rz(-1.0602927) q[2];
sx q[2];
rz(0.73869642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6826806) q[1];
sx q[1];
rz(-0.68411982) q[1];
sx q[1];
rz(-0.64351179) q[1];
rz(1.3579426) q[3];
sx q[3];
rz(-2.2567333) q[3];
sx q[3];
rz(2.6712772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7015486) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(-0.80553833) q[2];
rz(2.7696179) q[3];
sx q[3];
rz(-1.6476846) q[3];
sx q[3];
rz(2.7443547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032967903) q[0];
sx q[0];
rz(-1.911835) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(-2.4188614) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(-2.5236948) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5586097) q[0];
sx q[0];
rz(-1.8454058) q[0];
sx q[0];
rz(0.78297575) q[0];
rz(-0.63510908) q[2];
sx q[2];
rz(-0.76580566) q[2];
sx q[2];
rz(2.2857411) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7405211) q[1];
sx q[1];
rz(-0.24722543) q[1];
sx q[1];
rz(2.3257491) q[1];
rz(1.4862107) q[3];
sx q[3];
rz(-1.4874377) q[3];
sx q[3];
rz(-0.74970923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.939398) q[2];
sx q[2];
rz(-1.7859744) q[2];
sx q[2];
rz(-0.00016577684) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(0.59529006) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(-1.3407002) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(0.98967057) q[2];
sx q[2];
rz(-2.4151617) q[2];
sx q[2];
rz(-2.5834609) q[2];
rz(0.44092785) q[3];
sx q[3];
rz(-1.8971414) q[3];
sx q[3];
rz(-1.8446445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
