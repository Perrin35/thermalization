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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9865532) q[0];
sx q[0];
rz(-1.4615834) q[0];
sx q[0];
rz(0.10515736) q[0];
x q[1];
rz(2.7784155) q[2];
sx q[2];
rz(-0.70549772) q[2];
sx q[2];
rz(-0.35340912) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0103156) q[1];
sx q[1];
rz(-2.4020264) q[1];
sx q[1];
rz(-0.063994813) q[1];
x q[2];
rz(1.3702014) q[3];
sx q[3];
rz(-1.7852968) q[3];
sx q[3];
rz(0.89591276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9316445) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.9317365) q[2];
rz(3.0620388) q[3];
sx q[3];
rz(-0.39593655) q[3];
sx q[3];
rz(2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88456589) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(0.89029038) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(2.8382137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.035159) q[0];
sx q[0];
rz(-2.1037397) q[0];
sx q[0];
rz(-2.919801) q[0];
x q[1];
rz(2.0069507) q[2];
sx q[2];
rz(-2.580603) q[2];
sx q[2];
rz(2.2962388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0413826) q[1];
sx q[1];
rz(-3.0315371) q[1];
sx q[1];
rz(-0.96918126) q[1];
rz(-2.0201265) q[3];
sx q[3];
rz(-0.96689618) q[3];
sx q[3];
rz(-0.7627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7552135) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(0.73491043) q[2];
rz(-0.84878659) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-0.36062226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(-1.2944006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20939182) q[0];
sx q[0];
rz(-3.0552312) q[0];
sx q[0];
rz(2.3461653) q[0];
rz(2.2357777) q[2];
sx q[2];
rz(-0.73630263) q[2];
sx q[2];
rz(-1.5018963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.489276) q[1];
sx q[1];
rz(-1.4536263) q[1];
sx q[1];
rz(2.8588309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9031676) q[3];
sx q[3];
rz(-0.88078558) q[3];
sx q[3];
rz(2.391614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89982975) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(-1.5166616) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(-0.70037705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65275943) q[0];
sx q[0];
rz(-0.82010287) q[0];
sx q[0];
rz(-2.5732727) q[0];
rz(1.9991416) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(1.4521339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88599151) q[0];
sx q[0];
rz(-1.7924249) q[0];
sx q[0];
rz(0.63943531) q[0];
rz(-1.6644888) q[2];
sx q[2];
rz(-1.7498651) q[2];
sx q[2];
rz(0.25824091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85739112) q[1];
sx q[1];
rz(-2.439154) q[1];
sx q[1];
rz(-2.0676916) q[1];
rz(-pi) q[2];
rz(-1.9098572) q[3];
sx q[3];
rz(-1.7393469) q[3];
sx q[3];
rz(1.3685365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6292754) q[2];
sx q[2];
rz(-0.53048152) q[2];
sx q[2];
rz(-0.076889195) q[2];
rz(0.39938375) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870616) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(-2.8416908) q[0];
rz(-1.9918282) q[1];
sx q[1];
rz(-2.1297784) q[1];
sx q[1];
rz(-2.6531175) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91575275) q[0];
sx q[0];
rz(-1.3855278) q[0];
sx q[0];
rz(-0.033323296) q[0];
rz(-pi) q[1];
rz(-2.4534493) q[2];
sx q[2];
rz(-0.71758413) q[2];
sx q[2];
rz(2.6731499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8676973) q[1];
sx q[1];
rz(-0.58389837) q[1];
sx q[1];
rz(-0.03429596) q[1];
rz(-2.4998922) q[3];
sx q[3];
rz(-1.7880758) q[3];
sx q[3];
rz(-3.1171796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38146314) q[2];
sx q[2];
rz(-0.13801408) q[2];
sx q[2];
rz(-1.4786973) q[2];
rz(-2.0315157) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(0.64322513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4019302) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(0.31841835) q[0];
rz(3.1069801) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(-0.45077032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0863257) q[0];
sx q[0];
rz(-2.0811715) q[0];
sx q[0];
rz(-3.1253405) q[0];
rz(-pi) q[1];
rz(-2.320596) q[2];
sx q[2];
rz(-2.2911173) q[2];
sx q[2];
rz(-2.3483417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5391716) q[1];
sx q[1];
rz(-0.97214593) q[1];
sx q[1];
rz(2.6461698) q[1];
rz(-0.3698753) q[3];
sx q[3];
rz(-2.2910017) q[3];
sx q[3];
rz(0.50923075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.047711756) q[2];
sx q[2];
rz(-0.1897976) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(-0.098585248) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-0.70257598) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75673574) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(-3.1016438) q[0];
rz(2.3710251) q[1];
sx q[1];
rz(-2.4295085) q[1];
sx q[1];
rz(-0.22824731) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3632715) q[0];
sx q[0];
rz(-1.4822905) q[0];
sx q[0];
rz(-0.077544942) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1399916) q[2];
sx q[2];
rz(-1.2646156) q[2];
sx q[2];
rz(1.209335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.764594) q[1];
sx q[1];
rz(-1.5176763) q[1];
sx q[1];
rz(-0.082100987) q[1];
rz(0.41577783) q[3];
sx q[3];
rz(-0.97857394) q[3];
sx q[3];
rz(1.7433188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0780636) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(-2.883319) q[2];
rz(-2.9410948) q[3];
sx q[3];
rz(-2.343488) q[3];
sx q[3];
rz(2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-2.8343416) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(2.5603851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8756105) q[0];
sx q[0];
rz(-1.6020244) q[0];
sx q[0];
rz(1.1660378) q[0];
x q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4631066) q[1];
sx q[1];
rz(-2.3358279) q[1];
sx q[1];
rz(-1.8496978) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5804306) q[3];
sx q[3];
rz(-2.5303839) q[3];
sx q[3];
rz(-1.8451913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4622978) q[2];
sx q[2];
rz(-1.1335979) q[2];
sx q[2];
rz(3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.3616819) q[3];
sx q[3];
rz(-2.1112736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533326) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(-2.4326676) q[0];
rz(2.9290579) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(2.2841891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564726) q[0];
sx q[0];
rz(-1.5758638) q[0];
sx q[0];
rz(1.4410254) q[0];
rz(2.6231758) q[2];
sx q[2];
rz(-2.5441493) q[2];
sx q[2];
rz(3.0068676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58536968) q[1];
sx q[1];
rz(-1.959728) q[1];
sx q[1];
rz(2.5635368) q[1];
rz(-1.3579426) q[3];
sx q[3];
rz(-0.88485938) q[3];
sx q[3];
rz(2.6712772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4400441) q[2];
sx q[2];
rz(-2.7893119) q[2];
sx q[2];
rz(0.80553833) q[2];
rz(2.7696179) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(-2.7443547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032967903) q[0];
sx q[0];
rz(-1.911835) q[0];
sx q[0];
rz(-0.97880542) q[0];
rz(0.72273123) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-0.61789787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.582983) q[0];
sx q[0];
rz(-1.2961868) q[0];
sx q[0];
rz(-0.78297575) q[0];
rz(-pi) q[1];
rz(2.0892136) q[2];
sx q[2];
rz(-2.1627277) q[2];
sx q[2];
rz(0.059378864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4010716) q[1];
sx q[1];
rz(-2.8943672) q[1];
sx q[1];
rz(0.81584357) q[1];
rz(-pi) q[2];
rz(-2.3506451) q[3];
sx q[3];
rz(-3.0229048) q[3];
sx q[3];
rz(-0.044767901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.939398) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(0.00016577684) q[2];
rz(-0.58445066) q[3];
sx q[3];
rz(-1.0060468) q[3];
sx q[3];
rz(-2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38681876) q[0];
sx q[0];
rz(-1.5712354) q[0];
sx q[0];
rz(-1.5729217) q[0];
rz(1.3407002) q[1];
sx q[1];
rz(-2.0410213) q[1];
sx q[1];
rz(-1.6165728) q[1];
rz(0.45380886) q[2];
sx q[2];
rz(-2.1593675) q[2];
sx q[2];
rz(-0.16271954) q[2];
rz(-1.9289005) q[3];
sx q[3];
rz(-1.9869589) q[3];
sx q[3];
rz(-0.12369894) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
