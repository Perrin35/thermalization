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
rz(0.17231365) q[0];
sx q[0];
rz(-2.9336689) q[0];
sx q[0];
rz(1.2907668) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(-0.42833498) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243619) q[0];
sx q[0];
rz(-2.9901282) q[0];
sx q[0];
rz(-2.3343655) q[0];
x q[1];
rz(1.8646127) q[2];
sx q[2];
rz(-0.91962469) q[2];
sx q[2];
rz(-3.0319954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60778695) q[1];
sx q[1];
rz(-1.6139107) q[1];
sx q[1];
rz(-0.73854609) q[1];
rz(-pi) q[2];
rz(-0.21875225) q[3];
sx q[3];
rz(-1.7667337) q[3];
sx q[3];
rz(0.71813717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20994818) q[2];
sx q[2];
rz(-2.1817709) q[2];
sx q[2];
rz(1.9317365) q[2];
rz(-3.0620388) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(2.615926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88456589) q[0];
sx q[0];
rz(-2.631812) q[0];
sx q[0];
rz(0.38401815) q[0];
rz(-0.89029038) q[1];
sx q[1];
rz(-1.236981) q[1];
sx q[1];
rz(2.8382137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6174406) q[0];
sx q[0];
rz(-2.5684803) q[0];
sx q[0];
rz(-1.9277431) q[0];
rz(-pi) q[1];
rz(1.134642) q[2];
sx q[2];
rz(-0.56098962) q[2];
sx q[2];
rz(2.2962388) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0133923) q[1];
sx q[1];
rz(-1.5085932) q[1];
sx q[1];
rz(1.661646) q[1];
x q[2];
rz(1.1214662) q[3];
sx q[3];
rz(-2.1746965) q[3];
sx q[3];
rz(0.7627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7552135) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(-2.4066822) q[2];
rz(-2.2928061) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-2.7809704) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81247771) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(-0.078711674) q[0];
rz(2.3178237) q[1];
sx q[1];
rz(-2.0562101) q[1];
sx q[1];
rz(-1.8471921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9322008) q[0];
sx q[0];
rz(-0.086361445) q[0];
sx q[0];
rz(-2.3461653) q[0];
rz(-pi) q[1];
rz(0.95125385) q[2];
sx q[2];
rz(-1.9980556) q[2];
sx q[2];
rz(-0.45742971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0261198) q[1];
sx q[1];
rz(-1.8515665) q[1];
sx q[1];
rz(1.6927648) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2749316) q[3];
sx q[3];
rz(-1.7539644) q[3];
sx q[3];
rz(-0.66732348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2417629) q[2];
sx q[2];
rz(-2.7624942) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(-1.624931) q[3];
sx q[3];
rz(-0.95976019) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4888332) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(0.56831992) q[0];
rz(1.142451) q[1];
sx q[1];
rz(-1.3520974) q[1];
sx q[1];
rz(-1.4521339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2556011) q[0];
sx q[0];
rz(-1.7924249) q[0];
sx q[0];
rz(2.5021573) q[0];
rz(-pi) q[1];
rz(-2.9617519) q[2];
sx q[2];
rz(-1.4786063) q[2];
sx q[2];
rz(-1.8123019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6664847) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(0.38352769) q[1];
rz(-1.9098572) q[3];
sx q[3];
rz(-1.4022458) q[3];
sx q[3];
rz(1.7730561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(0.076889195) q[2];
rz(2.7422089) q[3];
sx q[3];
rz(-2.2279584) q[3];
sx q[3];
rz(0.54401773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15453108) q[0];
sx q[0];
rz(-0.35683826) q[0];
sx q[0];
rz(2.8416908) q[0];
rz(1.9918282) q[1];
sx q[1];
rz(-2.1297784) q[1];
sx q[1];
rz(-0.48847517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66118427) q[0];
sx q[0];
rz(-1.6035491) q[0];
sx q[0];
rz(1.3854273) q[0];
rz(1.064642) q[2];
sx q[2];
rz(-2.1035668) q[2];
sx q[2];
rz(-1.8440994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2327959) q[1];
sx q[1];
rz(-0.98728647) q[1];
sx q[1];
rz(-1.5481434) q[1];
x q[2];
rz(2.7882448) q[3];
sx q[3];
rz(-0.67253695) q[3];
sx q[3];
rz(1.3143244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7601295) q[2];
sx q[2];
rz(-0.13801408) q[2];
sx q[2];
rz(1.6628954) q[2];
rz(-2.0315157) q[3];
sx q[3];
rz(-2.1434982) q[3];
sx q[3];
rz(-2.4983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.7396624) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(-2.8231743) q[0];
rz(-3.1069801) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(0.45077032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0863257) q[0];
sx q[0];
rz(-2.0811715) q[0];
sx q[0];
rz(3.1253405) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66025205) q[2];
sx q[2];
rz(-0.98838943) q[2];
sx q[2];
rz(-2.9803515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3035982) q[1];
sx q[1];
rz(-0.7571836) q[1];
sx q[1];
rz(-2.1794469) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7717174) q[3];
sx q[3];
rz(-0.85059097) q[3];
sx q[3];
rz(2.6323619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.047711756) q[2];
sx q[2];
rz(-2.9517951) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(0.098585248) q[3];
sx q[3];
rz(-0.74461377) q[3];
sx q[3];
rz(-2.4390167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.75673574) q[0];
sx q[0];
rz(-1.1133794) q[0];
sx q[0];
rz(-3.1016438) q[0];
rz(2.3710251) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(0.22824731) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21439274) q[0];
sx q[0];
rz(-1.4935555) q[0];
sx q[0];
rz(1.4820251) q[0];
rz(-pi) q[1];
rz(1.5657317) q[2];
sx q[2];
rz(-0.3061848) q[2];
sx q[2];
rz(-1.2146467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3769987) q[1];
sx q[1];
rz(-1.6239163) q[1];
sx q[1];
rz(3.0594917) q[1];
x q[2];
rz(-1.030121) q[3];
sx q[3];
rz(-0.70899963) q[3];
sx q[3];
rz(-1.0741155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0780636) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(2.883319) q[2];
rz(-0.20049788) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.9789326) q[0];
sx q[0];
rz(-1.7698092) q[0];
sx q[0];
rz(-2.8343416) q[0];
rz(-1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(-0.58120751) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2659822) q[0];
sx q[0];
rz(-1.5395682) q[0];
sx q[0];
rz(-1.1660378) q[0];
rz(2.4318784) q[2];
sx q[2];
rz(-1.4462556) q[2];
sx q[2];
rz(2.1309851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30347428) q[1];
sx q[1];
rz(-1.7707157) q[1];
sx q[1];
rz(-0.78472991) q[1];
x q[2];
rz(0.53544541) q[3];
sx q[3];
rz(-1.8811444) q[3];
sx q[3];
rz(0.20099881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4622978) q[2];
sx q[2];
rz(-2.0079948) q[2];
sx q[2];
rz(-3.1104258) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(-1.0303191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088260055) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(0.70892507) q[0];
rz(-0.21253474) q[1];
sx q[1];
rz(-1.0147076) q[1];
sx q[1];
rz(-0.85740352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5564726) q[0];
sx q[0];
rz(-1.5657288) q[0];
sx q[0];
rz(1.7005672) q[0];
rz(-pi) q[1];
rz(-0.51841684) q[2];
sx q[2];
rz(-2.5441493) q[2];
sx q[2];
rz(3.0068676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58536968) q[1];
sx q[1];
rz(-1.1818647) q[1];
sx q[1];
rz(0.57805581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69714947) q[3];
sx q[3];
rz(-1.4065885) q[3];
sx q[3];
rz(-2.177161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032967903) q[0];
sx q[0];
rz(-1.911835) q[0];
sx q[0];
rz(-2.1627872) q[0];
rz(2.4188614) q[1];
sx q[1];
rz(-2.0131854) q[1];
sx q[1];
rz(2.5236948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5586097) q[0];
sx q[0];
rz(-1.2961868) q[0];
sx q[0];
rz(-2.3586169) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63510908) q[2];
sx q[2];
rz(-0.76580566) q[2];
sx q[2];
rz(-0.85585153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1714461) q[1];
sx q[1];
rz(-1.3916124) q[1];
sx q[1];
rz(-0.1712562) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4862107) q[3];
sx q[3];
rz(-1.4874377) q[3];
sx q[3];
rz(-0.74970923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2021947) q[2];
sx q[2];
rz(-1.7859744) q[2];
sx q[2];
rz(-3.1414269) q[2];
rz(0.58445066) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(-2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38681876) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(-1.3407002) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(-0.93201119) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(1.9289005) q[3];
sx q[3];
rz(-1.1546338) q[3];
sx q[3];
rz(3.0178937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
