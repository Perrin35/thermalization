OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(0.23101097) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(-0.45431554) q[1];
sx q[1];
rz(1.2872202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481349) q[0];
sx q[0];
rz(-1.0890111) q[0];
sx q[0];
rz(-0.30695494) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24536774) q[2];
sx q[2];
rz(-1.7984219) q[2];
sx q[2];
rz(2.1466704) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8227777) q[1];
sx q[1];
rz(-2.9479369) q[1];
sx q[1];
rz(0.53656399) q[1];
rz(1.9735967) q[3];
sx q[3];
rz(-2.9514545) q[3];
sx q[3];
rz(2.8791752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7370558) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(1.1245842) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(-1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(-0.65688175) q[0];
rz(2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(2.2064256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864185) q[0];
sx q[0];
rz(-1.6011229) q[0];
sx q[0];
rz(-3.0369989) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6627533) q[2];
sx q[2];
rz(-1.8269681) q[2];
sx q[2];
rz(2.0785477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3352901) q[1];
sx q[1];
rz(-2.9893901) q[1];
sx q[1];
rz(-2.8782513) q[1];
x q[2];
rz(-2.786119) q[3];
sx q[3];
rz(-1.7661816) q[3];
sx q[3];
rz(3.1110142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8344581) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(-1.2163986) q[2];
rz(-0.52792102) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(-1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6127748) q[0];
sx q[0];
rz(-0.82614326) q[0];
sx q[0];
rz(2.5566027) q[0];
rz(0.12475573) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(-1.4488719) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18066474) q[0];
sx q[0];
rz(-1.573473) q[0];
sx q[0];
rz(1.5737783) q[0];
rz(-pi) q[1];
rz(-0.18529202) q[2];
sx q[2];
rz(-3.0681562) q[2];
sx q[2];
rz(2.9543608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8068778) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(-0.54912864) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33749859) q[3];
sx q[3];
rz(-1.322116) q[3];
sx q[3];
rz(1.4583704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8831732) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(1.4460571) q[2];
rz(2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739173) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(-2.0539334) q[0];
rz(2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-2.1813724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86726928) q[0];
sx q[0];
rz(-1.835863) q[0];
sx q[0];
rz(-3.1004219) q[0];
rz(0.029898568) q[2];
sx q[2];
rz(-1.471386) q[2];
sx q[2];
rz(1.7123897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5615446) q[1];
sx q[1];
rz(-2.5001038) q[1];
sx q[1];
rz(-2.3381691) q[1];
rz(-pi) q[2];
rz(0.26953617) q[3];
sx q[3];
rz(-1.663066) q[3];
sx q[3];
rz(2.5470981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9356161) q[2];
sx q[2];
rz(-1.874186) q[2];
sx q[2];
rz(-1.4975632) q[2];
rz(2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(-2.5373051) q[0];
rz(1.9970278) q[1];
sx q[1];
rz(-1.4860169) q[1];
sx q[1];
rz(0.42246517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4504964) q[0];
sx q[0];
rz(-1.7393149) q[0];
sx q[0];
rz(3.0590579) q[0];
rz(-0.79549148) q[2];
sx q[2];
rz(-0.73372148) q[2];
sx q[2];
rz(-2.7249641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.656184) q[1];
sx q[1];
rz(-2.2586169) q[1];
sx q[1];
rz(1.4774035) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9129539) q[3];
sx q[3];
rz(-1.2099724) q[3];
sx q[3];
rz(0.35508852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(-2.4853415) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(0.58073616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-2.1298213) q[0];
sx q[0];
rz(0.62498012) q[0];
rz(-0.75195733) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(1.5350852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4706375) q[0];
sx q[0];
rz(-0.29330253) q[0];
sx q[0];
rz(-2.9706012) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67120303) q[2];
sx q[2];
rz(-0.83121383) q[2];
sx q[2];
rz(-2.9321456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9823204) q[1];
sx q[1];
rz(-0.96010027) q[1];
sx q[1];
rz(-1.8507694) q[1];
rz(-2.1702437) q[3];
sx q[3];
rz(-2.0347715) q[3];
sx q[3];
rz(0.52677192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(-0.9291741) q[2];
rz(0.82516986) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55649844) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.3744542) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(0.62320954) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334367) q[0];
sx q[0];
rz(-1.2386453) q[0];
sx q[0];
rz(2.9812814) q[0];
rz(2.8002732) q[2];
sx q[2];
rz(-0.53487294) q[2];
sx q[2];
rz(-0.90082263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0125) q[1];
sx q[1];
rz(-2.9575037) q[1];
sx q[1];
rz(-2.2232008) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74759746) q[3];
sx q[3];
rz(-1.7204525) q[3];
sx q[3];
rz(0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(3.1189611) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(2.8959494) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(-1.0173215) q[0];
rz(1.3941049) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(0.62754935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857972) q[0];
sx q[0];
rz(-1.1748519) q[0];
sx q[0];
rz(0.93314472) q[0];
rz(1.4862138) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(-1.6135474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1422954) q[1];
sx q[1];
rz(-0.31041708) q[1];
sx q[1];
rz(2.4503972) q[1];
rz(-pi) q[2];
rz(1.2305607) q[3];
sx q[3];
rz(-2.0215394) q[3];
sx q[3];
rz(2.4204202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5158186) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(3.0913894) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(-0.075411782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5589767) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(-2.8726752) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(-0.1300098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63142473) q[0];
sx q[0];
rz(-1.893259) q[0];
sx q[0];
rz(0.23309751) q[0];
rz(-1.9275877) q[2];
sx q[2];
rz(-1.2399106) q[2];
sx q[2];
rz(-0.39168229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42327729) q[1];
sx q[1];
rz(-2.1287103) q[1];
sx q[1];
rz(-1.7608791) q[1];
rz(-0.49874108) q[3];
sx q[3];
rz(-2.1864656) q[3];
sx q[3];
rz(1.2241253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93398634) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-0.096573528) q[2];
rz(-1.5038331) q[3];
sx q[3];
rz(-1.3373172) q[3];
sx q[3];
rz(2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0132975) q[0];
sx q[0];
rz(-2.9533563) q[0];
sx q[0];
rz(0.53303322) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(-1.689555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5851327) q[0];
sx q[0];
rz(-2.1341354) q[0];
sx q[0];
rz(0.4507395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31923652) q[2];
sx q[2];
rz(-1.4672654) q[2];
sx q[2];
rz(1.5543303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9033135) q[1];
sx q[1];
rz(-0.36064816) q[1];
sx q[1];
rz(-3.0566932) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.96569) q[3];
sx q[3];
rz(-2.6697192) q[3];
sx q[3];
rz(2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4286246) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(3.1089605) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.2266351) q[3];
sx q[3];
rz(-2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3451155) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-2.6061229) q[2];
sx q[2];
rz(-1.3146123) q[2];
sx q[2];
rz(-1.8563369) q[2];
rz(-3.0211021) q[3];
sx q[3];
rz(-2.3269666) q[3];
sx q[3];
rz(2.5757488) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
