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
rz(1.0032049) q[0];
sx q[0];
rz(-1.2004852) q[0];
sx q[0];
rz(1.0263654) q[0];
rz(0.29844555) q[1];
sx q[1];
rz(-1.6331853) q[1];
sx q[1];
rz(1.0025947) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93556773) q[0];
sx q[0];
rz(-1.5479857) q[0];
sx q[0];
rz(-1.3669694) q[0];
x q[1];
rz(-3.0709362) q[2];
sx q[2];
rz(-0.66412726) q[2];
sx q[2];
rz(2.9717367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1977495) q[1];
sx q[1];
rz(-1.437206) q[1];
sx q[1];
rz(2.751256) q[1];
rz(-pi) q[2];
rz(-0.63319541) q[3];
sx q[3];
rz(-1.060462) q[3];
sx q[3];
rz(2.0106493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5528494) q[2];
sx q[2];
rz(-1.0007977) q[2];
sx q[2];
rz(-0.0351077) q[2];
rz(1.9378174) q[3];
sx q[3];
rz(-1.3222008) q[3];
sx q[3];
rz(-0.28425851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069000706) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(-2.1373855) q[0];
rz(3.0191811) q[1];
sx q[1];
rz(-1.6740084) q[1];
sx q[1];
rz(2.4845128) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5733684) q[0];
sx q[0];
rz(-1.4889297) q[0];
sx q[0];
rz(-0.38827814) q[0];
x q[1];
rz(-1.051527) q[2];
sx q[2];
rz(-0.31683559) q[2];
sx q[2];
rz(-0.46665442) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.66533684) q[1];
sx q[1];
rz(-2.0447392) q[1];
sx q[1];
rz(-2.302235) q[1];
x q[2];
rz(1.8155926) q[3];
sx q[3];
rz(-2.3395633) q[3];
sx q[3];
rz(1.0753461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5848026) q[2];
sx q[2];
rz(-2.5636797) q[2];
sx q[2];
rz(-0.59188262) q[2];
rz(1.6651734) q[3];
sx q[3];
rz(-1.6074601) q[3];
sx q[3];
rz(-2.2679451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1378491) q[0];
sx q[0];
rz(-0.20895222) q[0];
sx q[0];
rz(-1.9269706) q[0];
rz(0.95642033) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(1.1716243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5342006) q[0];
sx q[0];
rz(-1.1063684) q[0];
sx q[0];
rz(-0.44363108) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2713216) q[2];
sx q[2];
rz(-1.7286125) q[2];
sx q[2];
rz(-1.6788788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4332648) q[1];
sx q[1];
rz(-1.2339795) q[1];
sx q[1];
rz(-0.59971209) q[1];
x q[2];
rz(1.5258777) q[3];
sx q[3];
rz(-1.4810331) q[3];
sx q[3];
rz(-2.1819851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9413635) q[2];
sx q[2];
rz(-1.0503294) q[2];
sx q[2];
rz(1.7639147) q[2];
rz(0.40669835) q[3];
sx q[3];
rz(-0.64928693) q[3];
sx q[3];
rz(-2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3520626) q[0];
sx q[0];
rz(-1.2593513) q[0];
sx q[0];
rz(1.9011185) q[0];
rz(-2.4234096) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(2.9026418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5137635) q[0];
sx q[0];
rz(-1.2830334) q[0];
sx q[0];
rz(-1.4445369) q[0];
rz(-pi) q[1];
rz(-2.3117606) q[2];
sx q[2];
rz(-0.90832635) q[2];
sx q[2];
rz(-1.5892513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0480369) q[1];
sx q[1];
rz(-0.29680064) q[1];
sx q[1];
rz(2.0778627) q[1];
rz(-pi) q[2];
rz(-1.7051058) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(-2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6255528) q[2];
sx q[2];
rz(-1.3935139) q[2];
sx q[2];
rz(1.3429406) q[2];
rz(2.2602153) q[3];
sx q[3];
rz(-0.16888976) q[3];
sx q[3];
rz(-0.78099293) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4137022) q[0];
sx q[0];
rz(-2.8975633) q[0];
sx q[0];
rz(1.6572886) q[0];
rz(-0.65385747) q[1];
sx q[1];
rz(-1.5212719) q[1];
sx q[1];
rz(0.13793764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0237831) q[0];
sx q[0];
rz(-1.217917) q[0];
sx q[0];
rz(2.9830052) q[0];
rz(-pi) q[1];
rz(-0.32021145) q[2];
sx q[2];
rz(-0.67812014) q[2];
sx q[2];
rz(1.7004287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5893464) q[1];
sx q[1];
rz(-1.4728496) q[1];
sx q[1];
rz(2.0022654) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7791747) q[3];
sx q[3];
rz(-0.26431686) q[3];
sx q[3];
rz(2.6770279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5518034) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(-3.0544082) q[2];
rz(-3.0269571) q[3];
sx q[3];
rz(-0.81511027) q[3];
sx q[3];
rz(-0.4172999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7755985) q[0];
sx q[0];
rz(-1.6021148) q[0];
sx q[0];
rz(-0.48496801) q[0];
rz(2.1980749) q[1];
sx q[1];
rz(-0.8129932) q[1];
sx q[1];
rz(-1.9224723) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6935558) q[0];
sx q[0];
rz(-0.72759923) q[0];
sx q[0];
rz(-1.4240525) q[0];
rz(-1.5218505) q[2];
sx q[2];
rz(-0.32238797) q[2];
sx q[2];
rz(1.173855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18029848) q[1];
sx q[1];
rz(-1.4438259) q[1];
sx q[1];
rz(3.0650223) q[1];
rz(-2.6146919) q[3];
sx q[3];
rz(-0.8325067) q[3];
sx q[3];
rz(-0.47952521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9014827) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(-0.94129747) q[2];
rz(2.8955722) q[3];
sx q[3];
rz(-1.6451719) q[3];
sx q[3];
rz(-1.3601607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0435903) q[0];
sx q[0];
rz(-2.0397546) q[0];
sx q[0];
rz(-1.164042) q[0];
rz(-1.7835435) q[1];
sx q[1];
rz(-0.86654228) q[1];
sx q[1];
rz(1.39303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2201033) q[0];
sx q[0];
rz(-1.1951726) q[0];
sx q[0];
rz(-0.079567841) q[0];
rz(2.5921408) q[2];
sx q[2];
rz(-1.3260815) q[2];
sx q[2];
rz(-0.050087226) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.021780304) q[1];
sx q[1];
rz(-2.0396898) q[1];
sx q[1];
rz(1.6890668) q[1];
rz(-pi) q[2];
rz(-1.8461269) q[3];
sx q[3];
rz(-2.4440401) q[3];
sx q[3];
rz(-1.5675074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4712269) q[2];
sx q[2];
rz(-1.4281851) q[2];
sx q[2];
rz(-0.29895374) q[2];
rz(0.40545884) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(0.56431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4295171) q[0];
sx q[0];
rz(-0.27292621) q[0];
sx q[0];
rz(-2.3522229) q[0];
rz(2.7366267) q[1];
sx q[1];
rz(-1.7908432) q[1];
sx q[1];
rz(2.6466323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69823658) q[0];
sx q[0];
rz(-2.3670417) q[0];
sx q[0];
rz(-2.2775035) q[0];
rz(-pi) q[1];
rz(-1.1210006) q[2];
sx q[2];
rz(-2.5230222) q[2];
sx q[2];
rz(0.44396675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72515762) q[1];
sx q[1];
rz(-0.817653) q[1];
sx q[1];
rz(2.8055951) q[1];
rz(2.1405419) q[3];
sx q[3];
rz(-2.3214139) q[3];
sx q[3];
rz(-2.6626183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1627545) q[2];
sx q[2];
rz(-1.9163722) q[2];
sx q[2];
rz(-0.2552574) q[2];
rz(-0.034959547) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(2.8953654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.404945) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(-2.431562) q[0];
rz(-2.1404449) q[1];
sx q[1];
rz(-0.89718693) q[1];
sx q[1];
rz(-2.5206916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17711711) q[0];
sx q[0];
rz(-1.510448) q[0];
sx q[0];
rz(-1.9326769) q[0];
rz(-pi) q[1];
rz(2.5502903) q[2];
sx q[2];
rz(-0.74858758) q[2];
sx q[2];
rz(0.32574124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1961851) q[1];
sx q[1];
rz(-1.2239478) q[1];
sx q[1];
rz(1.4566684) q[1];
x q[2];
rz(2.3647978) q[3];
sx q[3];
rz(-1.771254) q[3];
sx q[3];
rz(-1.3471229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3322477) q[2];
sx q[2];
rz(-1.5231909) q[2];
sx q[2];
rz(0.28918949) q[2];
rz(2.0293763) q[3];
sx q[3];
rz(-0.29505348) q[3];
sx q[3];
rz(-2.1212063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4777098) q[0];
sx q[0];
rz(-1.7998671) q[0];
sx q[0];
rz(-2.2148602) q[0];
rz(-2.0345188) q[1];
sx q[1];
rz(-1.6278382) q[1];
sx q[1];
rz(-0.56799299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3044514) q[0];
sx q[0];
rz(-0.97263179) q[0];
sx q[0];
rz(1.9564081) q[0];
rz(-pi) q[1];
rz(0.51189152) q[2];
sx q[2];
rz(-1.3395759) q[2];
sx q[2];
rz(1.1752886) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8958139) q[1];
sx q[1];
rz(-1.2061283) q[1];
sx q[1];
rz(-0.90943894) q[1];
rz(-pi) q[2];
rz(-2.5015478) q[3];
sx q[3];
rz(-2.0012104) q[3];
sx q[3];
rz(2.0295582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.052281436) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(1.4386162) q[2];
rz(-2.8449521) q[3];
sx q[3];
rz(-0.34160015) q[3];
sx q[3];
rz(2.6945485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68375568) q[0];
sx q[0];
rz(-1.0975657) q[0];
sx q[0];
rz(0.56957635) q[0];
rz(2.8029022) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(2.4581134) q[2];
sx q[2];
rz(-0.31961315) q[2];
sx q[2];
rz(1.4994958) q[2];
rz(0.16868563) q[3];
sx q[3];
rz(-1.2999693) q[3];
sx q[3];
rz(-0.6948605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
