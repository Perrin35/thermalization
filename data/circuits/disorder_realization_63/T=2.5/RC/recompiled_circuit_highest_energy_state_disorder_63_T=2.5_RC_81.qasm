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
rz(-2.1152273) q[0];
rz(-2.8431471) q[1];
sx q[1];
rz(-1.5084074) q[1];
sx q[1];
rz(-1.0025947) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63051414) q[0];
sx q[0];
rz(-1.7745695) q[0];
sx q[0];
rz(3.1183) q[0];
rz(-pi) q[1];
rz(-0.66291507) q[2];
sx q[2];
rz(-1.5272682) q[2];
sx q[2];
rz(-1.3452665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94384313) q[1];
sx q[1];
rz(-1.437206) q[1];
sx q[1];
rz(0.39033668) q[1];
rz(-0.63319541) q[3];
sx q[3];
rz(-1.060462) q[3];
sx q[3];
rz(-1.1309433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58874321) q[2];
sx q[2];
rz(-1.0007977) q[2];
sx q[2];
rz(0.0351077) q[2];
rz(-1.2037753) q[3];
sx q[3];
rz(-1.3222008) q[3];
sx q[3];
rz(-0.28425851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725919) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(-2.1373855) q[0];
rz(-0.12241157) q[1];
sx q[1];
rz(-1.4675843) q[1];
sx q[1];
rz(-2.4845128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9417038) q[0];
sx q[0];
rz(-2.7452069) q[0];
sx q[0];
rz(-0.21342166) q[0];
rz(1.2934712) q[2];
sx q[2];
rz(-1.4155626) q[2];
sx q[2];
rz(-0.6065795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66533684) q[1];
sx q[1];
rz(-1.0968535) q[1];
sx q[1];
rz(0.83935763) q[1];
x q[2];
rz(-1.3260001) q[3];
sx q[3];
rz(-0.80202937) q[3];
sx q[3];
rz(-1.0753461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5848026) q[2];
sx q[2];
rz(-2.5636797) q[2];
sx q[2];
rz(-0.59188262) q[2];
rz(-1.6651734) q[3];
sx q[3];
rz(-1.5341325) q[3];
sx q[3];
rz(-2.2679451) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0037435) q[0];
sx q[0];
rz(-0.20895222) q[0];
sx q[0];
rz(1.214622) q[0];
rz(-2.1851723) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(1.1716243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60739205) q[0];
sx q[0];
rz(-1.1063684) q[0];
sx q[0];
rz(-0.44363108) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20522988) q[2];
sx q[2];
rz(-2.2608888) q[2];
sx q[2];
rz(-3.1179259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4332648) q[1];
sx q[1];
rz(-1.9076132) q[1];
sx q[1];
rz(-2.5418806) q[1];
rz(-pi) q[2];
rz(1.5258777) q[3];
sx q[3];
rz(-1.4810331) q[3];
sx q[3];
rz(-2.1819851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9413635) q[2];
sx q[2];
rz(-1.0503294) q[2];
sx q[2];
rz(1.3776779) q[2];
rz(0.40669835) q[3];
sx q[3];
rz(-2.4923057) q[3];
sx q[3];
rz(2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7895301) q[0];
sx q[0];
rz(-1.2593513) q[0];
sx q[0];
rz(-1.9011185) q[0];
rz(2.4234096) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(-2.9026418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6278291) q[0];
sx q[0];
rz(-1.2830334) q[0];
sx q[0];
rz(-1.4445369) q[0];
x q[1];
rz(2.4282794) q[2];
sx q[2];
rz(-2.1916766) q[2];
sx q[2];
rz(2.5312405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.093555778) q[1];
sx q[1];
rz(-0.29680064) q[1];
sx q[1];
rz(2.0778627) q[1];
x q[2];
rz(1.4364868) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(-2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6255528) q[2];
sx q[2];
rz(-1.3935139) q[2];
sx q[2];
rz(1.7986521) q[2];
rz(2.2602153) q[3];
sx q[3];
rz(-2.9727029) q[3];
sx q[3];
rz(0.78099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.4137022) q[0];
sx q[0];
rz(-2.8975633) q[0];
sx q[0];
rz(1.6572886) q[0];
rz(-2.4877352) q[1];
sx q[1];
rz(-1.6203208) q[1];
sx q[1];
rz(-3.003655) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50820275) q[0];
sx q[0];
rz(-1.7195367) q[0];
sx q[0];
rz(1.2138019) q[0];
x q[1];
rz(-1.8191255) q[2];
sx q[2];
rz(-0.93298027) q[2];
sx q[2];
rz(-2.1030104) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.063542343) q[1];
sx q[1];
rz(-1.1415328) q[1];
sx q[1];
rz(3.0338366) q[1];
x q[2];
rz(-1.3624179) q[3];
sx q[3];
rz(-0.26431686) q[3];
sx q[3];
rz(-2.6770279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5518034) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(-3.0544082) q[2];
rz(-0.11463556) q[3];
sx q[3];
rz(-2.3264824) q[3];
sx q[3];
rz(2.7242928) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36599416) q[0];
sx q[0];
rz(-1.5394779) q[0];
sx q[0];
rz(-2.6566246) q[0];
rz(0.94351774) q[1];
sx q[1];
rz(-2.3285995) q[1];
sx q[1];
rz(1.2191204) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2526357) q[0];
sx q[0];
rz(-0.85272861) q[0];
sx q[0];
rz(0.12949334) q[0];
rz(-pi) q[1];
rz(1.6197422) q[2];
sx q[2];
rz(-2.8192047) q[2];
sx q[2];
rz(1.9677377) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4002126) q[1];
sx q[1];
rz(-1.6467491) q[1];
sx q[1];
rz(1.6981359) q[1];
rz(2.3818822) q[3];
sx q[3];
rz(-1.9518765) q[3];
sx q[3];
rz(-2.423513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9014827) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(0.94129747) q[2];
rz(2.8955722) q[3];
sx q[3];
rz(-1.6451719) q[3];
sx q[3];
rz(-1.3601607) q[3];
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
sx q[0];
rz(-pi/2) q[2];
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
rz(-1.7485626) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37993615) q[0];
sx q[0];
rz(-1.4967866) q[0];
sx q[0];
rz(-1.1940899) q[0];
rz(-pi) q[1];
rz(-1.2859401) q[2];
sx q[2];
rz(-2.1021038) q[2];
sx q[2];
rz(-1.667995) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6026615) q[1];
sx q[1];
rz(-1.4653413) q[1];
sx q[1];
rz(2.6698673) q[1];
rz(-pi) q[2];
rz(1.8461269) q[3];
sx q[3];
rz(-2.4440401) q[3];
sx q[3];
rz(-1.5740852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67036575) q[2];
sx q[2];
rz(-1.7134075) q[2];
sx q[2];
rz(-2.8426389) q[2];
rz(2.7361338) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(-0.56431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7120755) q[0];
sx q[0];
rz(-0.27292621) q[0];
sx q[0];
rz(-0.78936973) q[0];
rz(0.40496597) q[1];
sx q[1];
rz(-1.7908432) q[1];
sx q[1];
rz(0.49496034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32463418) q[0];
sx q[0];
rz(-1.0993891) q[0];
sx q[0];
rz(2.2105634) q[0];
x q[1];
rz(1.1210006) q[2];
sx q[2];
rz(-0.61857046) q[2];
sx q[2];
rz(-2.6976259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5304184) q[1];
sx q[1];
rz(-1.8137167) q[1];
sx q[1];
rz(-2.3526885) q[1];
rz(-pi) q[2];
rz(-2.6172753) q[3];
sx q[3];
rz(-0.90745196) q[3];
sx q[3];
rz(1.2330221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97883812) q[2];
sx q[2];
rz(-1.2252204) q[2];
sx q[2];
rz(0.2552574) q[2];
rz(3.1066331) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(-0.24622723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.404945) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(2.431562) q[0];
rz(1.0011477) q[1];
sx q[1];
rz(-0.89718693) q[1];
sx q[1];
rz(0.62090105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17711711) q[0];
sx q[0];
rz(-1.510448) q[0];
sx q[0];
rz(1.9326769) q[0];
rz(-pi) q[1];
rz(0.59130238) q[2];
sx q[2];
rz(-0.74858758) q[2];
sx q[2];
rz(2.8158514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2706333) q[1];
sx q[1];
rz(-2.7771725) q[1];
sx q[1];
rz(0.30521133) q[1];
rz(0.28212237) q[3];
sx q[3];
rz(-2.3446313) q[3];
sx q[3];
rz(-0.4235426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8093449) q[2];
sx q[2];
rz(-1.5231909) q[2];
sx q[2];
rz(2.8524032) q[2];
rz(-2.0293763) q[3];
sx q[3];
rz(-0.29505348) q[3];
sx q[3];
rz(2.1212063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6638829) q[0];
sx q[0];
rz(-1.3417256) q[0];
sx q[0];
rz(-0.92673242) q[0];
rz(-2.0345188) q[1];
sx q[1];
rz(-1.6278382) q[1];
sx q[1];
rz(2.5735997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999789) q[0];
sx q[0];
rz(-1.2547412) q[0];
sx q[0];
rz(0.63412447) q[0];
rz(-pi) q[1];
rz(-0.51189152) q[2];
sx q[2];
rz(-1.8020168) q[2];
sx q[2];
rz(1.1752886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2462594) q[1];
sx q[1];
rz(-2.3998108) q[1];
sx q[1];
rz(-2.1269024) q[1];
rz(2.5015478) q[3];
sx q[3];
rz(-1.1403822) q[3];
sx q[3];
rz(-1.1120344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.052281436) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(1.4386162) q[2];
rz(0.29664052) q[3];
sx q[3];
rz(-2.7999925) q[3];
sx q[3];
rz(-2.6945485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.457837) q[0];
sx q[0];
rz(-2.044027) q[0];
sx q[0];
rz(-2.5720163) q[0];
rz(0.33869047) q[1];
sx q[1];
rz(-1.5298005) q[1];
sx q[1];
rz(-1.6385967) q[1];
rz(-0.25119943) q[2];
sx q[2];
rz(-1.3710556) q[2];
sx q[2];
rz(-2.5547169) q[2];
rz(-2.114645) q[3];
sx q[3];
rz(-0.31796496) q[3];
sx q[3];
rz(1.8798469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
