OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(-0.93710605) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694323) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(-0.62299563) q[0];
rz(-2.063077) q[2];
sx q[2];
rz(-2.1247851) q[2];
sx q[2];
rz(-2.2433777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1569251) q[1];
sx q[1];
rz(-0.19913864) q[1];
sx q[1];
rz(1.9763293) q[1];
rz(-pi) q[2];
rz(0.097221656) q[3];
sx q[3];
rz(-1.3029649) q[3];
sx q[3];
rz(-0.83084805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15930882) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(0.75749767) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(-2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57698292) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(0.27145162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46389929) q[0];
sx q[0];
rz(-1.1743744) q[0];
sx q[0];
rz(1.0248313) q[0];
rz(-pi) q[1];
rz(-1.8797305) q[2];
sx q[2];
rz(-1.075282) q[2];
sx q[2];
rz(0.77906424) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6352306) q[1];
sx q[1];
rz(-0.87345424) q[1];
sx q[1];
rz(-0.49763775) q[1];
rz(0.35630393) q[3];
sx q[3];
rz(-0.84173991) q[3];
sx q[3];
rz(-0.19405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(-1.4097759) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975824) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(-0.14257167) q[0];
rz(-1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81461834) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(2.4787865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9833343) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(-1.0227026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28029728) q[1];
sx q[1];
rz(-1.1251083) q[1];
sx q[1];
rz(0.46511005) q[1];
x q[2];
rz(-1.7861869) q[3];
sx q[3];
rz(-1.2743909) q[3];
sx q[3];
rz(2.4733558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(-0.16734853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(-2.1787815) q[0];
rz(2.6722233) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(3.1406291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2036635) q[0];
sx q[0];
rz(-2.2712049) q[0];
sx q[0];
rz(-0.80313375) q[0];
x q[1];
rz(3.1374627) q[2];
sx q[2];
rz(-3.0175856) q[2];
sx q[2];
rz(2.4069402) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0469989) q[1];
sx q[1];
rz(-1.440289) q[1];
sx q[1];
rz(-0.99884896) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56941454) q[3];
sx q[3];
rz(-1.4953519) q[3];
sx q[3];
rz(-0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(-0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951293) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(2.8809663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31291744) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(-0.28767985) q[0];
x q[1];
rz(-0.79029681) q[2];
sx q[2];
rz(-2.7919263) q[2];
sx q[2];
rz(2.5555573) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6071636) q[1];
sx q[1];
rz(-0.57463127) q[1];
sx q[1];
rz(-0.032392153) q[1];
x q[2];
rz(1.4277677) q[3];
sx q[3];
rz(-1.4266532) q[3];
sx q[3];
rz(3.0343057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2300718) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(0.22932209) q[2];
rz(0.54245943) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77263537) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.649958) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6468069) q[0];
sx q[0];
rz(-1.8456869) q[0];
sx q[0];
rz(-1.5216212) q[0];
rz(-pi) q[1];
rz(2.6238407) q[2];
sx q[2];
rz(-1.3377829) q[2];
sx q[2];
rz(1.053996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9767178) q[1];
sx q[1];
rz(-2.4396982) q[1];
sx q[1];
rz(-0.81275069) q[1];
rz(-1.2392483) q[3];
sx q[3];
rz(-0.22833951) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-0.12167715) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50169045) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-0.34696607) q[0];
rz(-pi) q[1];
rz(-0.17860883) q[2];
sx q[2];
rz(-0.046769301) q[2];
sx q[2];
rz(2.5472819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.788159) q[1];
sx q[1];
rz(-2.806059) q[1];
sx q[1];
rz(-2.5787756) q[1];
rz(-pi) q[2];
rz(-1.2225371) q[3];
sx q[3];
rz(-1.5969443) q[3];
sx q[3];
rz(0.98046434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.014331269) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(2.2793615) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7858793) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(2.1070811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1326133) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10266281) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(1.3753124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58971436) q[1];
sx q[1];
rz(-1.0052181) q[1];
sx q[1];
rz(-3.0977071) q[1];
rz(-1.4626059) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(-3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2404279) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-2.365716) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(-0.016816703) q[0];
rz(0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1982225) q[0];
sx q[0];
rz(-1.3226042) q[0];
sx q[0];
rz(2.7944837) q[0];
rz(-2.0428033) q[2];
sx q[2];
rz(-1.325377) q[2];
sx q[2];
rz(2.5185042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(-2.1437953) q[1];
x q[2];
rz(1.2905144) q[3];
sx q[3];
rz(-2.9119769) q[3];
sx q[3];
rz(1.5667825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4497711) q[2];
sx q[2];
rz(-1.3042973) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(-1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(-2.418628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3176214) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(-2.1872107) q[0];
rz(-pi) q[1];
rz(1.9777708) q[2];
sx q[2];
rz(-2.7343035) q[2];
sx q[2];
rz(-2.2182857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5889659) q[1];
sx q[1];
rz(-0.50682658) q[1];
sx q[1];
rz(-1.6848906) q[1];
x q[2];
rz(-2.3530657) q[3];
sx q[3];
rz(-1.8619814) q[3];
sx q[3];
rz(-2.9332719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(2.6160713) q[2];
rz(-0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(-2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(1.0992959) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(-0.47808403) q[2];
sx q[2];
rz(-0.95778428) q[2];
sx q[2];
rz(-0.20245353) q[2];
rz(-1.7952193) q[3];
sx q[3];
rz(-1.8891469) q[3];
sx q[3];
rz(-0.3286152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
