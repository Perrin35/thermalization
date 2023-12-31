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
rz(2.2044866) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694323) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(-0.62299563) q[0];
rz(-pi) q[1];
rz(-2.063077) q[2];
sx q[2];
rz(-2.1247851) q[2];
sx q[2];
rz(0.89821494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.01552445) q[1];
sx q[1];
rz(-1.6489195) q[1];
sx q[1];
rz(-1.7541581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9107781) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(-1.9576548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-2.0479726) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(-0.27145162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6730246) q[0];
sx q[0];
rz(-2.4789171) q[0];
sx q[0];
rz(0.89232348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51241264) q[2];
sx q[2];
rz(-0.57704848) q[2];
sx q[2];
rz(1.3702099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2084864) q[1];
sx q[1];
rz(-0.83175627) q[1];
sx q[1];
rz(2.0887124) q[1];
x q[2];
rz(0.35630393) q[3];
sx q[3];
rz(-2.2998527) q[3];
sx q[3];
rz(0.19405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(-0.57717741) q[2];
rz(0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(-1.3525195) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81461834) q[0];
sx q[0];
rz(-1.0083535) q[0];
sx q[0];
rz(-2.4787865) q[0];
rz(-pi) q[1];
rz(-2.4019037) q[2];
sx q[2];
rz(-1.8834754) q[2];
sx q[2];
rz(-0.82175053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-2.0164844) q[1];
sx q[1];
rz(0.46511005) q[1];
x q[2];
rz(1.3554058) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(0.6682369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(2.9742441) q[3];
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
rz(1.7053112) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(-0.96281111) q[0];
rz(-2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(3.1406291) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9199333) q[0];
sx q[0];
rz(-0.98826212) q[0];
sx q[0];
rz(-0.68908738) q[0];
x q[1];
rz(-0.12400603) q[2];
sx q[2];
rz(-1.5702855) q[2];
sx q[2];
rz(-0.83204568) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8648659) q[1];
sx q[1];
rz(-2.556567) q[1];
sx q[1];
rz(1.8086955) q[1];
x q[2];
rz(0.56941454) q[3];
sx q[3];
rz(-1.6462407) q[3];
sx q[3];
rz(0.30573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0756388) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-0.11165079) q[2];
rz(-2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(2.7850889) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(2.8809663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31291744) q[0];
sx q[0];
rz(-0.14794359) q[0];
sx q[0];
rz(0.28767985) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8243276) q[2];
sx q[2];
rz(-1.3273444) q[2];
sx q[2];
rz(-1.4075116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6071636) q[1];
sx q[1];
rz(-0.57463127) q[1];
sx q[1];
rz(3.1092005) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[1];
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
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(-1.4916346) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0522239) q[0];
sx q[0];
rz(-1.5234689) q[0];
sx q[0];
rz(-2.8663859) q[0];
rz(-1.3041777) q[2];
sx q[2];
rz(-1.068371) q[2];
sx q[2];
rz(2.7555639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1097088) q[1];
sx q[1];
rz(-1.110853) q[1];
sx q[1];
rz(1.0201395) q[1];
x q[2];
rz(1.7870951) q[3];
sx q[3];
rz(-1.6445451) q[3];
sx q[3];
rz(2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(0.027739851) q[2];
rz(-0.49267832) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-0.12167715) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(2.7391403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284751) q[0];
sx q[0];
rz(-0.74281456) q[0];
sx q[0];
rz(1.166056) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0955663) q[2];
sx q[2];
rz(-1.5791025) q[2];
sx q[2];
rz(1.9866895) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.23558815) q[1];
sx q[1];
rz(-1.2885805) q[1];
sx q[1];
rz(1.7547592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6472858) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.014331269) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-1.0345116) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5401934) q[0];
sx q[0];
rz(-1.7978151) q[0];
sx q[0];
rz(-0.17417234) q[0];
x q[1];
rz(-2.5402252) q[2];
sx q[2];
rz(-3.0172918) q[2];
sx q[2];
rz(-0.79324317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58971436) q[1];
sx q[1];
rz(-1.0052181) q[1];
sx q[1];
rz(-3.0977071) q[1];
x q[2];
rz(-1.6789867) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(-0.070852208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(-0.016816703) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.3341981) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.461207) q[0];
sx q[0];
rz(-1.2347504) q[0];
sx q[0];
rz(-1.3075605) q[0];
x q[1];
rz(2.07431) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(-0.50349456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15021819) q[1];
sx q[1];
rz(-0.43979904) q[1];
sx q[1];
rz(-0.99779731) q[1];
rz(-pi) q[2];
rz(-1.2905144) q[3];
sx q[3];
rz(-2.9119769) q[3];
sx q[3];
rz(1.5748101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-0.71643913) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(0.72296468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82397126) q[0];
sx q[0];
rz(-2.9569607) q[0];
sx q[0];
rz(-2.1872107) q[0];
x q[1];
rz(1.1935913) q[2];
sx q[2];
rz(-1.4133487) q[2];
sx q[2];
rz(2.8709656) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55262676) q[1];
sx q[1];
rz(-0.50682658) q[1];
sx q[1];
rz(-1.456702) q[1];
x q[2];
rz(-2.74182) q[3];
sx q[3];
rz(-0.82953605) q[3];
sx q[3];
rz(1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(-2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-2.0422968) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(0.90080558) q[2];
sx q[2];
rz(-1.9566036) q[2];
sx q[2];
rz(-1.4835139) q[2];
rz(2.8156149) q[3];
sx q[3];
rz(-1.3578284) q[3];
sx q[3];
rz(1.1708543) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
