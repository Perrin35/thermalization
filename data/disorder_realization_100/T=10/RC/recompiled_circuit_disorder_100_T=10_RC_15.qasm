OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(2.7913845) q[0];
sx q[0];
rz(12.933001) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769387) q[0];
sx q[0];
rz(-1.9908449) q[0];
sx q[0];
rz(-2.7647892) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8800814) q[2];
sx q[2];
rz(-0.72578428) q[2];
sx q[2];
rz(2.6860565) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9283596) q[1];
sx q[1];
rz(-1.302049) q[1];
sx q[1];
rz(1.0006204) q[1];
rz(1.5815758) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(-0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.4166819) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.5753997) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5088168) q[0];
sx q[0];
rz(-0.59558376) q[0];
sx q[0];
rz(-2.3217208) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86512489) q[2];
sx q[2];
rz(-1.7127617) q[2];
sx q[2];
rz(2.1324468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.054489) q[1];
sx q[1];
rz(-2.7627355) q[1];
sx q[1];
rz(-2.8701251) q[1];
rz(2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(1.1225835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(-2.5235126) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(-1.9504257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55710775) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(2.0273897) q[0];
rz(-pi) q[1];
rz(-0.70706681) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(0.51800874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35347963) q[1];
sx q[1];
rz(-0.2971) q[1];
sx q[1];
rz(2.3658845) q[1];
x q[2];
rz(-1.1690087) q[3];
sx q[3];
rz(-1.0649293) q[3];
sx q[3];
rz(1.5798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(2.1276786) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(-1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-2.9372835) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2134593) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(0.34780985) q[0];
rz(2.8795305) q[2];
sx q[2];
rz(-0.045496551) q[2];
sx q[2];
rz(-0.87693518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7912485) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.1715602) q[1];
rz(0.022425671) q[3];
sx q[3];
rz(-0.86655819) q[3];
sx q[3];
rz(1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4541645) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(2.6320809) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(2.1557504) q[1];
sx q[1];
rz(-1.3580094) q[1];
sx q[1];
rz(-1.3607508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96201345) q[0];
sx q[0];
rz(-2.8498581) q[0];
sx q[0];
rz(2.4289262) q[0];
rz(-pi) q[1];
rz(2.3809793) q[2];
sx q[2];
rz(-2.0509655) q[2];
sx q[2];
rz(-1.7825356) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4701925) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(0.43882521) q[1];
rz(-pi) q[2];
rz(-3.0642062) q[3];
sx q[3];
rz(-2.1688528) q[3];
sx q[3];
rz(-0.070904562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(-0.45587513) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(-1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.8967569) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21101418) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(-0.37822322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2116824) q[2];
sx q[2];
rz(-1.1615796) q[2];
sx q[2];
rz(-0.90604679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6323159) q[1];
sx q[1];
rz(-2.5264611) q[1];
sx q[1];
rz(-2.5110911) q[1];
x q[2];
rz(-0.95008696) q[3];
sx q[3];
rz(-1.6380777) q[3];
sx q[3];
rz(0.34611191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-0.45864027) q[2];
rz(2.4300872) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28073072) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.77805) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(0.71969676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62455432) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(-2.4863003) q[0];
rz(-pi) q[1];
rz(2.3709488) q[2];
sx q[2];
rz(-1.0174123) q[2];
sx q[2];
rz(1.9945952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2232659) q[1];
sx q[1];
rz(-2.1152788) q[1];
sx q[1];
rz(1.5189927) q[1];
rz(-pi) q[2];
rz(1.3004488) q[3];
sx q[3];
rz(-0.90793228) q[3];
sx q[3];
rz(2.8519252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65934962) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(-1.6465181) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0197501) q[0];
sx q[0];
rz(-1.9103721) q[0];
sx q[0];
rz(-1.3440488) q[0];
rz(-pi) q[1];
rz(-1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(-1.9538823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9547792) q[1];
sx q[1];
rz(-0.28007945) q[1];
sx q[1];
rz(-0.38444744) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80963482) q[3];
sx q[3];
rz(-1.0757043) q[3];
sx q[3];
rz(-1.9593057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(-1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(-2.263608) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.7810129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1454865) q[0];
sx q[0];
rz(-1.5331475) q[0];
sx q[0];
rz(2.1426175) q[0];
x q[1];
rz(-0.046710308) q[2];
sx q[2];
rz(-0.94418664) q[2];
sx q[2];
rz(1.7358629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18371836) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(0.79224371) q[1];
x q[2];
rz(-0.15501539) q[3];
sx q[3];
rz(-1.985637) q[3];
sx q[3];
rz(2.993194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(-3.0950586) q[0];
rz(-0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-0.7199026) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324632) q[0];
sx q[0];
rz(-1.6502408) q[0];
sx q[0];
rz(-1.3264873) q[0];
x q[1];
rz(1.5762395) q[2];
sx q[2];
rz(-1.8796225) q[2];
sx q[2];
rz(0.18143166) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9626999) q[1];
sx q[1];
rz(-1.8176515) q[1];
sx q[1];
rz(-0.015951338) q[1];
rz(-0.87749691) q[3];
sx q[3];
rz(-2.174456) q[3];
sx q[3];
rz(0.391215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-0.36454501) q[2];
sx q[2];
rz(-2.400763) q[2];
sx q[2];
rz(2.7999067) q[2];
rz(-3.1090267) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];