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
rz(-2.7201535) q[0];
sx q[0];
rz(-1.8869737) q[0];
sx q[0];
rz(-0.52809554) q[0];
rz(2.7948607) q[1];
sx q[1];
rz(-1.3277227) q[1];
sx q[1];
rz(-2.858207) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91862667) q[0];
sx q[0];
rz(-0.98969995) q[0];
sx q[0];
rz(-1.6870935) q[0];
rz(-1.07017) q[2];
sx q[2];
rz(-2.3267022) q[2];
sx q[2];
rz(1.7931149) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49453255) q[1];
sx q[1];
rz(-1.8683829) q[1];
sx q[1];
rz(1.3045761) q[1];
x q[2];
rz(-0.10589604) q[3];
sx q[3];
rz(-1.9710324) q[3];
sx q[3];
rz(-0.36508158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9458719) q[2];
sx q[2];
rz(-1.4737782) q[2];
sx q[2];
rz(-2.8742068) q[2];
rz(-2.9684559) q[3];
sx q[3];
rz(-2.0045529) q[3];
sx q[3];
rz(-0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6807569) q[0];
sx q[0];
rz(-0.84592485) q[0];
sx q[0];
rz(-2.8435006) q[0];
rz(-0.2671034) q[1];
sx q[1];
rz(-2.2593081) q[1];
sx q[1];
rz(2.2411236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33997275) q[0];
sx q[0];
rz(-0.34944926) q[0];
sx q[0];
rz(0.79876113) q[0];
x q[1];
rz(2.3411269) q[2];
sx q[2];
rz(-2.9334407) q[2];
sx q[2];
rz(0.58074927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5067277) q[1];
sx q[1];
rz(-1.3758476) q[1];
sx q[1];
rz(-0.69201236) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44690981) q[3];
sx q[3];
rz(-0.70820184) q[3];
sx q[3];
rz(0.59659905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.71821) q[2];
sx q[2];
rz(-2.155828) q[2];
sx q[2];
rz(0.59744376) q[2];
rz(2.1256223) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(1.5731251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628767) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(0.065091982) q[0];
rz(2.0386631) q[1];
sx q[1];
rz(-1.8850336) q[1];
sx q[1];
rz(-2.7535313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86923118) q[0];
sx q[0];
rz(-2.3359406) q[0];
sx q[0];
rz(-2.8153573) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36696649) q[2];
sx q[2];
rz(-1.6410368) q[2];
sx q[2];
rz(-2.5052469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2414723) q[1];
sx q[1];
rz(-0.44229506) q[1];
sx q[1];
rz(-1.1957748) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8592836) q[3];
sx q[3];
rz(-1.3460717) q[3];
sx q[3];
rz(-2.2109179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97543657) q[2];
sx q[2];
rz(-2.1076951) q[2];
sx q[2];
rz(3.0530829) q[2];
rz(-0.65781188) q[3];
sx q[3];
rz(-0.43840539) q[3];
sx q[3];
rz(-1.2737761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.45519644) q[0];
sx q[0];
rz(-0.96298591) q[0];
sx q[0];
rz(2.2121867) q[0];
rz(-1.1691673) q[1];
sx q[1];
rz(-0.87375748) q[1];
sx q[1];
rz(0.1263617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39462659) q[0];
sx q[0];
rz(-1.639059) q[0];
sx q[0];
rz(1.9771876) q[0];
rz(-1.7772555) q[2];
sx q[2];
rz(-1.2480309) q[2];
sx q[2];
rz(-1.9545123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.026638) q[1];
sx q[1];
rz(-1.8787787) q[1];
sx q[1];
rz(-2.6796209) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8245623) q[3];
sx q[3];
rz(-1.9477748) q[3];
sx q[3];
rz(2.1351817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8440642) q[2];
sx q[2];
rz(-2.0487831) q[2];
sx q[2];
rz(-1.4275985) q[2];
rz(1.1045688) q[3];
sx q[3];
rz(-1.5966871) q[3];
sx q[3];
rz(-0.67232084) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70668689) q[0];
sx q[0];
rz(-0.66449419) q[0];
sx q[0];
rz(-1.2427166) q[0];
rz(0.13936123) q[1];
sx q[1];
rz(-0.31529537) q[1];
sx q[1];
rz(0.79162663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2682429) q[0];
sx q[0];
rz(-1.5333383) q[0];
sx q[0];
rz(1.7245049) q[0];
rz(1.3940195) q[2];
sx q[2];
rz(-1.0097754) q[2];
sx q[2];
rz(-2.0732806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6649896) q[1];
sx q[1];
rz(-1.2634189) q[1];
sx q[1];
rz(3.1165131) q[1];
rz(-pi) q[2];
rz(3.0221927) q[3];
sx q[3];
rz(-2.1704333) q[3];
sx q[3];
rz(1.3493845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7255154) q[2];
sx q[2];
rz(-1.6107586) q[2];
sx q[2];
rz(1.3610972) q[2];
rz(2.1040037) q[3];
sx q[3];
rz(-1.6890084) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1182275) q[0];
sx q[0];
rz(-0.26015493) q[0];
sx q[0];
rz(-2.9628229) q[0];
rz(-2.692692) q[1];
sx q[1];
rz(-1.1652378) q[1];
sx q[1];
rz(1.9662205) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7473874) q[0];
sx q[0];
rz(-0.76524599) q[0];
sx q[0];
rz(-0.09455643) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.825338) q[2];
sx q[2];
rz(-2.3582728) q[2];
sx q[2];
rz(-2.6412727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11131249) q[1];
sx q[1];
rz(-2.2759109) q[1];
sx q[1];
rz(2.6837803) q[1];
x q[2];
rz(1.337758) q[3];
sx q[3];
rz(-1.0653598) q[3];
sx q[3];
rz(-2.4580372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7281404) q[2];
sx q[2];
rz(-2.7588625) q[2];
sx q[2];
rz(2.9276796) q[2];
rz(1.7119857) q[3];
sx q[3];
rz(-1.6702646) q[3];
sx q[3];
rz(1.8179998) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0578617) q[0];
sx q[0];
rz(-1.2849176) q[0];
sx q[0];
rz(-0.65823746) q[0];
rz(-1.4312875) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(-1.5586982) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76097902) q[0];
sx q[0];
rz(-1.7581994) q[0];
sx q[0];
rz(2.6949469) q[0];
x q[1];
rz(-0.30679254) q[2];
sx q[2];
rz(-1.3500096) q[2];
sx q[2];
rz(0.38635269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.538529) q[1];
sx q[1];
rz(-2.6683806) q[1];
sx q[1];
rz(2.3296001) q[1];
rz(-1.962108) q[3];
sx q[3];
rz(-1.0265304) q[3];
sx q[3];
rz(2.5044092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3805286) q[2];
sx q[2];
rz(-1.5498127) q[2];
sx q[2];
rz(-1.9350249) q[2];
rz(-2.1207464) q[3];
sx q[3];
rz(-0.51488334) q[3];
sx q[3];
rz(-0.90768874) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9970488) q[0];
sx q[0];
rz(-0.21268614) q[0];
sx q[0];
rz(-2.8295243) q[0];
rz(2.8128305) q[1];
sx q[1];
rz(-1.5212675) q[1];
sx q[1];
rz(0.75334466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8341136) q[0];
sx q[0];
rz(-0.8571107) q[0];
sx q[0];
rz(0.77425787) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0814799) q[2];
sx q[2];
rz(-1.5897255) q[2];
sx q[2];
rz(-0.14046803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6623913) q[1];
sx q[1];
rz(-2.0298784) q[1];
sx q[1];
rz(0.69028141) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83266074) q[3];
sx q[3];
rz(-0.7023905) q[3];
sx q[3];
rz(-0.47846068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.83851492) q[2];
sx q[2];
rz(-1.3884156) q[2];
sx q[2];
rz(0.77084368) q[2];
rz(-0.54909697) q[3];
sx q[3];
rz(-1.2884459) q[3];
sx q[3];
rz(-2.2652266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751806) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(-1.3394248) q[0];
rz(-2.0088947) q[1];
sx q[1];
rz(-2.7481672) q[1];
sx q[1];
rz(0.045230953) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6409174) q[0];
sx q[0];
rz(-1.7637327) q[0];
sx q[0];
rz(2.0125602) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2140177) q[2];
sx q[2];
rz(-1.456481) q[2];
sx q[2];
rz(3.1014991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69581) q[1];
sx q[1];
rz(-2.8572866) q[1];
sx q[1];
rz(0.48281702) q[1];
rz(1.0186152) q[3];
sx q[3];
rz(-1.3613627) q[3];
sx q[3];
rz(-2.6554573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.905978) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(0.81602412) q[2];
rz(1.3927381) q[3];
sx q[3];
rz(-1.1002898) q[3];
sx q[3];
rz(-1.6284774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29869646) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(-1.5604875) q[0];
rz(-0.16079482) q[1];
sx q[1];
rz(-2.000688) q[1];
sx q[1];
rz(-2.2198832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.164249) q[0];
sx q[0];
rz(-1.5613587) q[0];
sx q[0];
rz(1.759985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99803136) q[2];
sx q[2];
rz(-1.0596399) q[2];
sx q[2];
rz(-2.9645408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8556577) q[1];
sx q[1];
rz(-1.8536356) q[1];
sx q[1];
rz(-0.15423473) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5215048) q[3];
sx q[3];
rz(-0.26066565) q[3];
sx q[3];
rz(1.0083782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2738652) q[2];
sx q[2];
rz(-0.92550698) q[2];
sx q[2];
rz(-1.072139) q[2];
rz(-2.9169361) q[3];
sx q[3];
rz(-1.6499237) q[3];
sx q[3];
rz(-0.91607654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199725) q[0];
sx q[0];
rz(-2.2157123) q[0];
sx q[0];
rz(-2.8257688) q[0];
rz(-0.074180457) q[1];
sx q[1];
rz(-2.0769495) q[1];
sx q[1];
rz(3.1202797) q[1];
rz(1.5431719) q[2];
sx q[2];
rz(-2.3541214) q[2];
sx q[2];
rz(0.036938112) q[2];
rz(-0.76618054) q[3];
sx q[3];
rz(-1.8591559) q[3];
sx q[3];
rz(0.55606203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
