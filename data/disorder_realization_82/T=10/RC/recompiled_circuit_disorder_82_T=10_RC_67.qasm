OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(-1.0441138) q[1];
sx q[1];
rz(1.27966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188388) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(1.3214146) q[0];
rz(1.3470634) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.3734693) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3763334) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(2.5024662) q[1];
rz(-2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(1.83028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0063643) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(2.6541236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23232846) q[0];
sx q[0];
rz(-0.10350138) q[0];
sx q[0];
rz(0.89719015) q[0];
x q[1];
rz(-2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(-2.2573059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5408052) q[1];
sx q[1];
rz(-1.729319) q[1];
sx q[1];
rz(-1.7608789) q[1];
rz(-pi) q[2];
rz(-1.8879714) q[3];
sx q[3];
rz(-1.653169) q[3];
sx q[3];
rz(1.9272643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.1797103) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-3.1029207) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4653141) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(-2.7200384) q[0];
rz(-pi) q[1];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(-0.97937102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-2.5249081) q[1];
sx q[1];
rz(-0.1174121) q[1];
x q[2];
rz(2.1915216) q[3];
sx q[3];
rz(-2.7347703) q[3];
sx q[3];
rz(-2.1086958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(2.0866701) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(-0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(2.2213675) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738149) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(-0.52315229) q[0];
x q[1];
rz(-0.22569457) q[2];
sx q[2];
rz(-0.36964551) q[2];
sx q[2];
rz(0.29633488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0238266) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(1.0002506) q[1];
rz(-pi) q[2];
rz(1.7761049) q[3];
sx q[3];
rz(-1.9234386) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(-1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(-0.76830307) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950726) q[0];
sx q[0];
rz(-3.0591024) q[0];
sx q[0];
rz(-2.0399658) q[0];
rz(-1.808851) q[2];
sx q[2];
rz(-2.9946054) q[2];
sx q[2];
rz(-0.57839314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0799775) q[1];
sx q[1];
rz(-1.6528168) q[1];
sx q[1];
rz(-1.7916937) q[1];
x q[2];
rz(2.1702607) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0129464) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(-2.9034555) q[0];
rz(2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(1.628081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742764) q[0];
sx q[0];
rz(-0.58642504) q[0];
sx q[0];
rz(-0.61074722) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.020535) q[2];
sx q[2];
rz(-1.626484) q[2];
sx q[2];
rz(0.72292098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(0.97126295) q[1];
rz(-0.87168872) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(-0.20760078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.9546753) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163527) q[0];
sx q[0];
rz(-1.2901187) q[0];
sx q[0];
rz(-1.0513845) q[0];
x q[1];
rz(0.45849623) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(2.3022848) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0901511) q[1];
sx q[1];
rz(-1.7045583) q[1];
sx q[1];
rz(-2.306288) q[1];
rz(-2.7508221) q[3];
sx q[3];
rz(-2.4739389) q[3];
sx q[3];
rz(2.8323176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223406) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-2.2858455) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(1.0303248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1925826) q[0];
sx q[0];
rz(-2.2826676) q[0];
sx q[0];
rz(-2.683995) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6542997) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(0.63937843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47146637) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(-2.9832277) q[1];
x q[2];
rz(0.86352591) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(-1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90056706) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(2.7450949) q[0];
rz(-3.1124434) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(2.0707891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70442048) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(-2.0541595) q[1];
rz(-pi) q[2];
rz(-3.1177403) q[3];
sx q[3];
rz(-0.94787041) q[3];
sx q[3];
rz(2.9305487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(0.35183516) q[0];
rz(0.31967638) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-0.19616729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88403945) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(2.2474225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9253795) q[2];
sx q[2];
rz(-0.14419975) q[2];
sx q[2];
rz(-1.9229638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13388053) q[1];
sx q[1];
rz(-2.4073283) q[1];
sx q[1];
rz(-1.5937362) q[1];
rz(-0.65532834) q[3];
sx q[3];
rz(-0.95643759) q[3];
sx q[3];
rz(-2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-0.081136726) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-0.20028533) q[2];
sx q[2];
rz(-0.23858783) q[2];
sx q[2];
rz(-1.5993652) q[2];
rz(0.046394596) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];