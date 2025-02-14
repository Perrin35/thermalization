OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(2.252993) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(-1.3808274) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91869575) q[0];
sx q[0];
rz(-1.7641506) q[0];
sx q[0];
rz(1.7599289) q[0];
rz(2.1031453) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(2.7388245) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9819698) q[1];
sx q[1];
rz(-1.2352984) q[1];
sx q[1];
rz(-1.5606176) q[1];
x q[2];
rz(1.0303792) q[3];
sx q[3];
rz(-0.4183397) q[3];
sx q[3];
rz(-2.6650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8218653) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(0.099893959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(-1.2167759) q[0];
rz(-1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(-2.7144576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1841284) q[0];
sx q[0];
rz(-1.9165358) q[0];
sx q[0];
rz(-0.60681245) q[0];
rz(-0.92556503) q[2];
sx q[2];
rz(-1.3750018) q[2];
sx q[2];
rz(-0.93709556) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1698902) q[1];
sx q[1];
rz(-2.8198543) q[1];
sx q[1];
rz(-2.02969) q[1];
rz(1.7747545) q[3];
sx q[3];
rz(-1.5195623) q[3];
sx q[3];
rz(-0.27328396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0624258) q[2];
sx q[2];
rz(-1.7288952) q[2];
sx q[2];
rz(1.3986577) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(-1.5980501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(1.1687763) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(2.5659836) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.19307) q[0];
sx q[0];
rz(-2.2098477) q[0];
sx q[0];
rz(1.5568514) q[0];
rz(-3.124442) q[2];
sx q[2];
rz(-2.2313925) q[2];
sx q[2];
rz(-1.2397546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.528426) q[1];
sx q[1];
rz(-1.7878782) q[1];
sx q[1];
rz(2.1291127) q[1];
x q[2];
rz(-3.1166385) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(0.11634532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1316954) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(-2.1843145) q[2];
rz(2.5668528) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(-1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-2.1888581) q[0];
sx q[0];
rz(-1.2611457) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(-1.7877158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1607504) q[0];
sx q[0];
rz(-2.7880362) q[0];
sx q[0];
rz(-2.4048231) q[0];
rz(0.48575966) q[2];
sx q[2];
rz(-1.7339212) q[2];
sx q[2];
rz(-2.8201134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9685843) q[1];
sx q[1];
rz(-0.49793303) q[1];
sx q[1];
rz(-3.0281316) q[1];
x q[2];
rz(0.98815496) q[3];
sx q[3];
rz(-2.372962) q[3];
sx q[3];
rz(2.5386794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42133078) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(-1.2565695) q[2];
rz(-2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(-0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.310815) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-0.34314439) q[0];
rz(0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.7247346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2767216) q[0];
sx q[0];
rz(-1.4007995) q[0];
sx q[0];
rz(1.8525339) q[0];
rz(-pi) q[1];
rz(-1.2447276) q[2];
sx q[2];
rz(-2.2448501) q[2];
sx q[2];
rz(-1.3671966) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.062202104) q[1];
sx q[1];
rz(-2.6401527) q[1];
sx q[1];
rz(-1.8556684) q[1];
rz(1.0846109) q[3];
sx q[3];
rz(-0.47322464) q[3];
sx q[3];
rz(1.5718232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5938277) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(-2.055638) q[2];
rz(-2.2222399) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74349657) q[0];
sx q[0];
rz(-1.9323823) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(-2.4010557) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(-0.83121306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75126264) q[0];
sx q[0];
rz(-2.7075276) q[0];
sx q[0];
rz(1.0573666) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1341427) q[2];
sx q[2];
rz(-1.2937355) q[2];
sx q[2];
rz(-1.7817093) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1676291) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(-2.7401398) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.020757787) q[3];
sx q[3];
rz(-1.100173) q[3];
sx q[3];
rz(0.94604674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0703766) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(-2.0224723) q[2];
rz(-1.8708771) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(-1.5420089) q[0];
rz(1.0150602) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-0.17280811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8625582) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(0.59086694) q[0];
rz(-pi) q[1];
rz(-0.68384275) q[2];
sx q[2];
rz(-1.8617861) q[2];
sx q[2];
rz(0.061372193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0902675) q[1];
sx q[1];
rz(-2.6410612) q[1];
sx q[1];
rz(-1.0330908) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0181581) q[3];
sx q[3];
rz(-0.88200906) q[3];
sx q[3];
rz(1.3017201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.479594) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(0.97770989) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(1.9112816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8136895) q[0];
sx q[0];
rz(-2.8459025) q[0];
sx q[0];
rz(-0.015722474) q[0];
rz(-1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(1.9421008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5003842) q[0];
sx q[0];
rz(-0.50725021) q[0];
sx q[0];
rz(2.0709867) q[0];
rz(-pi) q[1];
rz(0.052414465) q[2];
sx q[2];
rz(-2.0611066) q[2];
sx q[2];
rz(2.3971967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75592985) q[1];
sx q[1];
rz(-1.0781059) q[1];
sx q[1];
rz(2.1098968) q[1];
rz(1.1114798) q[3];
sx q[3];
rz(-0.62255854) q[3];
sx q[3];
rz(0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1477995) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.7564868) q[2];
rz(-2.5177453) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(-2.4483335) q[0];
rz(2.8765826) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(1.5230491) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023848195) q[0];
sx q[0];
rz(-1.5818412) q[0];
sx q[0];
rz(-2.3680229) q[0];
rz(2.1586645) q[2];
sx q[2];
rz(-2.1097906) q[2];
sx q[2];
rz(-1.7807775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0761145) q[1];
sx q[1];
rz(-1.5468742) q[1];
sx q[1];
rz(-2.8019399) q[1];
rz(0.070793666) q[3];
sx q[3];
rz(-0.99081836) q[3];
sx q[3];
rz(1.9656612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30274621) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(1.4576853) q[2];
rz(-2.1863106) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(0.83520755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569698) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(2.6589174) q[0];
rz(0.92974281) q[1];
sx q[1];
rz(-2.1371806) q[1];
sx q[1];
rz(-0.39628705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549462) q[0];
sx q[0];
rz(-0.56342331) q[0];
sx q[0];
rz(2.4231829) q[0];
rz(1.4192788) q[2];
sx q[2];
rz(-1.7477525) q[2];
sx q[2];
rz(-0.78864509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0490108) q[1];
sx q[1];
rz(-1.3900458) q[1];
sx q[1];
rz(0.84047079) q[1];
x q[2];
rz(-1.9307053) q[3];
sx q[3];
rz(-1.1904089) q[3];
sx q[3];
rz(-1.701327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(0.3271884) q[2];
rz(-0.78091019) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0384211) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(0.054389537) q[2];
sx q[2];
rz(-2.4475606) q[2];
sx q[2];
rz(-2.9771752) q[2];
rz(-0.38402186) q[3];
sx q[3];
rz(-2.1441318) q[3];
sx q[3];
rz(-1.287788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
