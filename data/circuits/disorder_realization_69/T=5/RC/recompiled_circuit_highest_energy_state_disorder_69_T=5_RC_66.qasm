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
rz(-2.7954571) q[0];
sx q[0];
rz(-1.6383642) q[0];
sx q[0];
rz(2.4099953) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(-0.43391689) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.159891) q[0];
sx q[0];
rz(-2.4004687) q[0];
sx q[0];
rz(-1.6228049) q[0];
x q[1];
rz(-2.3375422) q[2];
sx q[2];
rz(-0.75166273) q[2];
sx q[2];
rz(2.5676905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.078339442) q[1];
sx q[1];
rz(-1.1800449) q[1];
sx q[1];
rz(-2.2558252) q[1];
x q[2];
rz(-3.1103575) q[3];
sx q[3];
rz(-1.3797576) q[3];
sx q[3];
rz(-1.4174549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4291541) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(-2.1377371) q[2];
rz(2.6058274) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-3.1404449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67343229) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(0.72714192) q[0];
rz(0.06761059) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(-2.0179857) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.985874) q[0];
sx q[0];
rz(-0.56502461) q[0];
sx q[0];
rz(-2.5955517) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6933762) q[2];
sx q[2];
rz(-1.5059587) q[2];
sx q[2];
rz(1.4139869) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3930291) q[1];
sx q[1];
rz(-1.2911699) q[1];
sx q[1];
rz(-0.10348998) q[1];
x q[2];
rz(2.7026342) q[3];
sx q[3];
rz(-2.0062796) q[3];
sx q[3];
rz(1.2964695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2850538) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(-2.6598568) q[2];
rz(2.5887865) q[3];
sx q[3];
rz(-2.063664) q[3];
sx q[3];
rz(3.0520181) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32473096) q[0];
sx q[0];
rz(-0.4751927) q[0];
sx q[0];
rz(-3.0480296) q[0];
rz(1.7612673) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(-0.63724744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475315) q[0];
sx q[0];
rz(-2.0009181) q[0];
sx q[0];
rz(2.568666) q[0];
x q[1];
rz(1.3445417) q[2];
sx q[2];
rz(-0.97469869) q[2];
sx q[2];
rz(1.4486194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5384851) q[1];
sx q[1];
rz(-0.39834312) q[1];
sx q[1];
rz(1.6377691) q[1];
rz(-pi) q[2];
rz(-1.6138541) q[3];
sx q[3];
rz(-1.2803439) q[3];
sx q[3];
rz(1.8990979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9351585) q[2];
sx q[2];
rz(-2.5863402) q[2];
sx q[2];
rz(2.4578102) q[2];
rz(0.23009662) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(0.42207119) q[3];
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
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0919331) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(-1.7012713) q[0];
rz(0.47538844) q[1];
sx q[1];
rz(-0.63440591) q[1];
sx q[1];
rz(2.5604274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064539486) q[0];
sx q[0];
rz(-0.7270027) q[0];
sx q[0];
rz(-1.6995656) q[0];
rz(-pi) q[1];
rz(2.0175242) q[2];
sx q[2];
rz(-0.5047732) q[2];
sx q[2];
rz(3.1044132) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6789279) q[1];
sx q[1];
rz(-1.8724672) q[1];
sx q[1];
rz(0.48120166) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.851296) q[3];
sx q[3];
rz(-2.6321342) q[3];
sx q[3];
rz(-0.51028937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0541957) q[2];
sx q[2];
rz(-2.9554458) q[2];
sx q[2];
rz(-0.52097121) q[2];
rz(-1.4969131) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(-1.6269256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5882551) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(-2.0113373) q[0];
rz(2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(-2.030453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53752995) q[0];
sx q[0];
rz(-1.1141234) q[0];
sx q[0];
rz(-0.18484331) q[0];
rz(-pi) q[1];
rz(3.0658998) q[2];
sx q[2];
rz(-2.6106129) q[2];
sx q[2];
rz(0.39226433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4823522) q[1];
sx q[1];
rz(-0.9352881) q[1];
sx q[1];
rz(2.9457432) q[1];
rz(-pi) q[2];
rz(-2.2419315) q[3];
sx q[3];
rz(-2.3658345) q[3];
sx q[3];
rz(-0.44246261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40372103) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(-1.0408638) q[2];
rz(-2.3216085) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(0.6238873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2417004) q[0];
sx q[0];
rz(-2.5427759) q[0];
sx q[0];
rz(2.4851121) q[0];
rz(1.7680602) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(2.0097282) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0431932) q[0];
sx q[0];
rz(-0.88623673) q[0];
sx q[0];
rz(1.8266023) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0229576) q[2];
sx q[2];
rz(-2.1716087) q[2];
sx q[2];
rz(-3.09077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.060185) q[1];
sx q[1];
rz(-2.2814732) q[1];
sx q[1];
rz(-3.1050472) q[1];
x q[2];
rz(-1.5664212) q[3];
sx q[3];
rz(-1.4203912) q[3];
sx q[3];
rz(-2.7657417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6953096) q[2];
sx q[2];
rz(-0.60574564) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(-0.22932912) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730154) q[0];
sx q[0];
rz(-1.9972766) q[0];
sx q[0];
rz(-1.2314433) q[0];
rz(3.0629509) q[1];
sx q[1];
rz(-1.4631203) q[1];
sx q[1];
rz(0.15422779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.70792) q[0];
sx q[0];
rz(-1.7485973) q[0];
sx q[0];
rz(2.7975797) q[0];
x q[1];
rz(-2.8381589) q[2];
sx q[2];
rz(-1.1772441) q[2];
sx q[2];
rz(2.9224666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63111607) q[1];
sx q[1];
rz(-1.0923903) q[1];
sx q[1];
rz(-1.1523184) q[1];
rz(-2.6811872) q[3];
sx q[3];
rz(-2.6889927) q[3];
sx q[3];
rz(2.3456338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8057033) q[2];
sx q[2];
rz(-1.5259537) q[2];
sx q[2];
rz(1.6501144) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(-1.5707877) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042628057) q[0];
sx q[0];
rz(-2.0476116) q[0];
sx q[0];
rz(1.4549103) q[0];
rz(2.1658354) q[1];
sx q[1];
rz(-2.0819596) q[1];
sx q[1];
rz(1.8029433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2678626) q[0];
sx q[0];
rz(-2.3288245) q[0];
sx q[0];
rz(2.1186747) q[0];
rz(-pi) q[1];
rz(-1.067131) q[2];
sx q[2];
rz(-1.2758453) q[2];
sx q[2];
rz(-0.3857715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2690791) q[1];
sx q[1];
rz(-2.1912406) q[1];
sx q[1];
rz(-0.2618813) q[1];
x q[2];
rz(1.1986794) q[3];
sx q[3];
rz(-0.94148472) q[3];
sx q[3];
rz(-0.019364186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6300388) q[2];
sx q[2];
rz(-1.0147107) q[2];
sx q[2];
rz(2.3021728) q[2];
rz(-1.2342341) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3473174) q[0];
sx q[0];
rz(-1.3823771) q[0];
sx q[0];
rz(-0.41912249) q[0];
rz(1.4979111) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(2.7083414) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.142665) q[0];
sx q[0];
rz(-1.4523399) q[0];
sx q[0];
rz(-2.0303557) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.584132) q[2];
sx q[2];
rz(-2.2277701) q[2];
sx q[2];
rz(-1.8767534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8318787) q[1];
sx q[1];
rz(-2.0818748) q[1];
sx q[1];
rz(-2.8161777) q[1];
rz(-2.8178704) q[3];
sx q[3];
rz(-0.60327134) q[3];
sx q[3];
rz(2.3761185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7822632) q[2];
sx q[2];
rz(-1.9908345) q[2];
sx q[2];
rz(-2.9456054) q[2];
rz(-2.0187142) q[3];
sx q[3];
rz(-0.71176088) q[3];
sx q[3];
rz(-1.6897197) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48542431) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(-2.0458903) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(0.21024545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4304786) q[0];
sx q[0];
rz(-1.0000044) q[0];
sx q[0];
rz(2.9080703) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0265507) q[2];
sx q[2];
rz(-0.91918531) q[2];
sx q[2];
rz(-2.5082805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89344674) q[1];
sx q[1];
rz(-0.67091828) q[1];
sx q[1];
rz(-2.4144961) q[1];
x q[2];
rz(-0.59981646) q[3];
sx q[3];
rz(-1.2251523) q[3];
sx q[3];
rz(-1.7315854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11253396) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(-1.0851592) q[2];
rz(-2.8339913) q[3];
sx q[3];
rz(-0.31809536) q[3];
sx q[3];
rz(-2.4881081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45162421) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(2.4901509) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(-0.015301558) q[2];
sx q[2];
rz(-0.87991426) q[2];
sx q[2];
rz(2.8755398) q[2];
rz(-2.2023946) q[3];
sx q[3];
rz(-1.9892577) q[3];
sx q[3];
rz(1.6792959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
