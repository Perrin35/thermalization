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
rz(-2.4969555) q[0];
sx q[0];
rz(-0.58965373) q[0];
sx q[0];
rz(-1.0876422) q[0];
rz(-0.015406869) q[1];
sx q[1];
rz(3.5313731) q[1];
sx q[1];
rz(11.402147) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22237865) q[0];
sx q[0];
rz(-1.5660677) q[0];
sx q[0];
rz(1.7008002) q[0];
rz(1.2306326) q[2];
sx q[2];
rz(-0.90848604) q[2];
sx q[2];
rz(1.8669805) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77736311) q[1];
sx q[1];
rz(-2.573488) q[1];
sx q[1];
rz(1.9757759) q[1];
rz(0.075957493) q[3];
sx q[3];
rz(-1.6884558) q[3];
sx q[3];
rz(-0.75353384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9001793) q[2];
sx q[2];
rz(-2.5799077) q[2];
sx q[2];
rz(-0.64603311) q[2];
rz(-0.25520405) q[3];
sx q[3];
rz(-2.6845158) q[3];
sx q[3];
rz(-1.7469762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141041) q[0];
sx q[0];
rz(-2.8059967) q[0];
sx q[0];
rz(-2.3345729) q[0];
rz(1.457816) q[1];
sx q[1];
rz(-0.57676637) q[1];
sx q[1];
rz(1.797765) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3407077) q[0];
sx q[0];
rz(-1.721764) q[0];
sx q[0];
rz(2.8503391) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97118158) q[2];
sx q[2];
rz(-0.65752013) q[2];
sx q[2];
rz(-0.91577851) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8763153) q[1];
sx q[1];
rz(-2.1723299) q[1];
sx q[1];
rz(-0.6024036) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5206512) q[3];
sx q[3];
rz(-0.72948958) q[3];
sx q[3];
rz(1.4514841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.9606216) q[2];
sx q[2];
rz(-1.0701067) q[2];
sx q[2];
rz(0.19372678) q[2];
rz(-1.6173897) q[3];
sx q[3];
rz(-2.3834855) q[3];
sx q[3];
rz(0.1787506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206361) q[0];
sx q[0];
rz(-2.9036324) q[0];
sx q[0];
rz(-0.97610193) q[0];
rz(-2.2528265) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(-0.1730473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7905549) q[0];
sx q[0];
rz(-1.5332744) q[0];
sx q[0];
rz(0.0028126082) q[0];
x q[1];
rz(-2.9251418) q[2];
sx q[2];
rz(-0.5308639) q[2];
sx q[2];
rz(0.59950874) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0806607) q[1];
sx q[1];
rz(-2.073604) q[1];
sx q[1];
rz(1.7601556) q[1];
rz(-0.32239989) q[3];
sx q[3];
rz(-1.4677257) q[3];
sx q[3];
rz(-0.88503982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.139107) q[2];
sx q[2];
rz(-1.8128914) q[2];
sx q[2];
rz(-0.71387449) q[2];
rz(0.21126963) q[3];
sx q[3];
rz(-2.1043089) q[3];
sx q[3];
rz(-2.5762288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738889) q[0];
sx q[0];
rz(-1.9581032) q[0];
sx q[0];
rz(1.1327889) q[0];
rz(-2.6702113) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(2.7105791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4296161) q[0];
sx q[0];
rz(-1.683748) q[0];
sx q[0];
rz(-0.21265257) q[0];
rz(1.9507031) q[2];
sx q[2];
rz(-1.9089193) q[2];
sx q[2];
rz(-0.1491216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9412028) q[1];
sx q[1];
rz(-1.871849) q[1];
sx q[1];
rz(-2.5832504) q[1];
rz(-pi) q[2];
rz(0.94915821) q[3];
sx q[3];
rz(-1.7931058) q[3];
sx q[3];
rz(-2.3883826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5203349) q[2];
sx q[2];
rz(-0.73953491) q[2];
sx q[2];
rz(0.4062824) q[2];
rz(1.2064365) q[3];
sx q[3];
rz(-2.0483569) q[3];
sx q[3];
rz(2.5549197) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5093812) q[0];
sx q[0];
rz(-2.2659232) q[0];
sx q[0];
rz(-0.32633728) q[0];
rz(-0.95411602) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(3.1180678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81563595) q[0];
sx q[0];
rz(-1.1679497) q[0];
sx q[0];
rz(2.7340552) q[0];
rz(-pi) q[1];
rz(0.40074375) q[2];
sx q[2];
rz(-2.4432428) q[2];
sx q[2];
rz(-1.1782447) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70290804) q[1];
sx q[1];
rz(-2.6881235) q[1];
sx q[1];
rz(2.2681342) q[1];
x q[2];
rz(0.29554772) q[3];
sx q[3];
rz(-1.4369399) q[3];
sx q[3];
rz(-0.16242151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.894459) q[2];
sx q[2];
rz(-2.6829312) q[2];
sx q[2];
rz(-0.66391724) q[2];
rz(-0.89020056) q[3];
sx q[3];
rz(-1.5492487) q[3];
sx q[3];
rz(-0.012705407) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3506055) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(-0.33082333) q[0];
rz(-0.36258969) q[1];
sx q[1];
rz(-1.1793084) q[1];
sx q[1];
rz(-2.6258452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2543751) q[0];
sx q[0];
rz(-1.3332813) q[0];
sx q[0];
rz(-2.919048) q[0];
x q[1];
rz(-0.59402664) q[2];
sx q[2];
rz(-1.7286679) q[2];
sx q[2];
rz(-1.8094339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52236103) q[1];
sx q[1];
rz(-2.4454678) q[1];
sx q[1];
rz(-7*pi/13) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7031851) q[3];
sx q[3];
rz(-0.77927941) q[3];
sx q[3];
rz(-2.8900823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7911239) q[2];
sx q[2];
rz(-1.6914657) q[2];
sx q[2];
rz(-0.80751944) q[2];
rz(3.0430072) q[3];
sx q[3];
rz(-0.8980631) q[3];
sx q[3];
rz(2.0487823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7035141) q[0];
sx q[0];
rz(-2.5683537) q[0];
sx q[0];
rz(-2.692063) q[0];
rz(0.685177) q[1];
sx q[1];
rz(-2.7339869) q[1];
sx q[1];
rz(2.5221241) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0280058) q[0];
sx q[0];
rz(-1.8802693) q[0];
sx q[0];
rz(1.600515) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082498468) q[2];
sx q[2];
rz(-0.52972016) q[2];
sx q[2];
rz(-0.32495435) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7799311) q[1];
sx q[1];
rz(-1.2355991) q[1];
sx q[1];
rz(-0.72104309) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9569217) q[3];
sx q[3];
rz(-1.369056) q[3];
sx q[3];
rz(1.4278442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5770136) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(2.9435834) q[2];
rz(1.7791344) q[3];
sx q[3];
rz(-0.30006108) q[3];
sx q[3];
rz(-0.70039606) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144416) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(1.2855592) q[0];
rz(0.27664912) q[1];
sx q[1];
rz(-1.3686907) q[1];
sx q[1];
rz(3.0329774) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1254107) q[0];
sx q[0];
rz(-1.6526395) q[0];
sx q[0];
rz(-2.0731519) q[0];
x q[1];
rz(3.0622903) q[2];
sx q[2];
rz(-2.1933043) q[2];
sx q[2];
rz(-2.9626737) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6985982) q[1];
sx q[1];
rz(-2.3957771) q[1];
sx q[1];
rz(-1.2520176) q[1];
x q[2];
rz(-3.0158086) q[3];
sx q[3];
rz(-1.1871119) q[3];
sx q[3];
rz(-0.48641274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8331929) q[2];
sx q[2];
rz(-1.9449642) q[2];
sx q[2];
rz(-1.2742554) q[2];
rz(1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(2.275009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891147) q[0];
sx q[0];
rz(-0.72565961) q[0];
sx q[0];
rz(-0.12055483) q[0];
rz(2.4001135) q[1];
sx q[1];
rz(-1.9375786) q[1];
sx q[1];
rz(-3.0659058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2601608) q[0];
sx q[0];
rz(-0.3494702) q[0];
sx q[0];
rz(0.49728877) q[0];
rz(-0.077421256) q[2];
sx q[2];
rz(-1.110437) q[2];
sx q[2];
rz(1.3321239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2774128) q[1];
sx q[1];
rz(-2.2267739) q[1];
sx q[1];
rz(1.6026864) q[1];
rz(-pi) q[2];
rz(-2.5406557) q[3];
sx q[3];
rz(-1.4190201) q[3];
sx q[3];
rz(-2.5545718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.072448298) q[2];
sx q[2];
rz(-2.2587903) q[2];
sx q[2];
rz(-2.7952588) q[2];
rz(2.196178) q[3];
sx q[3];
rz(-1.3225222) q[3];
sx q[3];
rz(0.94714981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2746177) q[0];
sx q[0];
rz(-1.1981542) q[0];
sx q[0];
rz(-2.6937038) q[0];
rz(-0.92944324) q[1];
sx q[1];
rz(-1.7421236) q[1];
sx q[1];
rz(-2.4670752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68088082) q[0];
sx q[0];
rz(-1.370791) q[0];
sx q[0];
rz(2.3671211) q[0];
rz(-1.3849887) q[2];
sx q[2];
rz(-1.4685493) q[2];
sx q[2];
rz(-0.24413689) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0688388) q[1];
sx q[1];
rz(-1.3974032) q[1];
sx q[1];
rz(0.0011464107) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.075427051) q[3];
sx q[3];
rz(-1.9884681) q[3];
sx q[3];
rz(0.088069629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.931539) q[2];
sx q[2];
rz(-1.209582) q[2];
sx q[2];
rz(0.26739576) q[2];
rz(2.0343272) q[3];
sx q[3];
rz(-0.78247726) q[3];
sx q[3];
rz(0.84899181) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7312819) q[0];
sx q[0];
rz(-1.622643) q[0];
sx q[0];
rz(-1.0868764) q[0];
rz(0.80401737) q[1];
sx q[1];
rz(-1.6537279) q[1];
sx q[1];
rz(-2.0860685) q[1];
rz(1.7179391) q[2];
sx q[2];
rz(-1.4999119) q[2];
sx q[2];
rz(-2.9905408) q[2];
rz(3.0576684) q[3];
sx q[3];
rz(-2.2369583) q[3];
sx q[3];
rz(0.53085622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
