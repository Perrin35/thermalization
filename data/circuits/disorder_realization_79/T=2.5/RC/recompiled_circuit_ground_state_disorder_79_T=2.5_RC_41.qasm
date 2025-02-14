OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5496991) q[0];
sx q[0];
rz(-1.1130207) q[0];
sx q[0];
rz(0.23166238) q[0];
rz(-2.7910233) q[1];
sx q[1];
rz(-2.9313593) q[1];
sx q[1];
rz(0.76229131) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1660773) q[0];
sx q[0];
rz(-2.3535801) q[0];
sx q[0];
rz(-2.5322574) q[0];
rz(-2.0341124) q[2];
sx q[2];
rz(-1.1321403) q[2];
sx q[2];
rz(2.2194576) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8232441) q[1];
sx q[1];
rz(-1.6260514) q[1];
sx q[1];
rz(-2.9553508) q[1];
rz(-2.2985994) q[3];
sx q[3];
rz(-2.5268433) q[3];
sx q[3];
rz(-0.45045567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98213696) q[2];
sx q[2];
rz(-0.60672131) q[2];
sx q[2];
rz(0.43178001) q[2];
rz(0.86709705) q[3];
sx q[3];
rz(-2.1853787) q[3];
sx q[3];
rz(-1.6352656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1113488) q[0];
sx q[0];
rz(-0.83006492) q[0];
sx q[0];
rz(-1.1907955) q[0];
rz(-1.2442773) q[1];
sx q[1];
rz(-2.0121274) q[1];
sx q[1];
rz(2.5608889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3235) q[0];
sx q[0];
rz(-0.43002263) q[0];
sx q[0];
rz(-1.2318721) q[0];
rz(-0.9208619) q[2];
sx q[2];
rz(-0.63541401) q[2];
sx q[2];
rz(-2.970545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.001284842) q[1];
sx q[1];
rz(-1.8875889) q[1];
sx q[1];
rz(2.4771734) q[1];
rz(-2.689183) q[3];
sx q[3];
rz(-0.93660855) q[3];
sx q[3];
rz(-0.89434552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0018491) q[2];
sx q[2];
rz(-1.7364343) q[2];
sx q[2];
rz(-1.5025567) q[2];
rz(2.9338845) q[3];
sx q[3];
rz(-1.8307056) q[3];
sx q[3];
rz(-1.0072072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19596066) q[0];
sx q[0];
rz(-2.0646844) q[0];
sx q[0];
rz(0.14863241) q[0];
rz(-1.0736939) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(2.3587842) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40227213) q[0];
sx q[0];
rz(-1.625652) q[0];
sx q[0];
rz(-1.7412454) q[0];
rz(0.53289588) q[2];
sx q[2];
rz(-1.6072465) q[2];
sx q[2];
rz(0.25668555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69082123) q[1];
sx q[1];
rz(-1.3714002) q[1];
sx q[1];
rz(-2.6197207) q[1];
rz(-pi) q[2];
rz(2.3605669) q[3];
sx q[3];
rz(-1.7501004) q[3];
sx q[3];
rz(-0.73338311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9624376) q[2];
sx q[2];
rz(-2.16733) q[2];
sx q[2];
rz(-1.1243189) q[2];
rz(-3.1149241) q[3];
sx q[3];
rz(-2.4494438) q[3];
sx q[3];
rz(0.28287014) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35793316) q[0];
sx q[0];
rz(-1.986035) q[0];
sx q[0];
rz(-1.5643157) q[0];
rz(-1.0464926) q[1];
sx q[1];
rz(-2.107403) q[1];
sx q[1];
rz(-2.0139205) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8822206) q[0];
sx q[0];
rz(-1.5604856) q[0];
sx q[0];
rz(0.59132782) q[0];
rz(-1.2366946) q[2];
sx q[2];
rz(-1.646197) q[2];
sx q[2];
rz(1.4206574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31364031) q[1];
sx q[1];
rz(-1.8871739) q[1];
sx q[1];
rz(2.4183941) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.042759) q[3];
sx q[3];
rz(-1.4526794) q[3];
sx q[3];
rz(1.5631226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6005738) q[2];
sx q[2];
rz(-1.270673) q[2];
sx q[2];
rz(0.18826558) q[2];
rz(-2.6514734) q[3];
sx q[3];
rz(-1.6907588) q[3];
sx q[3];
rz(1.8671487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(1.984943) q[0];
sx q[0];
rz(-1.5771414) q[0];
sx q[0];
rz(-1.4666602) q[0];
rz(2.6690392) q[1];
sx q[1];
rz(-0.87635374) q[1];
sx q[1];
rz(0.98365012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0839) q[0];
sx q[0];
rz(-2.6760457) q[0];
sx q[0];
rz(-2.7280612) q[0];
rz(-0.3717513) q[2];
sx q[2];
rz(-2.6494827) q[2];
sx q[2];
rz(-0.12002698) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90685788) q[1];
sx q[1];
rz(-2.6277797) q[1];
sx q[1];
rz(-1.4545026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9463169) q[3];
sx q[3];
rz(-2.8569479) q[3];
sx q[3];
rz(1.9575021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6321565) q[2];
sx q[2];
rz(-1.3864653) q[2];
sx q[2];
rz(0.40718386) q[2];
rz(1.2716278) q[3];
sx q[3];
rz(-0.41912246) q[3];
sx q[3];
rz(2.4841888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206409) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(-1.2366914) q[0];
rz(-2.0115133) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(2.0225661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73858628) q[0];
sx q[0];
rz(-1.6386014) q[0];
sx q[0];
rz(-0.15371503) q[0];
rz(-pi) q[1];
rz(-0.97864464) q[2];
sx q[2];
rz(-2.2740764) q[2];
sx q[2];
rz(-0.58766174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0862217) q[1];
sx q[1];
rz(-1.1495591) q[1];
sx q[1];
rz(1.4122541) q[1];
rz(-pi) q[2];
rz(1.1743594) q[3];
sx q[3];
rz(-1.160495) q[3];
sx q[3];
rz(-2.045971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5469024) q[2];
sx q[2];
rz(-1.12744) q[2];
sx q[2];
rz(-0.31000578) q[2];
rz(1.7075432) q[3];
sx q[3];
rz(-1.5177582) q[3];
sx q[3];
rz(-1.8967418) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2731648) q[0];
sx q[0];
rz(-0.36968601) q[0];
sx q[0];
rz(-1.0721068) q[0];
rz(1.3605236) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(-1.2881813) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53258713) q[0];
sx q[0];
rz(-2.3095248) q[0];
sx q[0];
rz(-0.84858507) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6990876) q[2];
sx q[2];
rz(-2.1019843) q[2];
sx q[2];
rz(2.1489904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2350504) q[1];
sx q[1];
rz(-1.8120011) q[1];
sx q[1];
rz(-1.3969087) q[1];
x q[2];
rz(-1.0800519) q[3];
sx q[3];
rz(-0.86378686) q[3];
sx q[3];
rz(0.81719993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0518904) q[2];
sx q[2];
rz(-1.5081729) q[2];
sx q[2];
rz(1.7066329) q[2];
rz(0.47576225) q[3];
sx q[3];
rz(-0.69403726) q[3];
sx q[3];
rz(-0.3835052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39110228) q[0];
sx q[0];
rz(-0.35677156) q[0];
sx q[0];
rz(1.9619036) q[0];
rz(2.7791595) q[1];
sx q[1];
rz(-2.0795627) q[1];
sx q[1];
rz(1.5865883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3192958) q[0];
sx q[0];
rz(-1.8916115) q[0];
sx q[0];
rz(-1.7309624) q[0];
rz(-pi) q[1];
rz(-1.0969504) q[2];
sx q[2];
rz(-2.1948177) q[2];
sx q[2];
rz(2.129675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5240979) q[1];
sx q[1];
rz(-0.84071181) q[1];
sx q[1];
rz(0.9690271) q[1];
rz(0.15896564) q[3];
sx q[3];
rz(-1.0086802) q[3];
sx q[3];
rz(-1.5288439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72593752) q[2];
sx q[2];
rz(-2.026181) q[2];
sx q[2];
rz(1.9647145) q[2];
rz(-0.40500179) q[3];
sx q[3];
rz(-2.160725) q[3];
sx q[3];
rz(2.7485031) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28615752) q[0];
sx q[0];
rz(-0.09849184) q[0];
sx q[0];
rz(-1.3354907) q[0];
rz(-0.91241854) q[1];
sx q[1];
rz(-1.4211979) q[1];
sx q[1];
rz(-0.040249912) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4064317) q[0];
sx q[0];
rz(-1.7393775) q[0];
sx q[0];
rz(-0.038589007) q[0];
rz(-pi) q[1];
rz(2.3715437) q[2];
sx q[2];
rz(-2.382394) q[2];
sx q[2];
rz(-0.46261132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.16791028) q[1];
sx q[1];
rz(-1.0492322) q[1];
sx q[1];
rz(0.48814623) q[1];
rz(1.7829127) q[3];
sx q[3];
rz(-2.3276442) q[3];
sx q[3];
rz(0.98462405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6054298) q[2];
sx q[2];
rz(-1.1417814) q[2];
sx q[2];
rz(0.52967581) q[2];
rz(-1.7236494) q[3];
sx q[3];
rz(-2.6726186) q[3];
sx q[3];
rz(-0.50122195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0025075992) q[0];
sx q[0];
rz(-1.7739828) q[0];
sx q[0];
rz(-3.1097155) q[0];
rz(1.9335951) q[1];
sx q[1];
rz(-0.54787689) q[1];
sx q[1];
rz(-0.30778232) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2805177) q[0];
sx q[0];
rz(-1.3043748) q[0];
sx q[0];
rz(2.0841588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8372907) q[2];
sx q[2];
rz(-1.6937758) q[2];
sx q[2];
rz(-0.83713712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3584328) q[1];
sx q[1];
rz(-0.98020411) q[1];
sx q[1];
rz(-2.5161051) q[1];
rz(-pi) q[2];
rz(-1.4416306) q[3];
sx q[3];
rz(-1.4291479) q[3];
sx q[3];
rz(0.76412898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53177437) q[2];
sx q[2];
rz(-2.4419624) q[2];
sx q[2];
rz(-0.043665234) q[2];
rz(-1.4788491) q[3];
sx q[3];
rz(-0.50373977) q[3];
sx q[3];
rz(-1.3879205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762348) q[0];
sx q[0];
rz(-1.4690514) q[0];
sx q[0];
rz(1.5734191) q[0];
rz(0.028989446) q[1];
sx q[1];
rz(-2.034076) q[1];
sx q[1];
rz(-2.1709002) q[1];
rz(-1.3503648) q[2];
sx q[2];
rz(-1.6904808) q[2];
sx q[2];
rz(-1.3338989) q[2];
rz(-2.3787712) q[3];
sx q[3];
rz(-1.1991432) q[3];
sx q[3];
rz(-1.3292718) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
