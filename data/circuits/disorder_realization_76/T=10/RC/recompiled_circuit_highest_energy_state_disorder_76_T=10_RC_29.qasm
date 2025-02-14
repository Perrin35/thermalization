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
rz(-0.79987502) q[0];
sx q[0];
rz(-0.81698155) q[0];
sx q[0];
rz(-0.45720994) q[0];
rz(-2.7297821) q[1];
sx q[1];
rz(-1.6219985) q[1];
sx q[1];
rz(2.5817459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3362387) q[0];
sx q[0];
rz(-1.8913664) q[0];
sx q[0];
rz(0.47236116) q[0];
rz(-2.855515) q[2];
sx q[2];
rz(-2.0564579) q[2];
sx q[2];
rz(-2.3229675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35750439) q[1];
sx q[1];
rz(-2.5604446) q[1];
sx q[1];
rz(0.52304348) q[1];
x q[2];
rz(2.9212037) q[3];
sx q[3];
rz(-1.7145836) q[3];
sx q[3];
rz(1.4291562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7655699) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(1.8454856) q[2];
rz(-2.5206595) q[3];
sx q[3];
rz(-0.66578484) q[3];
sx q[3];
rz(-0.92500979) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6821482) q[0];
sx q[0];
rz(-0.58732533) q[0];
sx q[0];
rz(-2.4272954) q[0];
rz(2.4192877) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(1.3246271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52377578) q[0];
sx q[0];
rz(-2.0156228) q[0];
sx q[0];
rz(-0.33828783) q[0];
x q[1];
rz(-0.73280735) q[2];
sx q[2];
rz(-1.5298973) q[2];
sx q[2];
rz(-2.3529904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8216422) q[1];
sx q[1];
rz(-0.40377235) q[1];
sx q[1];
rz(-0.41492226) q[1];
x q[2];
rz(0.12217317) q[3];
sx q[3];
rz(-0.88910149) q[3];
sx q[3];
rz(-0.091191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.033919949) q[2];
sx q[2];
rz(-1.4227957) q[2];
sx q[2];
rz(2.1495492) q[2];
rz(-0.77573675) q[3];
sx q[3];
rz(-2.3596767) q[3];
sx q[3];
rz(-1.0824664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.42763212) q[0];
sx q[0];
rz(-1.7949224) q[0];
sx q[0];
rz(-2.2295075) q[0];
rz(-2.7569547) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(0.49547637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3289483) q[0];
sx q[0];
rz(-1.616713) q[0];
sx q[0];
rz(1.4981235) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74929535) q[2];
sx q[2];
rz(-1.4297419) q[2];
sx q[2];
rz(-0.8666477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66049536) q[1];
sx q[1];
rz(-2.4156791) q[1];
sx q[1];
rz(-3.0654415) q[1];
rz(0.94888249) q[3];
sx q[3];
rz(-1.5601741) q[3];
sx q[3];
rz(2.0257906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2428525) q[2];
sx q[2];
rz(-2.0169368) q[2];
sx q[2];
rz(-0.43453547) q[2];
rz(-3.0430326) q[3];
sx q[3];
rz(-1.4222654) q[3];
sx q[3];
rz(-2.3677473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74715215) q[0];
sx q[0];
rz(-0.23500615) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(-2.0647743) q[1];
sx q[1];
rz(-1.4067168) q[1];
sx q[1];
rz(-1.9487618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6859602) q[0];
sx q[0];
rz(-1.3979027) q[0];
sx q[0];
rz(-2.8894823) q[0];
rz(-2.3770726) q[2];
sx q[2];
rz(-2.9540375) q[2];
sx q[2];
rz(-2.9775567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4095278) q[1];
sx q[1];
rz(-1.3573779) q[1];
sx q[1];
rz(-1.7312538) q[1];
rz(-pi) q[2];
rz(1.3681196) q[3];
sx q[3];
rz(-2.2619288) q[3];
sx q[3];
rz(2.1512669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8863135) q[2];
sx q[2];
rz(-2.2463319) q[2];
sx q[2];
rz(-2.8423584) q[2];
rz(-2.8507161) q[3];
sx q[3];
rz(-1.8894922) q[3];
sx q[3];
rz(-2.4860184) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84151477) q[0];
sx q[0];
rz(-0.97989196) q[0];
sx q[0];
rz(-0.078068659) q[0];
rz(-0.95371753) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(1.5740707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2497319) q[0];
sx q[0];
rz(-1.390269) q[0];
sx q[0];
rz(2.2523227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7360953) q[2];
sx q[2];
rz(-0.73842305) q[2];
sx q[2];
rz(-0.15106311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6202653) q[1];
sx q[1];
rz(-0.6129188) q[1];
sx q[1];
rz(1.5144203) q[1];
x q[2];
rz(-3.0573274) q[3];
sx q[3];
rz(-2.4467496) q[3];
sx q[3];
rz(0.39277276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2786402) q[2];
sx q[2];
rz(-2.7095257) q[2];
sx q[2];
rz(0.26818177) q[2];
rz(2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(-1.4485654) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0443403) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(0.46446717) q[0];
rz(0.52344549) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(-2.4109667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7775041) q[0];
sx q[0];
rz(-0.83151649) q[0];
sx q[0];
rz(2.9502908) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36412698) q[2];
sx q[2];
rz(-2.4362323) q[2];
sx q[2];
rz(0.10157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9537415) q[1];
sx q[1];
rz(-1.9914728) q[1];
sx q[1];
rz(-0.71297808) q[1];
x q[2];
rz(-1.7028125) q[3];
sx q[3];
rz(-1.321133) q[3];
sx q[3];
rz(-2.3625847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0940493) q[2];
sx q[2];
rz(-1.0796248) q[2];
sx q[2];
rz(-3.003982) q[2];
rz(2.008647) q[3];
sx q[3];
rz(-1.9652941) q[3];
sx q[3];
rz(1.8784116) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9889744) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(0.033893943) q[0];
rz(2.7821817) q[1];
sx q[1];
rz(-1.6172599) q[1];
sx q[1];
rz(0.83438897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9890404) q[0];
sx q[0];
rz(-1.1942399) q[0];
sx q[0];
rz(1.9045619) q[0];
x q[1];
rz(-1.0831725) q[2];
sx q[2];
rz(-0.94218035) q[2];
sx q[2];
rz(-2.4982029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8270478) q[1];
sx q[1];
rz(-2.8706708) q[1];
sx q[1];
rz(-1.9344058) q[1];
rz(-pi) q[2];
rz(-1.5079751) q[3];
sx q[3];
rz(-0.90017002) q[3];
sx q[3];
rz(2.2986029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.413201) q[2];
sx q[2];
rz(-0.30343702) q[2];
sx q[2];
rz(2.0835853) q[2];
rz(-0.47438619) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(2.0856196) q[3];
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
rz(-pi) q[3];
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
rz(2.0673085) q[0];
sx q[0];
rz(-1.625076) q[0];
sx q[0];
rz(1.7838595) q[0];
rz(0.43074295) q[1];
sx q[1];
rz(-1.6669225) q[1];
sx q[1];
rz(0.46708435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3586836) q[0];
sx q[0];
rz(-3.0984146) q[0];
sx q[0];
rz(-1.6048649) q[0];
rz(-pi) q[1];
rz(1.378304) q[2];
sx q[2];
rz(-2.2644832) q[2];
sx q[2];
rz(-0.81499962) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.704485) q[1];
sx q[1];
rz(-0.52981716) q[1];
sx q[1];
rz(3.1049314) q[1];
rz(-pi) q[2];
rz(-0.93293087) q[3];
sx q[3];
rz(-1.6137505) q[3];
sx q[3];
rz(1.651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0507386) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(1.0806855) q[2];
rz(2.4596227) q[3];
sx q[3];
rz(-2.3365648) q[3];
sx q[3];
rz(-2.4431156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.4526116) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(2.3338351) q[0];
rz(-1.0001596) q[1];
sx q[1];
rz(-2.2554485) q[1];
sx q[1];
rz(-2.3449786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611872) q[0];
sx q[0];
rz(-3.0356501) q[0];
sx q[0];
rz(1.1619912) q[0];
rz(-pi) q[1];
rz(0.32795017) q[2];
sx q[2];
rz(-2.4900576) q[2];
sx q[2];
rz(-2.8095989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0137009) q[1];
sx q[1];
rz(-1.7096448) q[1];
sx q[1];
rz(1.9816887) q[1];
rz(-2.2439205) q[3];
sx q[3];
rz(-0.39413777) q[3];
sx q[3];
rz(1.2887312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5920068) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(0.31527147) q[2];
rz(1.5677876) q[3];
sx q[3];
rz(-1.9939634) q[3];
sx q[3];
rz(-2.018759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386117) q[0];
sx q[0];
rz(-2.4760315) q[0];
sx q[0];
rz(-0.82157201) q[0];
rz(2.6577677) q[1];
sx q[1];
rz(-1.4204104) q[1];
sx q[1];
rz(-0.22463591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53157887) q[0];
sx q[0];
rz(-1.2481127) q[0];
sx q[0];
rz(-2.9831397) q[0];
x q[1];
rz(2.3326229) q[2];
sx q[2];
rz(-2.3270049) q[2];
sx q[2];
rz(-0.23216173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37525081) q[1];
sx q[1];
rz(-2.2415889) q[1];
sx q[1];
rz(2.8743582) q[1];
rz(-pi) q[2];
rz(-1.3566689) q[3];
sx q[3];
rz(-1.4487293) q[3];
sx q[3];
rz(1.3510974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75767526) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(0.095001027) q[2];
rz(1.4875686) q[3];
sx q[3];
rz(-2.5520958) q[3];
sx q[3];
rz(2.3768363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8944396) q[0];
sx q[0];
rz(-0.40774397) q[0];
sx q[0];
rz(2.8751873) q[0];
rz(-15/(8*pi)) q[1];
sx q[1];
rz(-1.4529556) q[1];
sx q[1];
rz(-1.6642889) q[1];
rz(2.2483027) q[2];
sx q[2];
rz(-1.8532475) q[2];
sx q[2];
rz(-0.38459323) q[2];
rz(-0.0048051759) q[3];
sx q[3];
rz(-2.7660696) q[3];
sx q[3];
rz(0.39733359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
