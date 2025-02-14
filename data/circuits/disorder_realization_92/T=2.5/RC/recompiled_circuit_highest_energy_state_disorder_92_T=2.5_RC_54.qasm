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
rz(0.25843698) q[0];
sx q[0];
rz(-0.28764495) q[0];
sx q[0];
rz(-1.1076466) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(-2.5133361) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4391651) q[0];
sx q[0];
rz(-1.747923) q[0];
sx q[0];
rz(3.0217373) q[0];
rz(-pi) q[1];
rz(2.4187367) q[2];
sx q[2];
rz(-1.8713163) q[2];
sx q[2];
rz(-2.3006257) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3648277) q[1];
sx q[1];
rz(-1.3421632) q[1];
sx q[1];
rz(1.6079812) q[1];
x q[2];
rz(-0.89723779) q[3];
sx q[3];
rz(-0.64770401) q[3];
sx q[3];
rz(3.0239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85311741) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(-0.34118578) q[2];
rz(-0.8902542) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(-0.53513479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5219236) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(2.4271915) q[0];
rz(1.5419386) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(-0.094706789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9908087) q[0];
sx q[0];
rz(-1.8903194) q[0];
sx q[0];
rz(3.0421769) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0422499) q[2];
sx q[2];
rz(-0.65644342) q[2];
sx q[2];
rz(3.0489056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4112139) q[1];
sx q[1];
rz(-1.3347365) q[1];
sx q[1];
rz(-0.35309955) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50928874) q[3];
sx q[3];
rz(-1.4319766) q[3];
sx q[3];
rz(-1.2535439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9793205) q[2];
sx q[2];
rz(-1.8956192) q[2];
sx q[2];
rz(-0.47404131) q[2];
rz(-2.3794409) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(-2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-0.85292029) q[0];
rz(0.22311738) q[1];
sx q[1];
rz(-0.47839636) q[1];
sx q[1];
rz(0.51582897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0313161) q[0];
sx q[0];
rz(-3.0678684) q[0];
sx q[0];
rz(-0.81839968) q[0];
rz(0.25571172) q[2];
sx q[2];
rz(-0.94921321) q[2];
sx q[2];
rz(-2.5598124) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3577376) q[1];
sx q[1];
rz(-0.92354316) q[1];
sx q[1];
rz(-2.3308862) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9733377) q[3];
sx q[3];
rz(-0.73930891) q[3];
sx q[3];
rz(0.69956796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3333266) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(-0.094956368) q[2];
rz(2.700108) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(-1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74986356) q[0];
sx q[0];
rz(-2.5666105) q[0];
sx q[0];
rz(-1.0561426) q[0];
rz(-2.1978343) q[1];
sx q[1];
rz(-2.8902003) q[1];
sx q[1];
rz(1.4518849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34925845) q[0];
sx q[0];
rz(-0.99674388) q[0];
sx q[0];
rz(0.7416772) q[0];
x q[1];
rz(1.1805492) q[2];
sx q[2];
rz(-1.1392635) q[2];
sx q[2];
rz(0.081588946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31409594) q[1];
sx q[1];
rz(-0.24330595) q[1];
sx q[1];
rz(-2.9929586) q[1];
x q[2];
rz(0.11066505) q[3];
sx q[3];
rz(-1.5783211) q[3];
sx q[3];
rz(-1.8845817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(-2.441414) q[2];
rz(-0.87107301) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(-1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-2.6569271) q[0];
rz(2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(-2.7427618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8190552) q[0];
sx q[0];
rz(-0.56542552) q[0];
sx q[0];
rz(2.8877693) q[0];
rz(-pi) q[1];
rz(2.4381934) q[2];
sx q[2];
rz(-1.8336772) q[2];
sx q[2];
rz(0.92353283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5113885) q[1];
sx q[1];
rz(-0.92899073) q[1];
sx q[1];
rz(0.46496572) q[1];
rz(0.9449686) q[3];
sx q[3];
rz(-2.3381013) q[3];
sx q[3];
rz(3.0253493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47348076) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(-2.2637746) q[2];
rz(0.38464883) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28410742) q[0];
sx q[0];
rz(-2.6326038) q[0];
sx q[0];
rz(-2.8711163) q[0];
rz(-0.30666223) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(0.57672966) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33615112) q[0];
sx q[0];
rz(-2.9362943) q[0];
sx q[0];
rz(2.0144256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2606196) q[2];
sx q[2];
rz(-1.8554129) q[2];
sx q[2];
rz(1.5280444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2564078) q[1];
sx q[1];
rz(-1.9659623) q[1];
sx q[1];
rz(2.1705519) q[1];
rz(-pi) q[2];
rz(-0.33476013) q[3];
sx q[3];
rz(-1.5575512) q[3];
sx q[3];
rz(-0.56559004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(-0.97037399) q[2];
rz(-2.7905285) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2358667) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(0.50091499) q[0];
rz(-0.72055912) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(-0.93773425) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36762992) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(-3.0995661) q[0];
rz(-pi) q[1];
rz(0.65848041) q[2];
sx q[2];
rz(-2.5665433) q[2];
sx q[2];
rz(-2.5031896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8775505) q[1];
sx q[1];
rz(-2.1960947) q[1];
sx q[1];
rz(-0.45035119) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1295007) q[3];
sx q[3];
rz(-2.4547057) q[3];
sx q[3];
rz(-1.505132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0096036) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(0.84849882) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(-1.2257303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41035143) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(2.248726) q[0];
rz(0.061575312) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(0.16199131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6380986) q[0];
sx q[0];
rz(-1.5862521) q[0];
sx q[0];
rz(-1.4686122) q[0];
rz(3.0062469) q[2];
sx q[2];
rz(-0.70643808) q[2];
sx q[2];
rz(1.5069697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87016856) q[1];
sx q[1];
rz(-1.4335911) q[1];
sx q[1];
rz(-2.8297564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4668767) q[3];
sx q[3];
rz(-1.6641518) q[3];
sx q[3];
rz(0.68478497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.091247678) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(2.4166935) q[2];
rz(-0.27724087) q[3];
sx q[3];
rz(-0.96479779) q[3];
sx q[3];
rz(0.086517081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25170931) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(-3.1157893) q[0];
rz(1.39894) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(-0.10765156) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336595) q[0];
sx q[0];
rz(-1.7042394) q[0];
sx q[0];
rz(0.99876499) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3978762) q[2];
sx q[2];
rz(-0.60600835) q[2];
sx q[2];
rz(1.5569937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6580769) q[1];
sx q[1];
rz(-0.47464322) q[1];
sx q[1];
rz(-0.078981699) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4156451) q[3];
sx q[3];
rz(-0.11501139) q[3];
sx q[3];
rz(1.0105159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4699576) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(-2.5661902) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53330082) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(-2.4285512) q[0];
rz(-2.9167922) q[1];
sx q[1];
rz(-2.5495922) q[1];
sx q[1];
rz(-1.0932659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521234) q[0];
sx q[0];
rz(-0.88406815) q[0];
sx q[0];
rz(2.5144469) q[0];
rz(-pi) q[1];
rz(-2.9399038) q[2];
sx q[2];
rz(-1.8400157) q[2];
sx q[2];
rz(1.9112916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8902407) q[1];
sx q[1];
rz(-1.6820434) q[1];
sx q[1];
rz(2.6835068) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9841927) q[3];
sx q[3];
rz(-1.6781312) q[3];
sx q[3];
rz(1.8183551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1405868) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(2.801753) q[2];
rz(0.042424399) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(-2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18702678) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(-2.3595702) q[1];
sx q[1];
rz(-2.2047058) q[1];
sx q[1];
rz(2.4155736) q[1];
rz(2.4154631) q[2];
sx q[2];
rz(-0.48116044) q[2];
sx q[2];
rz(3.0143123) q[2];
rz(1.9966765) q[3];
sx q[3];
rz(-1.3419587) q[3];
sx q[3];
rz(1.6423196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
