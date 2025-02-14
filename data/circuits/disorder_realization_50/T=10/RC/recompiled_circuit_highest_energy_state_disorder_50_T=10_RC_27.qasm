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
rz(-2.6900238) q[0];
sx q[0];
rz(-1.4181925) q[0];
sx q[0];
rz(-0.42887846) q[0];
rz(2.9984527) q[1];
sx q[1];
rz(-1.0909456) q[1];
sx q[1];
rz(-0.080088869) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6038136) q[0];
sx q[0];
rz(-1.5575557) q[0];
sx q[0];
rz(-0.079188345) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2099138) q[2];
sx q[2];
rz(-1.4744802) q[2];
sx q[2];
rz(2.148984) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0804008) q[1];
sx q[1];
rz(-0.75729232) q[1];
sx q[1];
rz(1.1822834) q[1];
rz(-pi) q[2];
rz(-1.4065885) q[3];
sx q[3];
rz(-2.1421711) q[3];
sx q[3];
rz(0.77916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0668138) q[2];
sx q[2];
rz(-1.0513693) q[2];
sx q[2];
rz(-1.7195513) q[2];
rz(-1.7285796) q[3];
sx q[3];
rz(-1.541411) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-0.3438147) q[1];
sx q[1];
rz(-0.75210345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.89276) q[0];
sx q[0];
rz(-1.1247096) q[0];
sx q[0];
rz(-0.89527881) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0235748) q[2];
sx q[2];
rz(-0.22023295) q[2];
sx q[2];
rz(-1.756246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57422728) q[1];
sx q[1];
rz(-0.53814473) q[1];
sx q[1];
rz(1.7848915) q[1];
x q[2];
rz(-1.7957572) q[3];
sx q[3];
rz(-1.3153425) q[3];
sx q[3];
rz(-2.6663702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.077945262) q[2];
sx q[2];
rz(-2.2293978) q[2];
sx q[2];
rz(-2.4845541) q[2];
rz(-2.4296203) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(-2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7065358) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(1.0414498) q[0];
rz(-2.7667747) q[1];
sx q[1];
rz(-1.5033787) q[1];
sx q[1];
rz(0.55694881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9948363) q[0];
sx q[0];
rz(-2.2731011) q[0];
sx q[0];
rz(2.0392787) q[0];
rz(-pi) q[1];
rz(-1.1191423) q[2];
sx q[2];
rz(-1.5030393) q[2];
sx q[2];
rz(0.51167831) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8136602) q[1];
sx q[1];
rz(-0.90444649) q[1];
sx q[1];
rz(1.7009237) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30681614) q[3];
sx q[3];
rz(-2.4673438) q[3];
sx q[3];
rz(2.0381272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(2.7064145) q[2];
rz(2.6168881) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(3.0068126) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10114577) q[0];
sx q[0];
rz(-1.0643767) q[0];
sx q[0];
rz(-1.5144298) q[0];
rz(-1.0487522) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.705816) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.993282) q[0];
sx q[0];
rz(-1.7459501) q[0];
sx q[0];
rz(-1.5272642) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9210774) q[2];
sx q[2];
rz(-1.0547914) q[2];
sx q[2];
rz(-0.28921858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8183349) q[1];
sx q[1];
rz(-1.0498974) q[1];
sx q[1];
rz(3.0677615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9463598) q[3];
sx q[3];
rz(-0.55259669) q[3];
sx q[3];
rz(-2.432807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40020308) q[2];
sx q[2];
rz(-1.783739) q[2];
sx q[2];
rz(-0.025156585) q[2];
rz(-1.6545506) q[3];
sx q[3];
rz(-0.39792037) q[3];
sx q[3];
rz(-2.3253843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0461034) q[0];
sx q[0];
rz(-2.4920721) q[0];
sx q[0];
rz(-0.15897861) q[0];
rz(-2.036463) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(2.9172858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8858356) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(1.095039) q[0];
rz(-1.6231027) q[2];
sx q[2];
rz(-2.5807267) q[2];
sx q[2];
rz(0.069151783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.060886) q[1];
sx q[1];
rz(-1.4554441) q[1];
sx q[1];
rz(-1.0447698) q[1];
rz(-pi) q[2];
rz(-2.5644825) q[3];
sx q[3];
rz(-1.6770067) q[3];
sx q[3];
rz(-2.4662227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.20532456) q[2];
sx q[2];
rz(-1.3621796) q[2];
sx q[2];
rz(0.74860191) q[2];
rz(0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.1614138) q[0];
sx q[0];
rz(-0.36853376) q[0];
sx q[0];
rz(-0.73032105) q[0];
rz(-1.5018564) q[1];
sx q[1];
rz(-1.7514936) q[1];
sx q[1];
rz(-2.9072442) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4627107) q[0];
sx q[0];
rz(-2.2214799) q[0];
sx q[0];
rz(2.1367055) q[0];
rz(-0.5989093) q[2];
sx q[2];
rz(-1.0157758) q[2];
sx q[2];
rz(-2.2428494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6987353) q[1];
sx q[1];
rz(-2.9720829) q[1];
sx q[1];
rz(0.28556602) q[1];
rz(-pi) q[2];
rz(-1.430351) q[3];
sx q[3];
rz(-1.5101119) q[3];
sx q[3];
rz(-1.3827326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.948287) q[2];
sx q[2];
rz(-2.327658) q[2];
sx q[2];
rz(-1.8602271) q[2];
rz(0.60503259) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(1.9090778) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(0.61304098) q[0];
rz(2.4609861) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.8162762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26215689) q[0];
sx q[0];
rz(-2.5234294) q[0];
sx q[0];
rz(-2.6690041) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6540131) q[2];
sx q[2];
rz(-1.3583379) q[2];
sx q[2];
rz(-2.0709289) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.054106709) q[1];
sx q[1];
rz(-1.0102059) q[1];
sx q[1];
rz(-1.7203548) q[1];
rz(1.9716827) q[3];
sx q[3];
rz(-1.0260149) q[3];
sx q[3];
rz(-2.0044897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94720542) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(0.80805937) q[2];
rz(-0.64588532) q[3];
sx q[3];
rz(-2.8736727) q[3];
sx q[3];
rz(-0.3033692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5658257) q[0];
sx q[0];
rz(-0.72661) q[0];
sx q[0];
rz(2.6988244) q[0];
rz(-2.5495461) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(-2.9659042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4018049) q[0];
sx q[0];
rz(-1.2953838) q[0];
sx q[0];
rz(2.7729183) q[0];
rz(-pi) q[1];
rz(1.6106583) q[2];
sx q[2];
rz(-2.022036) q[2];
sx q[2];
rz(0.892837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1877039) q[1];
sx q[1];
rz(-0.64262455) q[1];
sx q[1];
rz(3.0695555) q[1];
rz(-pi) q[2];
rz(0.91276987) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(2.4778544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(2.234484) q[2];
rz(0.90726566) q[3];
sx q[3];
rz(-1.3943358) q[3];
sx q[3];
rz(-0.11769122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9573145) q[0];
sx q[0];
rz(-0.90183455) q[0];
sx q[0];
rz(-2.9515475) q[0];
rz(1.5826781) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(1.3072026) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35018539) q[0];
sx q[0];
rz(-2.7481347) q[0];
sx q[0];
rz(1.1883452) q[0];
x q[1];
rz(-1.9668504) q[2];
sx q[2];
rz(-2.4017334) q[2];
sx q[2];
rz(-2.6137874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70997256) q[1];
sx q[1];
rz(-0.11259544) q[1];
sx q[1];
rz(-1.8667875) q[1];
rz(1.6862737) q[3];
sx q[3];
rz(-1.1388123) q[3];
sx q[3];
rz(0.40528709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.068453161) q[2];
sx q[2];
rz(-2.0544724) q[2];
sx q[2];
rz(-2.9404409) q[2];
rz(-2.1189832) q[3];
sx q[3];
rz(-2.4590731) q[3];
sx q[3];
rz(-1.6020017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1744743) q[0];
sx q[0];
rz(-2.0956464) q[0];
sx q[0];
rz(2.2651267) q[0];
rz(-2.5190952) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(-0.34271398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4349608) q[0];
sx q[0];
rz(-1.6366658) q[0];
sx q[0];
rz(0.94191282) q[0];
rz(0.12212716) q[2];
sx q[2];
rz(-1.2363648) q[2];
sx q[2];
rz(1.328652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6500998) q[1];
sx q[1];
rz(-2.0180185) q[1];
sx q[1];
rz(-0.69112055) q[1];
rz(-2.7539135) q[3];
sx q[3];
rz(-1.8053384) q[3];
sx q[3];
rz(-2.925433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9524625) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-2.6978037) q[2];
rz(2.022838) q[3];
sx q[3];
rz(-1.5464562) q[3];
sx q[3];
rz(-0.55383033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2832058) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(-3.0781147) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-0.061894682) q[2];
sx q[2];
rz(-1.2409004) q[2];
sx q[2];
rz(0.093766669) q[2];
rz(-1.1984428) q[3];
sx q[3];
rz(-0.59881864) q[3];
sx q[3];
rz(-1.3074085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
