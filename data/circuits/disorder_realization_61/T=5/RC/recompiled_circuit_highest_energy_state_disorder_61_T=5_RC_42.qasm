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
rz(-2.3075624) q[0];
sx q[0];
rz(-0.40606719) q[0];
sx q[0];
rz(-1.9598444) q[0];
rz(-1.4015247) q[1];
sx q[1];
rz(-2.6670246) q[1];
sx q[1];
rz(-0.13303703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8048153) q[0];
sx q[0];
rz(-2.204737) q[0];
sx q[0];
rz(-0.32685771) q[0];
rz(-pi) q[1];
rz(2.4394717) q[2];
sx q[2];
rz(-0.63533855) q[2];
sx q[2];
rz(-1.8566657) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.592748) q[1];
sx q[1];
rz(-1.2557552) q[1];
sx q[1];
rz(2.9058019) q[1];
x q[2];
rz(2.1429012) q[3];
sx q[3];
rz(-1.8810728) q[3];
sx q[3];
rz(-2.7039122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0066321) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(0.68438619) q[2];
rz(2.0002174) q[3];
sx q[3];
rz(-2.8751774) q[3];
sx q[3];
rz(0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.77966493) q[0];
sx q[0];
rz(-2.4517038) q[0];
sx q[0];
rz(-0.46098125) q[0];
rz(-2.0224679) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(0.31203312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0587766) q[0];
sx q[0];
rz(-3.1110373) q[0];
sx q[0];
rz(-2.729134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.076032555) q[2];
sx q[2];
rz(-2.3302148) q[2];
sx q[2];
rz(1.8234314) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2300517) q[1];
sx q[1];
rz(-1.576606) q[1];
sx q[1];
rz(1.7464386) q[1];
rz(-pi) q[2];
rz(0.9189689) q[3];
sx q[3];
rz(-1.7953835) q[3];
sx q[3];
rz(2.74204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.95388874) q[2];
sx q[2];
rz(-2.0362594) q[2];
sx q[2];
rz(0.40404955) q[2];
rz(-2.7814501) q[3];
sx q[3];
rz(-1.6481684) q[3];
sx q[3];
rz(2.4586316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78524154) q[0];
sx q[0];
rz(-2.1389565) q[0];
sx q[0];
rz(2.8175765) q[0];
rz(-2.0815381) q[1];
sx q[1];
rz(-0.29393229) q[1];
sx q[1];
rz(2.4411328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4361079) q[0];
sx q[0];
rz(-1.4027081) q[0];
sx q[0];
rz(-1.7343434) q[0];
x q[1];
rz(2.1164139) q[2];
sx q[2];
rz(-2.7433925) q[2];
sx q[2];
rz(0.37665483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6571275) q[1];
sx q[1];
rz(-1.336369) q[1];
sx q[1];
rz(2.3695494) q[1];
rz(-pi) q[2];
rz(1.9869119) q[3];
sx q[3];
rz(-2.6857102) q[3];
sx q[3];
rz(-1.0030163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5141653) q[2];
sx q[2];
rz(-0.24719396) q[2];
sx q[2];
rz(2.3449281) q[2];
rz(1.944815) q[3];
sx q[3];
rz(-1.6007042) q[3];
sx q[3];
rz(2.5483907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6881123) q[0];
sx q[0];
rz(-1.0013094) q[0];
sx q[0];
rz(1.3230811) q[0];
rz(2.6755013) q[1];
sx q[1];
rz(-2.2345462) q[1];
sx q[1];
rz(-1.1493433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1150165) q[0];
sx q[0];
rz(-1.1372024) q[0];
sx q[0];
rz(1.0267016) q[0];
rz(-pi) q[1];
rz(1.0428814) q[2];
sx q[2];
rz(-2.4529874) q[2];
sx q[2];
rz(1.3618038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.016952) q[1];
sx q[1];
rz(-1.5257383) q[1];
sx q[1];
rz(-2.9393815) q[1];
x q[2];
rz(-1.9051311) q[3];
sx q[3];
rz(-1.7133681) q[3];
sx q[3];
rz(-2.4041374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.484802) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(0.82526866) q[2];
rz(1.8615865) q[3];
sx q[3];
rz(-1.620159) q[3];
sx q[3];
rz(-0.72117225) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.698302) q[0];
sx q[0];
rz(-1.8809141) q[0];
sx q[0];
rz(-1.3662421) q[0];
rz(-1.0900991) q[1];
sx q[1];
rz(-1.5049728) q[1];
sx q[1];
rz(1.3390138) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3109717) q[0];
sx q[0];
rz(-1.3960103) q[0];
sx q[0];
rz(2.4808148) q[0];
rz(-0.62241222) q[2];
sx q[2];
rz(-0.96660766) q[2];
sx q[2];
rz(-0.46350086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3303393) q[1];
sx q[1];
rz(-0.24345466) q[1];
sx q[1];
rz(-0.50916839) q[1];
rz(-pi) q[2];
rz(1.9129778) q[3];
sx q[3];
rz(-1.2552731) q[3];
sx q[3];
rz(0.70759892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1947386) q[2];
sx q[2];
rz(-1.4118492) q[2];
sx q[2];
rz(2.8080688) q[2];
rz(2.5528095) q[3];
sx q[3];
rz(-0.47313658) q[3];
sx q[3];
rz(-2.3401882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(2.3280585) q[0];
sx q[0];
rz(-1.3575587) q[0];
sx q[0];
rz(-1.300977) q[0];
rz(1.0282372) q[1];
sx q[1];
rz(-0.49003777) q[1];
sx q[1];
rz(-0.44233826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72318769) q[0];
sx q[0];
rz(-0.53779477) q[0];
sx q[0];
rz(-2.6502953) q[0];
rz(-pi) q[1];
rz(2.6643941) q[2];
sx q[2];
rz(-2.0522567) q[2];
sx q[2];
rz(1.9419326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42980151) q[1];
sx q[1];
rz(-0.45908976) q[1];
sx q[1];
rz(-2.7613822) q[1];
rz(-pi) q[2];
rz(2.836727) q[3];
sx q[3];
rz(-2.8072522) q[3];
sx q[3];
rz(-1.6337194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77171317) q[2];
sx q[2];
rz(-1.5934725) q[2];
sx q[2];
rz(-0.4717353) q[2];
rz(-1.6996023) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(2.877318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(3.0278552) q[0];
sx q[0];
rz(-1.4677445) q[0];
sx q[0];
rz(2.9748528) q[0];
rz(2.3475504) q[1];
sx q[1];
rz(-1.9414976) q[1];
sx q[1];
rz(-0.099743191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31064311) q[0];
sx q[0];
rz(-1.8825873) q[0];
sx q[0];
rz(-0.98592178) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95364596) q[2];
sx q[2];
rz(-1.8815478) q[2];
sx q[2];
rz(1.0342395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.21799) q[1];
sx q[1];
rz(-2.1384708) q[1];
sx q[1];
rz(1.7204549) q[1];
x q[2];
rz(2.7230303) q[3];
sx q[3];
rz(-2.1893614) q[3];
sx q[3];
rz(1.5582066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2793067) q[2];
sx q[2];
rz(-2.0760832) q[2];
sx q[2];
rz(3.0118937) q[2];
rz(1.8779523) q[3];
sx q[3];
rz(-1.940515) q[3];
sx q[3];
rz(2.9960846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.170914) q[0];
sx q[0];
rz(-2.2190974) q[0];
sx q[0];
rz(-0.42651919) q[0];
rz(-2.049394) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(-2.138864) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9364081) q[0];
sx q[0];
rz(-0.6615839) q[0];
sx q[0];
rz(-1.1422045) q[0];
rz(-0.72288798) q[2];
sx q[2];
rz(-1.6405511) q[2];
sx q[2];
rz(2.7107014) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0074689) q[1];
sx q[1];
rz(-1.136131) q[1];
sx q[1];
rz(2.77209) q[1];
x q[2];
rz(1.1158607) q[3];
sx q[3];
rz(-1.3580033) q[3];
sx q[3];
rz(1.4404595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87454522) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(-0.12346539) q[2];
rz(0.63001436) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(-0.71648487) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2444839) q[0];
sx q[0];
rz(-2.4429595) q[0];
sx q[0];
rz(0.83936349) q[0];
rz(0.15946236) q[1];
sx q[1];
rz(-2.1493252) q[1];
sx q[1];
rz(-2.2119567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2004505) q[0];
sx q[0];
rz(-1.8243891) q[0];
sx q[0];
rz(-0.7483878) q[0];
rz(-0.59105525) q[2];
sx q[2];
rz(-1.4690217) q[2];
sx q[2];
rz(-2.7391694) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1222555) q[1];
sx q[1];
rz(-2.4474065) q[1];
sx q[1];
rz(-1.3744785) q[1];
x q[2];
rz(-0.22945424) q[3];
sx q[3];
rz(-0.94564518) q[3];
sx q[3];
rz(-1.8895686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1975501) q[2];
sx q[2];
rz(-1.1669179) q[2];
sx q[2];
rz(-0.76773947) q[2];
rz(2.7770216) q[3];
sx q[3];
rz(-1.7578099) q[3];
sx q[3];
rz(-3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0070873) q[0];
sx q[0];
rz(-0.61408478) q[0];
sx q[0];
rz(-2.2579204) q[0];
rz(-0.20835749) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(1.8066033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7430089) q[0];
sx q[0];
rz(-1.5129876) q[0];
sx q[0];
rz(3.0763118) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7572549) q[2];
sx q[2];
rz(-2.0907195) q[2];
sx q[2];
rz(0.35042074) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9197651) q[1];
sx q[1];
rz(-2.1073256) q[1];
sx q[1];
rz(0.51862688) q[1];
x q[2];
rz(-2.317071) q[3];
sx q[3];
rz(-1.0484277) q[3];
sx q[3];
rz(2.8461412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9779382) q[2];
sx q[2];
rz(-1.5074707) q[2];
sx q[2];
rz(-3.0955691) q[2];
rz(-0.16511354) q[3];
sx q[3];
rz(-0.29296994) q[3];
sx q[3];
rz(1.2142115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56275) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(0.1651172) q[1];
sx q[1];
rz(-1.6830291) q[1];
sx q[1];
rz(2.8370636) q[1];
rz(-2.6982735) q[2];
sx q[2];
rz(-1.2740612) q[2];
sx q[2];
rz(-3.086003) q[2];
rz(-2.6611609) q[3];
sx q[3];
rz(-2.0365002) q[3];
sx q[3];
rz(-1.0040803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
