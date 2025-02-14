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
rz(-1.0919548) q[0];
sx q[0];
rz(5.2009861) q[0];
sx q[0];
rz(7.4823147) q[0];
rz(0.24425976) q[1];
sx q[1];
rz(-1.3173988) q[1];
sx q[1];
rz(2.2720845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6586745) q[0];
sx q[0];
rz(-2.344083) q[0];
sx q[0];
rz(1.2748383) q[0];
rz(-pi) q[1];
rz(2.2734322) q[2];
sx q[2];
rz(-2.5294315) q[2];
sx q[2];
rz(-0.56589076) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7092822) q[1];
sx q[1];
rz(-0.73768931) q[1];
sx q[1];
rz(-2.7861315) q[1];
x q[2];
rz(2.364089) q[3];
sx q[3];
rz(-0.99626675) q[3];
sx q[3];
rz(-0.053701775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96183744) q[2];
sx q[2];
rz(-1.9840252) q[2];
sx q[2];
rz(-1.8400787) q[2];
rz(0.88422042) q[3];
sx q[3];
rz(-2.0562462) q[3];
sx q[3];
rz(-1.4475383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41534153) q[0];
sx q[0];
rz(-1.911442) q[0];
sx q[0];
rz(-3.0905261) q[0];
rz(2.6452433) q[1];
sx q[1];
rz(-1.2800848) q[1];
sx q[1];
rz(0.11402093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43085264) q[0];
sx q[0];
rz(-0.56409696) q[0];
sx q[0];
rz(-2.7377955) q[0];
rz(-pi) q[1];
rz(-1.586726) q[2];
sx q[2];
rz(-0.6932879) q[2];
sx q[2];
rz(-0.42809286) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84987243) q[1];
sx q[1];
rz(-2.0270621) q[1];
sx q[1];
rz(0.48584818) q[1];
rz(-2.5258371) q[3];
sx q[3];
rz(-1.0395219) q[3];
sx q[3];
rz(-0.61946509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8337635) q[2];
sx q[2];
rz(-1.9139863) q[2];
sx q[2];
rz(-2.7755136) q[2];
rz(-0.45575538) q[3];
sx q[3];
rz(-0.78211203) q[3];
sx q[3];
rz(-2.8428049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5251821) q[0];
sx q[0];
rz(-0.95935482) q[0];
sx q[0];
rz(-0.33732238) q[0];
rz(1.6711383) q[1];
sx q[1];
rz(-2.2098139) q[1];
sx q[1];
rz(0.09387389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3009143) q[0];
sx q[0];
rz(-0.72193679) q[0];
sx q[0];
rz(-0.11322279) q[0];
x q[1];
rz(-0.63957165) q[2];
sx q[2];
rz(-0.54023114) q[2];
sx q[2];
rz(-0.28499261) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9843238) q[1];
sx q[1];
rz(-0.4254057) q[1];
sx q[1];
rz(-1.5509863) q[1];
rz(0.29993124) q[3];
sx q[3];
rz(-2.3818204) q[3];
sx q[3];
rz(-0.60985669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9570534) q[2];
sx q[2];
rz(-0.55593714) q[2];
sx q[2];
rz(2.8957193) q[2];
rz(-0.18694123) q[3];
sx q[3];
rz(-1.4176466) q[3];
sx q[3];
rz(0.42417446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71640715) q[0];
sx q[0];
rz(-2.462429) q[0];
sx q[0];
rz(0.19234046) q[0];
rz(1.964407) q[1];
sx q[1];
rz(-2.3277551) q[1];
sx q[1];
rz(-0.40843931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41685661) q[0];
sx q[0];
rz(-1.724986) q[0];
sx q[0];
rz(1.5132928) q[0];
x q[1];
rz(-2.2871509) q[2];
sx q[2];
rz(-2.7506094) q[2];
sx q[2];
rz(1.4650844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7725984) q[1];
sx q[1];
rz(-0.97681681) q[1];
sx q[1];
rz(1.0812378) q[1];
x q[2];
rz(2.2564234) q[3];
sx q[3];
rz(-1.7333631) q[3];
sx q[3];
rz(-3.0522961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.283215) q[2];
sx q[2];
rz(-2.242531) q[2];
sx q[2];
rz(0.176972) q[2];
rz(0.58961287) q[3];
sx q[3];
rz(-1.4859345) q[3];
sx q[3];
rz(-1.0346712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4347587) q[0];
sx q[0];
rz(-2.28573) q[0];
sx q[0];
rz(-0.51704299) q[0];
rz(-2.928858) q[1];
sx q[1];
rz(-1.0849625) q[1];
sx q[1];
rz(-0.83232602) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34236429) q[0];
sx q[0];
rz(-1.5629074) q[0];
sx q[0];
rz(-0.042021463) q[0];
rz(0.8950142) q[2];
sx q[2];
rz(-2.7594406) q[2];
sx q[2];
rz(2.7016751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7169818) q[1];
sx q[1];
rz(-1.2094133) q[1];
sx q[1];
rz(-0.34504621) q[1];
rz(-pi) q[2];
rz(-2.832507) q[3];
sx q[3];
rz(-0.49843542) q[3];
sx q[3];
rz(0.70685742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0011255) q[2];
sx q[2];
rz(-0.8320063) q[2];
sx q[2];
rz(2.7531085) q[2];
rz(1.2515986) q[3];
sx q[3];
rz(-1.3585217) q[3];
sx q[3];
rz(-2.0775011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851704) q[0];
sx q[0];
rz(-0.52260411) q[0];
sx q[0];
rz(-1.1473468) q[0];
rz(-1.7583678) q[1];
sx q[1];
rz(-1.5478094) q[1];
sx q[1];
rz(-2.9668818) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6679113) q[0];
sx q[0];
rz(-1.5690049) q[0];
sx q[0];
rz(0.67497336) q[0];
rz(-1.5632929) q[2];
sx q[2];
rz(-0.96133366) q[2];
sx q[2];
rz(0.96169072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6683674) q[1];
sx q[1];
rz(-0.36512926) q[1];
sx q[1];
rz(-2.5561246) q[1];
x q[2];
rz(-0.21214738) q[3];
sx q[3];
rz(-1.6039994) q[3];
sx q[3];
rz(-1.3591059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2061578) q[2];
sx q[2];
rz(-1.9400699) q[2];
sx q[2];
rz(-0.60714444) q[2];
rz(-0.28572765) q[3];
sx q[3];
rz(-0.61107475) q[3];
sx q[3];
rz(-0.40209517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6242591) q[0];
sx q[0];
rz(-1.3579955) q[0];
sx q[0];
rz(0.72878033) q[0];
rz(2.0023316) q[1];
sx q[1];
rz(-1.1092721) q[1];
sx q[1];
rz(-1.8702102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8761638) q[0];
sx q[0];
rz(-0.012152925) q[0];
sx q[0];
rz(1.3738213) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5554496) q[2];
sx q[2];
rz(-1.8019466) q[2];
sx q[2];
rz(1.4755206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1955586) q[1];
sx q[1];
rz(-0.6094774) q[1];
sx q[1];
rz(2.7835664) q[1];
x q[2];
rz(-1.3649105) q[3];
sx q[3];
rz(-0.48501462) q[3];
sx q[3];
rz(-0.58615326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83069688) q[2];
sx q[2];
rz(-1.4163914) q[2];
sx q[2];
rz(-2.8719416) q[2];
rz(-1.0722748) q[3];
sx q[3];
rz(-2.8282073) q[3];
sx q[3];
rz(-1.0543893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521745) q[0];
sx q[0];
rz(-1.1913238) q[0];
sx q[0];
rz(0.3840951) q[0];
rz(2.0317888) q[1];
sx q[1];
rz(-2.2566819) q[1];
sx q[1];
rz(-1.6222515) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42253112) q[0];
sx q[0];
rz(-1.5572955) q[0];
sx q[0];
rz(-0.11106981) q[0];
rz(0.39471925) q[2];
sx q[2];
rz(-0.92967722) q[2];
sx q[2];
rz(-1.2441928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1311595) q[1];
sx q[1];
rz(-1.466806) q[1];
sx q[1];
rz(1.3532182) q[1];
x q[2];
rz(0.48990814) q[3];
sx q[3];
rz(-1.7246609) q[3];
sx q[3];
rz(-2.7780848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9263837) q[2];
sx q[2];
rz(-1.0101725) q[2];
sx q[2];
rz(-2.9877648) q[2];
rz(-2.202863) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(1.4332829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7275823) q[0];
sx q[0];
rz(-1.753597) q[0];
sx q[0];
rz(-0.54186064) q[0];
rz(1.9210723) q[1];
sx q[1];
rz(-2.4668756) q[1];
sx q[1];
rz(3.0774679) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42196938) q[0];
sx q[0];
rz(-0.066264205) q[0];
sx q[0];
rz(0.55304773) q[0];
x q[1];
rz(1.7226965) q[2];
sx q[2];
rz(-2.4338401) q[2];
sx q[2];
rz(-0.47063702) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4628479) q[1];
sx q[1];
rz(-1.5539377) q[1];
sx q[1];
rz(-2.4900311) q[1];
rz(-pi) q[2];
rz(0.4145741) q[3];
sx q[3];
rz(-1.7547973) q[3];
sx q[3];
rz(-2.6340226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3037783) q[2];
sx q[2];
rz(-2.1393675) q[2];
sx q[2];
rz(3.0617867) q[2];
rz(0.59085733) q[3];
sx q[3];
rz(-0.30954027) q[3];
sx q[3];
rz(-0.045056067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87126842) q[0];
sx q[0];
rz(-2.7239983) q[0];
sx q[0];
rz(0.25750345) q[0];
rz(2.2389257) q[1];
sx q[1];
rz(-0.84044424) q[1];
sx q[1];
rz(-2.2534175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64110735) q[0];
sx q[0];
rz(-1.3543924) q[0];
sx q[0];
rz(-1.2016625) q[0];
x q[1];
rz(2.4492017) q[2];
sx q[2];
rz(-0.93442011) q[2];
sx q[2];
rz(2.6264735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0400454) q[1];
sx q[1];
rz(-2.2253621) q[1];
sx q[1];
rz(-2.1957869) q[1];
rz(-pi) q[2];
rz(0.8646263) q[3];
sx q[3];
rz(-0.99643222) q[3];
sx q[3];
rz(2.7780617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73457926) q[2];
sx q[2];
rz(-2.1740156) q[2];
sx q[2];
rz(0.61545294) q[2];
rz(-2.9493799) q[3];
sx q[3];
rz(-0.8513611) q[3];
sx q[3];
rz(-1.6063469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2007582) q[0];
sx q[0];
rz(-1.8466908) q[0];
sx q[0];
rz(1.2373663) q[0];
rz(-1.0279961) q[1];
sx q[1];
rz(-2.4528687) q[1];
sx q[1];
rz(-2.5797896) q[1];
rz(3.0801283) q[2];
sx q[2];
rz(-2.2066084) q[2];
sx q[2];
rz(1.6359272) q[2];
rz(1.2067704) q[3];
sx q[3];
rz(-2.0828928) q[3];
sx q[3];
rz(-0.29844007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
