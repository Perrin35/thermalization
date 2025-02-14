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
rz(-1.8972549) q[0];
sx q[0];
rz(-2.3490348) q[0];
sx q[0];
rz(-2.8886524) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(-2.7546496) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16432504) q[0];
sx q[0];
rz(-1.8935793) q[0];
sx q[0];
rz(-0.27077814) q[0];
rz(0.57780452) q[2];
sx q[2];
rz(-1.4981604) q[2];
sx q[2];
rz(-0.60608038) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36648892) q[1];
sx q[1];
rz(-1.6548043) q[1];
sx q[1];
rz(-1.4909527) q[1];
rz(1.0798595) q[3];
sx q[3];
rz(-1.0359952) q[3];
sx q[3];
rz(0.25568572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.046254961) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(1.7389899) q[2];
rz(1.4143573) q[3];
sx q[3];
rz(-2.4284913) q[3];
sx q[3];
rz(0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.6801179) q[0];
sx q[0];
rz(-0.48826996) q[0];
sx q[0];
rz(0.71037355) q[0];
rz(2.12517) q[1];
sx q[1];
rz(-1.2998394) q[1];
sx q[1];
rz(2.3764835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0046526) q[0];
sx q[0];
rz(-2.1015133) q[0];
sx q[0];
rz(-1.1999446) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4508453) q[2];
sx q[2];
rz(-1.0260858) q[2];
sx q[2];
rz(0.70554698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5580235) q[1];
sx q[1];
rz(-0.92881948) q[1];
sx q[1];
rz(2.7318673) q[1];
rz(-pi) q[2];
rz(0.14987544) q[3];
sx q[3];
rz(-1.2399763) q[3];
sx q[3];
rz(-0.56043032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1171099) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(0.70466858) q[2];
rz(-0.17413983) q[3];
sx q[3];
rz(-2.2691059) q[3];
sx q[3];
rz(-2.7599938) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0891377) q[0];
sx q[0];
rz(-1.313504) q[0];
sx q[0];
rz(2.254159) q[0];
rz(-1.2499836) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(1.3608305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3583547) q[0];
sx q[0];
rz(-1.1941507) q[0];
sx q[0];
rz(-2.5953617) q[0];
rz(-pi) q[1];
rz(-1.7086231) q[2];
sx q[2];
rz(-1.6500452) q[2];
sx q[2];
rz(0.85286959) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8831142) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(1.4264631) q[1];
rz(-pi) q[2];
rz(1.9599846) q[3];
sx q[3];
rz(-1.6103611) q[3];
sx q[3];
rz(0.54007441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7736194) q[2];
sx q[2];
rz(-2.7757288) q[2];
sx q[2];
rz(-1.6486637) q[2];
rz(-2.6922373) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(-1.9213283) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0538977) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(0.98211163) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(-1.7353479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236693) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(-1.7864321) q[0];
rz(-pi) q[1];
rz(2.2135229) q[2];
sx q[2];
rz(-0.065970369) q[2];
sx q[2];
rz(0.82974354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92038735) q[1];
sx q[1];
rz(-1.3246598) q[1];
sx q[1];
rz(2.1046776) q[1];
x q[2];
rz(-2.8729183) q[3];
sx q[3];
rz(-0.93963913) q[3];
sx q[3];
rz(2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1534319) q[2];
sx q[2];
rz(-1.016523) q[2];
sx q[2];
rz(-0.15596685) q[2];
rz(-2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-3.0273738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0680189) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(2.9408348) q[0];
rz(-1.9643895) q[1];
sx q[1];
rz(-0.8526082) q[1];
sx q[1];
rz(1.9815365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070712) q[0];
sx q[0];
rz(-1.7385027) q[0];
sx q[0];
rz(1.3434065) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54520901) q[2];
sx q[2];
rz(-3.0070947) q[2];
sx q[2];
rz(-2.4328977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39587089) q[1];
sx q[1];
rz(-1.9211287) q[1];
sx q[1];
rz(-1.8941903) q[1];
rz(-pi) q[2];
rz(1.0409768) q[3];
sx q[3];
rz(-1.3966832) q[3];
sx q[3];
rz(-1.8903738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3843627) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(-2.6833656) q[2];
rz(0.630817) q[3];
sx q[3];
rz(-1.2354555) q[3];
sx q[3];
rz(-0.49921504) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.683627) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-3.1347347) q[0];
rz(0.231617) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(-1.0622271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0322965) q[0];
sx q[0];
rz(-1.2589021) q[0];
sx q[0];
rz(-1.4410785) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7659645) q[2];
sx q[2];
rz(-0.53002073) q[2];
sx q[2];
rz(-2.9935868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4783027) q[1];
sx q[1];
rz(-1.8243102) q[1];
sx q[1];
rz(0.67914825) q[1];
rz(-pi) q[2];
rz(2.9763664) q[3];
sx q[3];
rz(-2.2731785) q[3];
sx q[3];
rz(-3.0000697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0007533) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(-1.2678649) q[2];
rz(0.86152348) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2086901) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-2.9260337) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(0.15942474) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529786) q[0];
sx q[0];
rz(-0.63487494) q[0];
sx q[0];
rz(1.5663516) q[0];
rz(-0.36783571) q[2];
sx q[2];
rz(-2.6421037) q[2];
sx q[2];
rz(0.50423451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21531184) q[1];
sx q[1];
rz(-2.9027853) q[1];
sx q[1];
rz(1.396149) q[1];
x q[2];
rz(-2.255329) q[3];
sx q[3];
rz(-0.66826398) q[3];
sx q[3];
rz(1.9192426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40519199) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(0.51652017) q[2];
rz(1.2107595) q[3];
sx q[3];
rz(-1.6885933) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844168) q[0];
sx q[0];
rz(-3.0653937) q[0];
sx q[0];
rz(-0.54022378) q[0];
rz(-1.5243439) q[1];
sx q[1];
rz(-2.1397619) q[1];
sx q[1];
rz(1.2670955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0277242) q[0];
sx q[0];
rz(-2.0696215) q[0];
sx q[0];
rz(-1.8201102) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3643963) q[2];
sx q[2];
rz(-1.0269477) q[2];
sx q[2];
rz(2.0573307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6798576) q[1];
sx q[1];
rz(-1.9377794) q[1];
sx q[1];
rz(-1.7253897) q[1];
rz(-0.9293988) q[3];
sx q[3];
rz(-1.2817973) q[3];
sx q[3];
rz(2.3232587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8871062) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(-2.9642588) q[2];
rz(-2.0717428) q[3];
sx q[3];
rz(-1.0515352) q[3];
sx q[3];
rz(0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808481) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(-0.10072197) q[0];
rz(-1.992647) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(1.1740059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3927395) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(-1.8003182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3753433) q[2];
sx q[2];
rz(-1.6403997) q[2];
sx q[2];
rz(2.0947667) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8496597) q[1];
sx q[1];
rz(-1.5330452) q[1];
sx q[1];
rz(0.019456073) q[1];
rz(-pi) q[2];
rz(-0.10281201) q[3];
sx q[3];
rz(-2.2592756) q[3];
sx q[3];
rz(1.4172872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9775057) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(2.738415) q[2];
rz(-0.66344231) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(-3.1373533) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4106301) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(0.078911111) q[0];
rz(0.90947378) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(-0.39631072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3003259) q[0];
sx q[0];
rz(-1.5551651) q[0];
sx q[0];
rz(0.0090740694) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.134507) q[2];
sx q[2];
rz(-0.79372294) q[2];
sx q[2];
rz(-2.7929579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7381607) q[1];
sx q[1];
rz(-2.6034174) q[1];
sx q[1];
rz(2.6713283) q[1];
rz(-pi) q[2];
rz(-2.2811398) q[3];
sx q[3];
rz(-0.75569587) q[3];
sx q[3];
rz(-2.6967654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69184715) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(-2.3425102) q[2];
rz(-3.0633022) q[3];
sx q[3];
rz(-1.2543863) q[3];
sx q[3];
rz(2.4962795) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66452022) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(-1.1935344) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(-0.1706201) q[2];
sx q[2];
rz(-1.1191875) q[2];
sx q[2];
rz(2.1714581) q[2];
rz(0.54936784) q[3];
sx q[3];
rz(-0.80437575) q[3];
sx q[3];
rz(1.5301216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
