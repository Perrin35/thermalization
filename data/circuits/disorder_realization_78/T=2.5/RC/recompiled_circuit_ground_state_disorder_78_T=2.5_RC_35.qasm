OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4888157) q[0];
sx q[0];
rz(-0.96537679) q[0];
sx q[0];
rz(2.2247347) q[0];
rz(1.4404453) q[1];
sx q[1];
rz(2.1135766) q[1];
sx q[1];
rz(10.120575) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.643508) q[0];
sx q[0];
rz(-1.1445481) q[0];
sx q[0];
rz(-3.0820937) q[0];
rz(-pi) q[1];
rz(2.2356625) q[2];
sx q[2];
rz(-3.0810789) q[2];
sx q[2];
rz(-2.0384332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6219139) q[1];
sx q[1];
rz(-1.7413472) q[1];
sx q[1];
rz(2.7615553) q[1];
rz(-pi) q[2];
rz(-0.37000044) q[3];
sx q[3];
rz(-1.1048855) q[3];
sx q[3];
rz(2.6107058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15771389) q[2];
sx q[2];
rz(-2.3215051) q[2];
sx q[2];
rz(0.90041089) q[2];
rz(-2.3484717) q[3];
sx q[3];
rz(-1.6499062) q[3];
sx q[3];
rz(0.69860631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124436) q[0];
sx q[0];
rz(-1.1684893) q[0];
sx q[0];
rz(2.4560665) q[0];
rz(0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(-1.3139075) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066931574) q[0];
sx q[0];
rz(-1.0767796) q[0];
sx q[0];
rz(-0.95290174) q[0];
rz(-pi) q[1];
rz(2.3150748) q[2];
sx q[2];
rz(-1.2213301) q[2];
sx q[2];
rz(0.068880388) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5470299) q[1];
sx q[1];
rz(-2.4395151) q[1];
sx q[1];
rz(2.3365993) q[1];
rz(-pi) q[2];
rz(-2.7502188) q[3];
sx q[3];
rz(-2.5487566) q[3];
sx q[3];
rz(0.22030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1757397) q[2];
sx q[2];
rz(-0.41030914) q[2];
sx q[2];
rz(-2.8941594) q[2];
rz(-0.11582173) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(2.5628832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9847357) q[0];
sx q[0];
rz(-2.7655089) q[0];
sx q[0];
rz(2.0767186) q[0];
rz(-0.6330601) q[1];
sx q[1];
rz(-1.9828911) q[1];
sx q[1];
rz(0.13654581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2257803) q[0];
sx q[0];
rz(-1.5951711) q[0];
sx q[0];
rz(-0.73631411) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34935475) q[2];
sx q[2];
rz(-2.0094479) q[2];
sx q[2];
rz(0.037540066) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3202793) q[1];
sx q[1];
rz(-1.1007231) q[1];
sx q[1];
rz(2.4072985) q[1];
rz(-pi) q[2];
rz(-0.8053369) q[3];
sx q[3];
rz(-2.5743765) q[3];
sx q[3];
rz(0.50931168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.442405) q[2];
sx q[2];
rz(-1.2523315) q[2];
sx q[2];
rz(1.5999954) q[2];
rz(1.0934746) q[3];
sx q[3];
rz(-1.5019865) q[3];
sx q[3];
rz(-2.7961965) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4271451) q[0];
sx q[0];
rz(-2.8801425) q[0];
sx q[0];
rz(2.5647822) q[0];
rz(-2.414074) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(2.2732546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6431516) q[0];
sx q[0];
rz(-0.77385573) q[0];
sx q[0];
rz(-2.5572148) q[0];
rz(-pi) q[1];
rz(-2.6170591) q[2];
sx q[2];
rz(-2.0887592) q[2];
sx q[2];
rz(0.7690767) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5975004) q[1];
sx q[1];
rz(-0.75832931) q[1];
sx q[1];
rz(2.6306398) q[1];
x q[2];
rz(2.2277539) q[3];
sx q[3];
rz(-2.4025177) q[3];
sx q[3];
rz(-0.55348671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19520983) q[2];
sx q[2];
rz(-1.752744) q[2];
sx q[2];
rz(2.4460804) q[2];
rz(0.21466151) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.6014086) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(-1.1001128) q[0];
rz(0.21057883) q[1];
sx q[1];
rz(-2.5081878) q[1];
sx q[1];
rz(1.6483773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8755668) q[0];
sx q[0];
rz(-1.1023836) q[0];
sx q[0];
rz(2.1126956) q[0];
x q[1];
rz(1.1334992) q[2];
sx q[2];
rz(-2.9455212) q[2];
sx q[2];
rz(2.9626949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82147089) q[1];
sx q[1];
rz(-1.4152598) q[1];
sx q[1];
rz(2.0633153) q[1];
x q[2];
rz(-2.8272618) q[3];
sx q[3];
rz(-2.4498457) q[3];
sx q[3];
rz(-1.421613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7275927) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(2.2445938) q[2];
rz(1.2796848) q[3];
sx q[3];
rz(-2.1309659) q[3];
sx q[3];
rz(0.047209386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419264) q[0];
sx q[0];
rz(-0.62586313) q[0];
sx q[0];
rz(2.8114787) q[0];
rz(2.5008028) q[1];
sx q[1];
rz(-1.375744) q[1];
sx q[1];
rz(-0.97553387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9331688) q[0];
sx q[0];
rz(-2.5481249) q[0];
sx q[0];
rz(0.10153254) q[0];
x q[1];
rz(2.8728666) q[2];
sx q[2];
rz(-1.8789504) q[2];
sx q[2];
rz(1.2040963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9793518) q[1];
sx q[1];
rz(-1.8514412) q[1];
sx q[1];
rz(-0.45307093) q[1];
rz(-pi) q[2];
rz(1.1694856) q[3];
sx q[3];
rz(-0.45404497) q[3];
sx q[3];
rz(1.3523449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2834125) q[2];
sx q[2];
rz(-0.25412574) q[2];
sx q[2];
rz(-0.88968366) q[2];
rz(-2.5455918) q[3];
sx q[3];
rz(-2.0723876) q[3];
sx q[3];
rz(-0.10665756) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058873) q[0];
sx q[0];
rz(-2.0054673) q[0];
sx q[0];
rz(1.1329875) q[0];
rz(0.63944447) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(-0.9674255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15617642) q[0];
sx q[0];
rz(-1.7207017) q[0];
sx q[0];
rz(1.0557337) q[0];
x q[1];
rz(-0.70761159) q[2];
sx q[2];
rz(-2.65776) q[2];
sx q[2];
rz(0.31602032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6220807) q[1];
sx q[1];
rz(-1.4556307) q[1];
sx q[1];
rz(2.9154816) q[1];
rz(0.65770517) q[3];
sx q[3];
rz(-1.8175565) q[3];
sx q[3];
rz(1.2031499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3966668) q[2];
sx q[2];
rz(-1.0090088) q[2];
sx q[2];
rz(-0.94433707) q[2];
rz(0.62263292) q[3];
sx q[3];
rz(-2.1532652) q[3];
sx q[3];
rz(-0.78701377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-3.0866163) q[0];
sx q[0];
rz(-2.3353307) q[0];
sx q[0];
rz(3.1373613) q[0];
rz(1.6731693) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(-0.17446336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3136351) q[0];
sx q[0];
rz(-2.9330755) q[0];
sx q[0];
rz(-2.0329518) q[0];
rz(-pi) q[1];
rz(0.44850839) q[2];
sx q[2];
rz(-1.4786198) q[2];
sx q[2];
rz(2.0192469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9618865) q[1];
sx q[1];
rz(-2.440036) q[1];
sx q[1];
rz(1.4683734) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.065062398) q[3];
sx q[3];
rz(-1.3351249) q[3];
sx q[3];
rz(0.5103569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2506432) q[2];
sx q[2];
rz(-0.093891056) q[2];
sx q[2];
rz(-2.6012361) q[2];
rz(-3.12449) q[3];
sx q[3];
rz(-2.159724) q[3];
sx q[3];
rz(2.4798149) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509725) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(-2.354994) q[0];
rz(2.0358918) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(-0.22122637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.008667058) q[0];
sx q[0];
rz(-1.2930388) q[0];
sx q[0];
rz(0.8322258) q[0];
rz(-pi) q[1];
rz(1.3972261) q[2];
sx q[2];
rz(-1.4772869) q[2];
sx q[2];
rz(-1.6261404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58824268) q[1];
sx q[1];
rz(-2.6094631) q[1];
sx q[1];
rz(-0.41779907) q[1];
x q[2];
rz(0.28085162) q[3];
sx q[3];
rz(-2.4316868) q[3];
sx q[3];
rz(0.8381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9959539) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(-1.9560187) q[2];
rz(2.0354185) q[3];
sx q[3];
rz(-1.0217051) q[3];
sx q[3];
rz(-0.98418981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7800061) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(0.7097882) q[0];
rz(2.2478726) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(-2.6745904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119962) q[0];
sx q[0];
rz(-1.3043367) q[0];
sx q[0];
rz(-0.14273739) q[0];
x q[1];
rz(-2.8728629) q[2];
sx q[2];
rz(-2.1788283) q[2];
sx q[2];
rz(-3.1399474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2645006) q[1];
sx q[1];
rz(-2.9000726) q[1];
sx q[1];
rz(1.7272374) q[1];
x q[2];
rz(-2.267546) q[3];
sx q[3];
rz(-1.0206887) q[3];
sx q[3];
rz(2.3693905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3422157) q[2];
sx q[2];
rz(-2.2365384) q[2];
sx q[2];
rz(0.13599914) q[2];
rz(-0.32986766) q[3];
sx q[3];
rz(-2.0743275) q[3];
sx q[3];
rz(-0.29712591) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8316523) q[0];
sx q[0];
rz(-1.1158054) q[0];
sx q[0];
rz(0.62624395) q[0];
rz(0.88824875) q[1];
sx q[1];
rz(-1.2445969) q[1];
sx q[1];
rz(0.1319763) q[1];
rz(1.3462979) q[2];
sx q[2];
rz(-1.1830876) q[2];
sx q[2];
rz(3.0814688) q[2];
rz(0.72867758) q[3];
sx q[3];
rz(-1.8986438) q[3];
sx q[3];
rz(1.8713634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
