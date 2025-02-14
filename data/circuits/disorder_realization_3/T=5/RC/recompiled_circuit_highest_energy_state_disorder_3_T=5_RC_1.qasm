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
rz(-3.1138104) q[0];
sx q[0];
rz(-1.8960928) q[0];
sx q[0];
rz(-1.0863289) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(1.2135222) q[1];
sx q[1];
rz(10.57099) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47482309) q[0];
sx q[0];
rz(-2.7753029) q[0];
sx q[0];
rz(-2.0369056) q[0];
x q[1];
rz(0.5733344) q[2];
sx q[2];
rz(-2.4681915) q[2];
sx q[2];
rz(-1.8512938) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1172197) q[1];
sx q[1];
rz(-2.0659323) q[1];
sx q[1];
rz(-0.58039033) q[1];
rz(-pi) q[2];
rz(-2.0729154) q[3];
sx q[3];
rz(-1.3716591) q[3];
sx q[3];
rz(-0.56875436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(-2.7794465) q[2];
rz(-0.95237887) q[3];
sx q[3];
rz(-2.5824472) q[3];
sx q[3];
rz(-2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.1700965) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(-1.6098518) q[0];
rz(-1.2486628) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-2.2889287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0097875) q[0];
sx q[0];
rz(-0.62776792) q[0];
sx q[0];
rz(-2.1814697) q[0];
x q[1];
rz(-3.018969) q[2];
sx q[2];
rz(-1.2957591) q[2];
sx q[2];
rz(-3.1240535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1601847) q[1];
sx q[1];
rz(-1.7716737) q[1];
sx q[1];
rz(2.1698879) q[1];
rz(-2.1308822) q[3];
sx q[3];
rz(-1.2398071) q[3];
sx q[3];
rz(0.50332848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93069211) q[2];
sx q[2];
rz(-2.1810668) q[2];
sx q[2];
rz(1.9453913) q[2];
rz(-0.026737468) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(-1.4727717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4836327) q[0];
sx q[0];
rz(-0.25068972) q[0];
sx q[0];
rz(-0.87580097) q[0];
rz(0.35975131) q[1];
sx q[1];
rz(-1.7774372) q[1];
sx q[1];
rz(-1.0035427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10570603) q[0];
sx q[0];
rz(-1.0125475) q[0];
sx q[0];
rz(-1.654675) q[0];
rz(1.0613173) q[2];
sx q[2];
rz(-0.70827863) q[2];
sx q[2];
rz(-0.88513155) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7454804) q[1];
sx q[1];
rz(-0.89876938) q[1];
sx q[1];
rz(1.7538096) q[1];
x q[2];
rz(-2.4963673) q[3];
sx q[3];
rz(-0.78515437) q[3];
sx q[3];
rz(0.19715362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2063736) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(0.75500542) q[2];
rz(3.0435009) q[3];
sx q[3];
rz(-2.5710929) q[3];
sx q[3];
rz(2.7437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(0.64495069) q[0];
rz(1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(-1.1114978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1350461) q[0];
sx q[0];
rz(-1.6301951) q[0];
sx q[0];
rz(0.063755449) q[0];
rz(-pi) q[1];
rz(0.45504163) q[2];
sx q[2];
rz(-0.18559449) q[2];
sx q[2];
rz(1.1374823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1810236) q[1];
sx q[1];
rz(-1.7180859) q[1];
sx q[1];
rz(2.9295314) q[1];
rz(-pi) q[2];
x q[2];
rz(2.901279) q[3];
sx q[3];
rz(-0.75551582) q[3];
sx q[3];
rz(-0.66080581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5081386) q[2];
sx q[2];
rz(-0.78712574) q[2];
sx q[2];
rz(2.656929) q[2];
rz(0.44114068) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.2539366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6345374) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(-1.7150568) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-2.8883002) q[1];
sx q[1];
rz(-2.3111129) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4035506) q[0];
sx q[0];
rz(-1.0071686) q[0];
sx q[0];
rz(0.72569816) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057382345) q[2];
sx q[2];
rz(-2.2028366) q[2];
sx q[2];
rz(-1.1626279) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75447202) q[1];
sx q[1];
rz(-0.74229147) q[1];
sx q[1];
rz(-0.20739389) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7375505) q[3];
sx q[3];
rz(-2.7528304) q[3];
sx q[3];
rz(2.5498912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44437235) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(-2.3946136) q[2];
rz(-0.30484453) q[3];
sx q[3];
rz(-1.7458785) q[3];
sx q[3];
rz(-1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1759258) q[0];
sx q[0];
rz(-1.3015863) q[0];
sx q[0];
rz(-2.9081705) q[0];
rz(-2.6790791) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(2.7375284) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83749002) q[0];
sx q[0];
rz(-1.7902052) q[0];
sx q[0];
rz(-0.93774937) q[0];
rz(-1.839371) q[2];
sx q[2];
rz(-1.982455) q[2];
sx q[2];
rz(-1.9869164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5517496) q[1];
sx q[1];
rz(-1.0115563) q[1];
sx q[1];
rz(-1.6980805) q[1];
rz(-pi) q[2];
rz(-0.29849739) q[3];
sx q[3];
rz(-2.4892471) q[3];
sx q[3];
rz(0.17078313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3290951) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(2.4883032) q[2];
rz(1.6434044) q[3];
sx q[3];
rz(-0.61763063) q[3];
sx q[3];
rz(0.60630715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46635982) q[0];
sx q[0];
rz(-1.9283858) q[0];
sx q[0];
rz(-0.32077041) q[0];
rz(0.65387154) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(-3.0810862) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.190028) q[0];
sx q[0];
rz(-1.4037644) q[0];
sx q[0];
rz(1.8442979) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94580067) q[2];
sx q[2];
rz(-2.201718) q[2];
sx q[2];
rz(0.45880246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6937032) q[1];
sx q[1];
rz(-0.9606908) q[1];
sx q[1];
rz(0.59887664) q[1];
rz(-2.2809847) q[3];
sx q[3];
rz(-0.51653701) q[3];
sx q[3];
rz(2.9252657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(-0.72236577) q[2];
rz(0.28313053) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(2.6319671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.79001456) q[0];
sx q[0];
rz(-2.5836662) q[0];
sx q[0];
rz(-2.956692) q[0];
rz(2.6376873) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(2.6769743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16581841) q[0];
sx q[0];
rz(-1.4041472) q[0];
sx q[0];
rz(-1.7383132) q[0];
rz(0.41733317) q[2];
sx q[2];
rz(-1.6586896) q[2];
sx q[2];
rz(-2.7869239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2459348) q[1];
sx q[1];
rz(-2.5634951) q[1];
sx q[1];
rz(-0.74191414) q[1];
x q[2];
rz(-2.9161386) q[3];
sx q[3];
rz(-0.73235336) q[3];
sx q[3];
rz(1.6868397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9087312) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(-2.0880584) q[2];
rz(-0.72707027) q[3];
sx q[3];
rz(-1.9126242) q[3];
sx q[3];
rz(2.5963636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19186774) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(0.18549347) q[0];
rz(0.56974757) q[1];
sx q[1];
rz(-0.92420095) q[1];
sx q[1];
rz(1.9843019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7737024) q[0];
sx q[0];
rz(-1.7130525) q[0];
sx q[0];
rz(0.41928798) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6389388) q[2];
sx q[2];
rz(-0.63387094) q[2];
sx q[2];
rz(2.9958452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2409199) q[1];
sx q[1];
rz(-1.0176588) q[1];
sx q[1];
rz(-1.8227804) q[1];
rz(0.80482153) q[3];
sx q[3];
rz(-1.2535044) q[3];
sx q[3];
rz(2.1498093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40413228) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(1.4004716) q[2];
rz(-0.36462668) q[3];
sx q[3];
rz(-0.76321634) q[3];
sx q[3];
rz(-2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51065651) q[0];
sx q[0];
rz(-2.3031504) q[0];
sx q[0];
rz(3.1070218) q[0];
rz(0.336054) q[1];
sx q[1];
rz(-2.2269109) q[1];
sx q[1];
rz(-2.5539982) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60025763) q[0];
sx q[0];
rz(-0.44352874) q[0];
sx q[0];
rz(0.33716538) q[0];
rz(2.7828092) q[2];
sx q[2];
rz(-2.3780895) q[2];
sx q[2];
rz(2.7051008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42660824) q[1];
sx q[1];
rz(-2.9054925) q[1];
sx q[1];
rz(-1.7515504) q[1];
x q[2];
rz(-2.5925007) q[3];
sx q[3];
rz(-0.7585511) q[3];
sx q[3];
rz(-0.18699924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8079638) q[2];
sx q[2];
rz(-2.2449292) q[2];
sx q[2];
rz(0.72820747) q[2];
rz(2.7306469) q[3];
sx q[3];
rz(-2.9185037) q[3];
sx q[3];
rz(1.8691011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644607) q[0];
sx q[0];
rz(-1.5560173) q[0];
sx q[0];
rz(-1.7898855) q[0];
rz(-2.7017055) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(-0.71285474) q[2];
sx q[2];
rz(-0.44217449) q[2];
sx q[2];
rz(0.90978734) q[2];
rz(-0.41070134) q[3];
sx q[3];
rz(-0.11631415) q[3];
sx q[3];
rz(0.61396413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
