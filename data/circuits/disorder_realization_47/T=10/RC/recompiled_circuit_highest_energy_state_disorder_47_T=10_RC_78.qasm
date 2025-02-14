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
rz(1.0951618) q[0];
sx q[0];
rz(-2.996063) q[0];
sx q[0];
rz(2.6382883) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(1.4454747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70899773) q[0];
sx q[0];
rz(-1.1721969) q[0];
sx q[0];
rz(1.523827) q[0];
x q[1];
rz(-0.66464773) q[2];
sx q[2];
rz(-0.87021135) q[2];
sx q[2];
rz(1.0657016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18516416) q[1];
sx q[1];
rz(-1.899029) q[1];
sx q[1];
rz(-0.64394711) q[1];
rz(-pi) q[2];
rz(2.8347391) q[3];
sx q[3];
rz(-3.0899991) q[3];
sx q[3];
rz(1.8966228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6426927) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(-1.6373681) q[2];
rz(-0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7351643) q[0];
sx q[0];
rz(-2.6728215) q[0];
sx q[0];
rz(0.51097393) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.7402486) q[1];
sx q[1];
rz(0.32646349) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8998979) q[0];
sx q[0];
rz(-2.7621671) q[0];
sx q[0];
rz(-2.3510758) q[0];
rz(-pi) q[1];
rz(-0.20997491) q[2];
sx q[2];
rz(-0.89824235) q[2];
sx q[2];
rz(-0.58174101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4703456) q[1];
sx q[1];
rz(-1.874028) q[1];
sx q[1];
rz(0.95945759) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9580919) q[3];
sx q[3];
rz(-1.3069469) q[3];
sx q[3];
rz(-1.6843759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9713126) q[2];
sx q[2];
rz(-0.22394094) q[2];
sx q[2];
rz(1.8453321) q[2];
rz(-2.7235459) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(-2.4372098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0728077) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(-0.54247722) q[0];
rz(1.7659278) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(3.1105522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2707996) q[0];
sx q[0];
rz(-2.7794771) q[0];
sx q[0];
rz(-2.5848081) q[0];
rz(-2.9152053) q[2];
sx q[2];
rz(-2.2157266) q[2];
sx q[2];
rz(1.3340536) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87796383) q[1];
sx q[1];
rz(-0.97743703) q[1];
sx q[1];
rz(-0.29747648) q[1];
x q[2];
rz(3.1150041) q[3];
sx q[3];
rz(-0.55517653) q[3];
sx q[3];
rz(-2.7257277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5151908) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(2.314563) q[2];
rz(2.6229897) q[3];
sx q[3];
rz(-2.0293472) q[3];
sx q[3];
rz(1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.9943635) q[0];
sx q[0];
rz(-1.8356859) q[0];
sx q[0];
rz(-1.9299141) q[0];
rz(-2.0934824) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(-1.3406219) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692684) q[0];
sx q[0];
rz(-0.71592531) q[0];
sx q[0];
rz(-2.6628519) q[0];
rz(-1.421809) q[2];
sx q[2];
rz(-1.3366615) q[2];
sx q[2];
rz(2.9643167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1010434) q[1];
sx q[1];
rz(-1.2770997) q[1];
sx q[1];
rz(2.5459176) q[1];
x q[2];
rz(2.8302835) q[3];
sx q[3];
rz(-1.177586) q[3];
sx q[3];
rz(1.1953199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-2.9643855) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.6107669) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(0.78260666) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57946) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(2.539047) q[0];
rz(0.58019477) q[1];
sx q[1];
rz(-1.5126901) q[1];
sx q[1];
rz(-2.3604732) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21602042) q[0];
sx q[0];
rz(-1.7402152) q[0];
sx q[0];
rz(1.8623307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4915472) q[2];
sx q[2];
rz(-2.263265) q[2];
sx q[2];
rz(1.6672857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.994532) q[1];
sx q[1];
rz(-0.72685469) q[1];
sx q[1];
rz(1.4406073) q[1];
rz(0.1991462) q[3];
sx q[3];
rz(-2.4059787) q[3];
sx q[3];
rz(-2.3922331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96131229) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(-2.9126634) q[2];
rz(-1.8217575) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(-0.84407097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9127386) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(1.9629021) q[0];
rz(1.9121869) q[1];
sx q[1];
rz(-2.7472159) q[1];
sx q[1];
rz(1.2961402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154436) q[0];
sx q[0];
rz(-2.8279732) q[0];
sx q[0];
rz(-1.1152033) q[0];
x q[1];
rz(-2.4086094) q[2];
sx q[2];
rz(-0.99501164) q[2];
sx q[2];
rz(-0.95966649) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.072159616) q[1];
sx q[1];
rz(-2.0516703) q[1];
sx q[1];
rz(2.062501) q[1];
x q[2];
rz(-0.48719897) q[3];
sx q[3];
rz(-1.4811084) q[3];
sx q[3];
rz(1.7023757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.487315) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(0.26228341) q[2];
rz(2.8413963) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(-2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009785) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(-1.8390919) q[0];
rz(0.66849661) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(1.0350593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9205017) q[0];
sx q[0];
rz(-1.1701487) q[0];
sx q[0];
rz(-0.23200881) q[0];
rz(-1.2012173) q[2];
sx q[2];
rz(-1.7453004) q[2];
sx q[2];
rz(-2.7480887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3429195) q[1];
sx q[1];
rz(-1.2274776) q[1];
sx q[1];
rz(0.080282057) q[1];
x q[2];
rz(2.3978698) q[3];
sx q[3];
rz(-1.111314) q[3];
sx q[3];
rz(-0.93839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0237026) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(-1.6560076) q[2];
rz(-2.8591136) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3561803) q[0];
sx q[0];
rz(-2.0996576) q[0];
sx q[0];
rz(0.27454141) q[0];
rz(1.9369594) q[1];
sx q[1];
rz(-2.5614673) q[1];
sx q[1];
rz(2.5801632) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954079) q[0];
sx q[0];
rz(-2.2642379) q[0];
sx q[0];
rz(-2.2827143) q[0];
x q[1];
rz(2.4653685) q[2];
sx q[2];
rz(-1.3315787) q[2];
sx q[2];
rz(-2.4382811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9778487) q[1];
sx q[1];
rz(-1.9789961) q[1];
sx q[1];
rz(-2.8099437) q[1];
rz(-2.8772914) q[3];
sx q[3];
rz(-1.2339379) q[3];
sx q[3];
rz(1.5764396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56208912) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(1.4037464) q[2];
rz(1.7956519) q[3];
sx q[3];
rz(-2.3965049) q[3];
sx q[3];
rz(-2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76029921) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(0.90989939) q[0];
rz(1.9445885) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(0.88814703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67535454) q[0];
sx q[0];
rz(-2.9976237) q[0];
sx q[0];
rz(-2.2813802) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65983628) q[2];
sx q[2];
rz(-1.5232695) q[2];
sx q[2];
rz(0.30872503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8348114) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(-2.4624834) q[1];
rz(0.03753438) q[3];
sx q[3];
rz(-1.2753295) q[3];
sx q[3];
rz(-0.42643828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0137279) q[2];
sx q[2];
rz(-1.4415386) q[2];
sx q[2];
rz(-1.0825276) q[2];
rz(-2.9108289) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312663) q[0];
sx q[0];
rz(-2.3895097) q[0];
sx q[0];
rz(0.58151522) q[0];
rz(2.5260018) q[1];
sx q[1];
rz(-0.94771996) q[1];
sx q[1];
rz(-1.8448578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55330861) q[0];
sx q[0];
rz(-0.020961449) q[0];
sx q[0];
rz(-2.0901438) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7974772) q[2];
sx q[2];
rz(-1.5798347) q[2];
sx q[2];
rz(1.2492385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77380005) q[1];
sx q[1];
rz(-1.3145216) q[1];
sx q[1];
rz(2.0092177) q[1];
rz(1.9711791) q[3];
sx q[3];
rz(-0.93850157) q[3];
sx q[3];
rz(-0.99611547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7221308) q[2];
sx q[2];
rz(-3.0007134) q[2];
sx q[2];
rz(-2.6922743) q[2];
rz(-0.32014534) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(-1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.768059) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(-1.3846579) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(-1.2133219) q[2];
sx q[2];
rz(-0.64897474) q[2];
sx q[2];
rz(-2.8233768) q[2];
rz(2.2438335) q[3];
sx q[3];
rz(-0.57716341) q[3];
sx q[3];
rz(-2.9152277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
