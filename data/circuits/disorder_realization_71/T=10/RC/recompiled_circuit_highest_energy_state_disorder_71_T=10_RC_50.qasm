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
rz(3.0687301) q[0];
sx q[0];
rz(-1.0412359) q[0];
sx q[0];
rz(-0.6676724) q[0];
rz(3.0296037) q[1];
sx q[1];
rz(-1.4854687) q[1];
sx q[1];
rz(2.8356584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65016215) q[0];
sx q[0];
rz(-1.1592277) q[0];
sx q[0];
rz(1.6603907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0245695) q[2];
sx q[2];
rz(-1.0033022) q[2];
sx q[2];
rz(-0.62376444) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48253912) q[1];
sx q[1];
rz(-1.6371563) q[1];
sx q[1];
rz(-0.37239667) q[1];
rz(-pi) q[2];
rz(-0.70641993) q[3];
sx q[3];
rz(-1.6230246) q[3];
sx q[3];
rz(-2.4097648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2970807) q[2];
sx q[2];
rz(-1.7252555) q[2];
sx q[2];
rz(2.1616006) q[2];
rz(-2.8196107) q[3];
sx q[3];
rz(-1.201509) q[3];
sx q[3];
rz(2.9856248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148773) q[0];
sx q[0];
rz(-1.5276271) q[0];
sx q[0];
rz(2.2514586) q[0];
rz(1.1164249) q[1];
sx q[1];
rz(-2.2869488) q[1];
sx q[1];
rz(1.1847624) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6044004) q[0];
sx q[0];
rz(-1.0774776) q[0];
sx q[0];
rz(-1.7866578) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0204732) q[2];
sx q[2];
rz(-1.2492078) q[2];
sx q[2];
rz(0.85105291) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1175307) q[1];
sx q[1];
rz(-1.4595929) q[1];
sx q[1];
rz(-3.0255513) q[1];
rz(-pi) q[2];
rz(1.5873648) q[3];
sx q[3];
rz(-1.145716) q[3];
sx q[3];
rz(-2.1147902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6637791) q[2];
sx q[2];
rz(-0.16872831) q[2];
sx q[2];
rz(-1.7185877) q[2];
rz(-2.6723828) q[3];
sx q[3];
rz(-1.4956632) q[3];
sx q[3];
rz(2.7928228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94217268) q[0];
sx q[0];
rz(-1.9439789) q[0];
sx q[0];
rz(2.0592101) q[0];
rz(0.17467817) q[1];
sx q[1];
rz(-1.7593316) q[1];
sx q[1];
rz(-2.6686525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4648165) q[0];
sx q[0];
rz(-1.7308141) q[0];
sx q[0];
rz(-1.141044) q[0];
x q[1];
rz(-0.86955118) q[2];
sx q[2];
rz(-1.8461178) q[2];
sx q[2];
rz(0.67050951) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14607695) q[1];
sx q[1];
rz(-2.5716166) q[1];
sx q[1];
rz(1.6126627) q[1];
rz(-pi) q[2];
rz(0.93227902) q[3];
sx q[3];
rz(-0.74018708) q[3];
sx q[3];
rz(0.96961752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6579154) q[2];
sx q[2];
rz(-1.4853073) q[2];
sx q[2];
rz(1.8243054) q[2];
rz(-0.45051908) q[3];
sx q[3];
rz(-0.91166383) q[3];
sx q[3];
rz(1.7821504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6555742) q[0];
sx q[0];
rz(-1.065932) q[0];
sx q[0];
rz(-2.4701212) q[0];
rz(0.94054049) q[1];
sx q[1];
rz(-0.42010072) q[1];
sx q[1];
rz(2.1817575) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.335054) q[0];
sx q[0];
rz(-1.949652) q[0];
sx q[0];
rz(-0.23083487) q[0];
rz(-pi) q[1];
rz(0.28428712) q[2];
sx q[2];
rz(-2.4552567) q[2];
sx q[2];
rz(-0.057675408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18273045) q[1];
sx q[1];
rz(-0.36793837) q[1];
sx q[1];
rz(0.63506215) q[1];
x q[2];
rz(-1.6790326) q[3];
sx q[3];
rz(-2.0104791) q[3];
sx q[3];
rz(-0.81793438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81022108) q[2];
sx q[2];
rz(-1.081531) q[2];
sx q[2];
rz(-0.01037154) q[2];
rz(-1.8593908) q[3];
sx q[3];
rz(-0.60716647) q[3];
sx q[3];
rz(2.8289294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4717167) q[0];
sx q[0];
rz(-0.80729055) q[0];
sx q[0];
rz(-0.5759936) q[0];
rz(-2.5895789) q[1];
sx q[1];
rz(-1.8355398) q[1];
sx q[1];
rz(0.920151) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3993307) q[0];
sx q[0];
rz(-0.98679155) q[0];
sx q[0];
rz(-0.76764099) q[0];
rz(-pi) q[1];
rz(2.8497177) q[2];
sx q[2];
rz(-1.5988962) q[2];
sx q[2];
rz(-0.71189881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11959521) q[1];
sx q[1];
rz(-1.9174077) q[1];
sx q[1];
rz(-2.892811) q[1];
x q[2];
rz(-2.4534485) q[3];
sx q[3];
rz(-2.9811908) q[3];
sx q[3];
rz(-0.80996603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0612023) q[2];
sx q[2];
rz(-1.7549425) q[2];
sx q[2];
rz(2.6913255) q[2];
rz(3.0007512) q[3];
sx q[3];
rz(-1.2248657) q[3];
sx q[3];
rz(2.5950477) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30215728) q[0];
sx q[0];
rz(-1.7385229) q[0];
sx q[0];
rz(-0.19254011) q[0];
rz(-0.64388609) q[1];
sx q[1];
rz(-1.5778678) q[1];
sx q[1];
rz(3.0582757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57184982) q[0];
sx q[0];
rz(-2.5170287) q[0];
sx q[0];
rz(-1.3032105) q[0];
rz(3.0599454) q[2];
sx q[2];
rz(-1.6892216) q[2];
sx q[2];
rz(-1.7497964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6308438) q[1];
sx q[1];
rz(-0.35707316) q[1];
sx q[1];
rz(-2.0791173) q[1];
rz(2.5531876) q[3];
sx q[3];
rz(-2.6987966) q[3];
sx q[3];
rz(-1.9857085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4087499) q[2];
sx q[2];
rz(-2.7538444) q[2];
sx q[2];
rz(0.55633083) q[2];
rz(-2.2522816) q[3];
sx q[3];
rz(-1.6272864) q[3];
sx q[3];
rz(0.48243943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675875) q[0];
sx q[0];
rz(-2.6370625) q[0];
sx q[0];
rz(2.0489847) q[0];
rz(-2.4448474) q[1];
sx q[1];
rz(-2.4431591) q[1];
sx q[1];
rz(-2.4901857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4428042) q[0];
sx q[0];
rz(-0.68969986) q[0];
sx q[0];
rz(2.882363) q[0];
rz(-2.2397141) q[2];
sx q[2];
rz(-1.9142391) q[2];
sx q[2];
rz(0.65367389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16800616) q[1];
sx q[1];
rz(-1.028722) q[1];
sx q[1];
rz(-2.8219221) q[1];
x q[2];
rz(1.2907372) q[3];
sx q[3];
rz(-0.055563888) q[3];
sx q[3];
rz(-2.0357214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6394028) q[2];
sx q[2];
rz(-2.4263589) q[2];
sx q[2];
rz(-0.56360847) q[2];
rz(3.0230076) q[3];
sx q[3];
rz(-1.010681) q[3];
sx q[3];
rz(2.3732869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80657715) q[0];
sx q[0];
rz(-1.4310358) q[0];
sx q[0];
rz(3.0203876) q[0];
rz(0.14713261) q[1];
sx q[1];
rz(-1.1558665) q[1];
sx q[1];
rz(-0.8228696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38050851) q[0];
sx q[0];
rz(-0.26603577) q[0];
sx q[0];
rz(2.7903275) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.046209411) q[2];
sx q[2];
rz(-2.5760057) q[2];
sx q[2];
rz(-2.2639745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40892021) q[1];
sx q[1];
rz(-2.3358445) q[1];
sx q[1];
rz(-0.21269704) q[1];
rz(-0.081667451) q[3];
sx q[3];
rz(-1.3586391) q[3];
sx q[3];
rz(-1.3426677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0269028) q[2];
sx q[2];
rz(-2.5451886) q[2];
sx q[2];
rz(-0.58843311) q[2];
rz(2.9450997) q[3];
sx q[3];
rz(-1.7076219) q[3];
sx q[3];
rz(2.4832723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999417) q[0];
sx q[0];
rz(-0.64084941) q[0];
sx q[0];
rz(0.95630056) q[0];
rz(-2.941533) q[1];
sx q[1];
rz(-1.450489) q[1];
sx q[1];
rz(0.53093451) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0414687) q[0];
sx q[0];
rz(-1.6371797) q[0];
sx q[0];
rz(1.081335) q[0];
rz(-pi) q[1];
rz(-2.5052257) q[2];
sx q[2];
rz(-1.8312935) q[2];
sx q[2];
rz(0.72450984) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8968205) q[1];
sx q[1];
rz(-1.4076809) q[1];
sx q[1];
rz(-2.5965124) q[1];
rz(-pi) q[2];
rz(-0.28778971) q[3];
sx q[3];
rz(-0.46817259) q[3];
sx q[3];
rz(-2.5361907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7134573) q[2];
sx q[2];
rz(-2.0811446) q[2];
sx q[2];
rz(1.1327845) q[2];
rz(0.97635859) q[3];
sx q[3];
rz(-2.4553757) q[3];
sx q[3];
rz(1.8648225) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3499632) q[0];
sx q[0];
rz(-3.0296453) q[0];
sx q[0];
rz(1.5891225) q[0];
rz(-2.7516229) q[1];
sx q[1];
rz(-1.5682805) q[1];
sx q[1];
rz(-0.76036298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92780471) q[0];
sx q[0];
rz(-0.47190168) q[0];
sx q[0];
rz(2.9899421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5422212) q[2];
sx q[2];
rz(-1.1090389) q[2];
sx q[2];
rz(0.28959488) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2901569) q[1];
sx q[1];
rz(-2.7692215) q[1];
sx q[1];
rz(1.6769767) q[1];
rz(2.0605956) q[3];
sx q[3];
rz(-2.5464129) q[3];
sx q[3];
rz(-0.90200952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4752263) q[2];
sx q[2];
rz(-2.1658289) q[2];
sx q[2];
rz(2.547612) q[2];
rz(2.5178759) q[3];
sx q[3];
rz(-1.8513605) q[3];
sx q[3];
rz(2.311603) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8754616) q[0];
sx q[0];
rz(-3.0103191) q[0];
sx q[0];
rz(-1.3651146) q[0];
rz(-0.020137067) q[1];
sx q[1];
rz(-2.8485841) q[1];
sx q[1];
rz(1.5905406) q[1];
rz(1.5503442) q[2];
sx q[2];
rz(-0.3849138) q[2];
sx q[2];
rz(0.32017564) q[2];
rz(0.094811335) q[3];
sx q[3];
rz(-2.9651793) q[3];
sx q[3];
rz(1.7549979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
