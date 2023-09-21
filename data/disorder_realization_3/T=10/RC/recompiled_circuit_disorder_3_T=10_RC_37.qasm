OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(-2.3303633) q[0];
rz(1.7180874) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(-1.5073843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1297076) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(3.0711864) q[1];
rz(-0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(-0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-2.6020715) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8115494) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(0.051785843) q[0];
rz(-0.72210724) q[2];
sx q[2];
rz(-1.2566084) q[2];
sx q[2];
rz(-2.9177641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9425548) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(1.1953137) q[1];
x q[2];
rz(2.1971365) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-2.8163547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494333) q[0];
sx q[0];
rz(-1.8349577) q[0];
sx q[0];
rz(2.507472) q[0];
rz(0.44552866) q[2];
sx q[2];
rz(-1.1330714) q[2];
sx q[2];
rz(0.23756269) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8032854) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(2.0558946) q[1];
rz(0.20817169) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8024575) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(-0.25948778) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(2.1156103) q[0];
x q[1];
rz(-0.133693) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(-2.0267817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(2.494032) q[1];
rz(2.6287574) q[3];
sx q[3];
rz(-0.25179112) q[3];
sx q[3];
rz(2.4076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(1.4506838) q[0];
rz(-2.5227929) q[2];
sx q[2];
rz(-2.7745562) q[2];
sx q[2];
rz(2.7235081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0610173) q[1];
sx q[1];
rz(-1.4772381) q[1];
sx q[1];
rz(-2.1993125) q[1];
rz(-pi) q[2];
rz(-0.4425211) q[3];
sx q[3];
rz(-1.8707152) q[3];
sx q[3];
rz(-0.72344852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(-1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87175831) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(1.7437177) q[0];
rz(0.78380615) q[2];
sx q[2];
rz(-1.1597826) q[2];
sx q[2];
rz(2.59936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1706714) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(1.7621653) q[1];
rz(-pi) q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(3.133657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7769653) q[0];
sx q[0];
rz(-2.1205175) q[0];
sx q[0];
rz(1.8718029) q[0];
rz(-pi) q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.9108859) q[2];
sx q[2];
rz(-2.0932587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24212813) q[1];
sx q[1];
rz(-2.2367396) q[1];
sx q[1];
rz(-0.66023402) q[1];
rz(1.1778529) q[3];
sx q[3];
rz(-1.279139) q[3];
sx q[3];
rz(1.1100811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68740326) q[0];
sx q[0];
rz(-2.8912376) q[0];
sx q[0];
rz(-0.040391163) q[0];
rz(1.0566696) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(2.4664997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2316206) q[1];
sx q[1];
rz(-0.93712229) q[1];
sx q[1];
rz(-0.85068591) q[1];
x q[2];
rz(1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(-0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(2.705943) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1736261) q[0];
sx q[0];
rz(-2.2244503) q[0];
sx q[0];
rz(3.0941512) q[0];
rz(-0.68848227) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(-2.288523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7855362) q[1];
sx q[1];
rz(-2.0691263) q[1];
sx q[1];
rz(-0.15798012) q[1];
x q[2];
rz(1.7008408) q[3];
sx q[3];
rz(-1.4308883) q[3];
sx q[3];
rz(-3.1062982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(2.9157675) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(-0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.8189925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117301) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25090353) q[2];
sx q[2];
rz(-1.13399) q[2];
sx q[2];
rz(-1.3383588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1111856) q[1];
sx q[1];
rz(-2.2135995) q[1];
sx q[1];
rz(-0.91099693) q[1];
x q[2];
rz(1.7694468) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(-2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.50921847) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
