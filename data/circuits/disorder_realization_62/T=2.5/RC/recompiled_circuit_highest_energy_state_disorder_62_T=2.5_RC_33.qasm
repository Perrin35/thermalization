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
rz(-0.38464889) q[0];
sx q[0];
rz(-2.3031213) q[0];
sx q[0];
rz(0.60929259) q[0];
rz(-2.4203909) q[1];
sx q[1];
rz(-2.2130122) q[1];
sx q[1];
rz(-1.0963) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168092) q[0];
sx q[0];
rz(-1.4583298) q[0];
sx q[0];
rz(2.2676668) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2112261) q[2];
sx q[2];
rz(-1.1784679) q[2];
sx q[2];
rz(2.0517672) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23725736) q[1];
sx q[1];
rz(-1.0544027) q[1];
sx q[1];
rz(-1.2468491) q[1];
x q[2];
rz(1.6814548) q[3];
sx q[3];
rz(-1.7250463) q[3];
sx q[3];
rz(-1.8024886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91813749) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(1.9523331) q[2];
rz(-2.7506645) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(-1.1322359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696058) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(2.4300785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92828099) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(-0.9071945) q[0];
rz(-1.8078126) q[2];
sx q[2];
rz(-2.0227602) q[2];
sx q[2];
rz(2.1353561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7603858) q[1];
sx q[1];
rz(-0.35522705) q[1];
sx q[1];
rz(2.9095501) q[1];
rz(-0.031458843) q[3];
sx q[3];
rz(-2.194187) q[3];
sx q[3];
rz(-1.7707847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4352162) q[2];
sx q[2];
rz(-1.3560359) q[2];
sx q[2];
rz(-1.0630652) q[2];
rz(1.4937909) q[3];
sx q[3];
rz(-0.46896514) q[3];
sx q[3];
rz(-2.4013605) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73055926) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(-2.539769) q[0];
rz(2.9959294) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(0.62785968) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3224761) q[0];
sx q[0];
rz(-1.5014582) q[0];
sx q[0];
rz(-0.11868166) q[0];
rz(0.35618393) q[2];
sx q[2];
rz(-1.0477475) q[2];
sx q[2];
rz(-1.5580391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8257432) q[1];
sx q[1];
rz(-1.2176732) q[1];
sx q[1];
rz(-0.24314116) q[1];
x q[2];
rz(-2.1842923) q[3];
sx q[3];
rz(-1.9578551) q[3];
sx q[3];
rz(-1.330803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38498983) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(-0.21558726) q[2];
rz(-2.0032517) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(-2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1348006) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(1.6670594) q[1];
sx q[1];
rz(-0.88913616) q[1];
sx q[1];
rz(2.3337591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0197163) q[0];
sx q[0];
rz(-0.84152475) q[0];
sx q[0];
rz(-2.3725933) q[0];
rz(-1.4363244) q[2];
sx q[2];
rz(-2.8155012) q[2];
sx q[2];
rz(1.6539314) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4588759) q[1];
sx q[1];
rz(-1.4035051) q[1];
sx q[1];
rz(0.95750995) q[1];
rz(-pi) q[2];
rz(2.4633154) q[3];
sx q[3];
rz(-1.5828352) q[3];
sx q[3];
rz(1.2348242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5011751) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(-2.625722) q[2];
rz(-3.0473895) q[3];
sx q[3];
rz(-2.3661864) q[3];
sx q[3];
rz(1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.1871673) q[0];
rz(-0.59133235) q[1];
sx q[1];
rz(-1.8449123) q[1];
sx q[1];
rz(2.3147413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7477404) q[0];
sx q[0];
rz(-1.4497787) q[0];
sx q[0];
rz(2.7583102) q[0];
rz(2.128878) q[2];
sx q[2];
rz(-1.5019755) q[2];
sx q[2];
rz(-2.0134913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0650658) q[1];
sx q[1];
rz(-0.69628104) q[1];
sx q[1];
rz(-1.3672921) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86381819) q[3];
sx q[3];
rz(-0.70138273) q[3];
sx q[3];
rz(2.6849298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(2.1264326) q[2];
rz(-2.2319131) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(2.4801586) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(0.30995187) q[0];
rz(-0.018772086) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(0.087336691) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7531484) q[0];
sx q[0];
rz(-0.57373057) q[0];
sx q[0];
rz(3.0818207) q[0];
x q[1];
rz(0.87126731) q[2];
sx q[2];
rz(-1.0443029) q[2];
sx q[2];
rz(0.73824182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2948871) q[1];
sx q[1];
rz(-2.069983) q[1];
sx q[1];
rz(0.45392068) q[1];
x q[2];
rz(2.099192) q[3];
sx q[3];
rz(-1.434552) q[3];
sx q[3];
rz(-0.50627499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7901018) q[2];
sx q[2];
rz(-0.47372207) q[2];
sx q[2];
rz(-0.52004415) q[2];
rz(-3.0067048) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(2.7556162) q[0];
rz(1.0999673) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(0.78752548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58819492) q[0];
sx q[0];
rz(-0.46632354) q[0];
sx q[0];
rz(0.88769261) q[0];
rz(-pi) q[1];
rz(1.196127) q[2];
sx q[2];
rz(-2.154899) q[2];
sx q[2];
rz(2.0541592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3508151) q[1];
sx q[1];
rz(-2.3190089) q[1];
sx q[1];
rz(2.050553) q[1];
rz(-pi) q[2];
rz(0.16997108) q[3];
sx q[3];
rz(-1.8026226) q[3];
sx q[3];
rz(-1.9915723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1197352) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(1.5137399) q[2];
rz(-1.2210023) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.8319703) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(-2.9085462) q[0];
rz(1.6538992) q[1];
sx q[1];
rz(-2.1079886) q[1];
sx q[1];
rz(2.1344562) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876935) q[0];
sx q[0];
rz(-1.0570044) q[0];
sx q[0];
rz(-1.9625241) q[0];
rz(-0.39733823) q[2];
sx q[2];
rz(-2.4724744) q[2];
sx q[2];
rz(2.4934409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5598232) q[1];
sx q[1];
rz(-1.9672462) q[1];
sx q[1];
rz(2.9828986) q[1];
rz(3.0606224) q[3];
sx q[3];
rz(-0.24484466) q[3];
sx q[3];
rz(3.0449344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2795589) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(-0.88279185) q[2];
rz(-0.27859303) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(0.45894233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19926628) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(-1.5698154) q[0];
rz(-2.3528631) q[1];
sx q[1];
rz(-0.9659583) q[1];
sx q[1];
rz(0.68005651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71841371) q[0];
sx q[0];
rz(-0.81162757) q[0];
sx q[0];
rz(1.4885097) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81794195) q[2];
sx q[2];
rz(-2.9688058) q[2];
sx q[2];
rz(0.30889749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9100807) q[1];
sx q[1];
rz(-0.68889131) q[1];
sx q[1];
rz(-2.5206294) q[1];
rz(0.84362883) q[3];
sx q[3];
rz(-1.1427726) q[3];
sx q[3];
rz(0.53800636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.9752768) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(0.62674633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.523943) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(-2.7045265) q[0];
rz(0.79824671) q[1];
sx q[1];
rz(-2.6415446) q[1];
sx q[1];
rz(0.22013586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58071346) q[0];
sx q[0];
rz(-2.3395943) q[0];
sx q[0];
rz(0.064994354) q[0];
rz(1.2557477) q[2];
sx q[2];
rz(-1.7654668) q[2];
sx q[2];
rz(2.6597629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0221631) q[1];
sx q[1];
rz(-0.8406175) q[1];
sx q[1];
rz(0.37067159) q[1];
rz(-pi) q[2];
rz(2.3399971) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(-0.012594087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2685214) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.2344454) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81454043) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(3.0771599) q[1];
sx q[1];
rz(-1.671052) q[1];
sx q[1];
rz(-0.64238092) q[1];
rz(2.1656373) q[2];
sx q[2];
rz(-1.5968762) q[2];
sx q[2];
rz(2.1322875) q[2];
rz(0.53917428) q[3];
sx q[3];
rz(-1.3414345) q[3];
sx q[3];
rz(0.10673005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
