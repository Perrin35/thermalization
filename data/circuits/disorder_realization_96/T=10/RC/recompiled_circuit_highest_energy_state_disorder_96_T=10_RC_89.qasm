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
rz(-2.8120256) q[0];
sx q[0];
rz(-0.94533935) q[0];
sx q[0];
rz(3.1402631) q[0];
rz(0.20345649) q[1];
sx q[1];
rz(-0.23509547) q[1];
sx q[1];
rz(0.57125616) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093319915) q[0];
sx q[0];
rz(-1.1128582) q[0];
sx q[0];
rz(0.67955525) q[0];
rz(1.4185888) q[2];
sx q[2];
rz(-1.7585351) q[2];
sx q[2];
rz(1.3433045) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9115736) q[1];
sx q[1];
rz(-1.0805478) q[1];
sx q[1];
rz(1.5515627) q[1];
x q[2];
rz(-0.68526973) q[3];
sx q[3];
rz(-0.56495404) q[3];
sx q[3];
rz(0.23811114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0433537) q[2];
sx q[2];
rz(-1.2231772) q[2];
sx q[2];
rz(-2.5377048) q[2];
rz(2.1308925) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(-2.219626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589703) q[0];
sx q[0];
rz(-1.0930467) q[0];
sx q[0];
rz(0.0086722886) q[0];
rz(-1.1588833) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(0.11223665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771292) q[0];
sx q[0];
rz(-0.57687974) q[0];
sx q[0];
rz(-2.5974899) q[0];
rz(1.5331763) q[2];
sx q[2];
rz(-0.92447399) q[2];
sx q[2];
rz(-2.1802156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6735001) q[1];
sx q[1];
rz(-1.348494) q[1];
sx q[1];
rz(-0.50428524) q[1];
rz(0.26437624) q[3];
sx q[3];
rz(-1.8473139) q[3];
sx q[3];
rz(0.087527601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9021641) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(2.2018137) q[2];
rz(-2.2043665) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(-0.37041131) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9409907) q[0];
sx q[0];
rz(-2.2414099) q[0];
sx q[0];
rz(2.5584333) q[0];
rz(1.0519625) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(-3.1254056) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72618851) q[0];
sx q[0];
rz(-2.4016018) q[0];
sx q[0];
rz(-1.9053024) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3397555) q[2];
sx q[2];
rz(-0.41473284) q[2];
sx q[2];
rz(-1.557359) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5213373) q[1];
sx q[1];
rz(-2.9456821) q[1];
sx q[1];
rz(2.1571062) q[1];
rz(-0.82456581) q[3];
sx q[3];
rz(-2.0494132) q[3];
sx q[3];
rz(-0.625911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(-0.072546093) q[2];
rz(1.0091311) q[3];
sx q[3];
rz(-2.3435209) q[3];
sx q[3];
rz(0.23536853) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9811454) q[0];
sx q[0];
rz(-0.35780847) q[0];
sx q[0];
rz(2.2571046) q[0];
rz(-1.5011939) q[1];
sx q[1];
rz(-1.4720935) q[1];
sx q[1];
rz(0.59008682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4761596) q[0];
sx q[0];
rz(-0.64393015) q[0];
sx q[0];
rz(-2.2225999) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8832189) q[2];
sx q[2];
rz(-0.35053634) q[2];
sx q[2];
rz(-2.3141298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1780488) q[1];
sx q[1];
rz(-0.60979382) q[1];
sx q[1];
rz(-2.1765373) q[1];
rz(-pi) q[2];
rz(-1.45118) q[3];
sx q[3];
rz(-1.1362459) q[3];
sx q[3];
rz(1.0512607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13714743) q[2];
sx q[2];
rz(-1.659617) q[2];
sx q[2];
rz(-1.8972634) q[2];
rz(1.7826049) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(2.9984503) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2537848) q[0];
sx q[0];
rz(-2.4847327) q[0];
sx q[0];
rz(-2.7746871) q[0];
rz(-1.7985571) q[1];
sx q[1];
rz(-2.8470706) q[1];
sx q[1];
rz(1.3142745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524805) q[0];
sx q[0];
rz(-2.3398031) q[0];
sx q[0];
rz(-2.5125971) q[0];
rz(-pi) q[1];
rz(2.6268705) q[2];
sx q[2];
rz(-0.54385563) q[2];
sx q[2];
rz(-2.0946696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1340577) q[1];
sx q[1];
rz(-2.4278054) q[1];
sx q[1];
rz(-2.0573045) q[1];
rz(0.35700123) q[3];
sx q[3];
rz(-1.9244266) q[3];
sx q[3];
rz(-2.2382527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.542995) q[2];
sx q[2];
rz(-0.44671217) q[2];
sx q[2];
rz(-2.561595) q[2];
rz(-2.9567806) q[3];
sx q[3];
rz(-1.5671268) q[3];
sx q[3];
rz(-2.4549761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0762894) q[0];
sx q[0];
rz(-0.46173254) q[0];
sx q[0];
rz(-2.8522458) q[0];
rz(-0.099960001) q[1];
sx q[1];
rz(-2.3048765) q[1];
sx q[1];
rz(-2.3438556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983618) q[0];
sx q[0];
rz(-1.7055178) q[0];
sx q[0];
rz(0.93905296) q[0];
rz(-0.50578588) q[2];
sx q[2];
rz(-0.2295851) q[2];
sx q[2];
rz(-1.9167655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3926852) q[1];
sx q[1];
rz(-1.2756383) q[1];
sx q[1];
rz(1.5172076) q[1];
rz(0.84074214) q[3];
sx q[3];
rz(-1.4952097) q[3];
sx q[3];
rz(-0.10160916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37819698) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(2.3217616) q[2];
rz(-1.4929006) q[3];
sx q[3];
rz(-0.35455743) q[3];
sx q[3];
rz(-2.7259887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019526871) q[0];
sx q[0];
rz(-0.30240914) q[0];
sx q[0];
rz(-0.31461883) q[0];
rz(-2.4954691) q[1];
sx q[1];
rz(-1.6982634) q[1];
sx q[1];
rz(-3.019928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3356161) q[0];
sx q[0];
rz(-1.8467963) q[0];
sx q[0];
rz(1.1742623) q[0];
rz(-0.88603772) q[2];
sx q[2];
rz(-2.2412063) q[2];
sx q[2];
rz(1.836402) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77326425) q[1];
sx q[1];
rz(-0.46302893) q[1];
sx q[1];
rz(-2.8431588) q[1];
x q[2];
rz(1.2950778) q[3];
sx q[3];
rz(-2.5150617) q[3];
sx q[3];
rz(-2.3008418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4612215) q[2];
sx q[2];
rz(-1.17522) q[2];
sx q[2];
rz(-2.1628974) q[2];
rz(0.77224246) q[3];
sx q[3];
rz(-2.6148836) q[3];
sx q[3];
rz(1.0129119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0270281) q[0];
sx q[0];
rz(-1.1440729) q[0];
sx q[0];
rz(1.8373328) q[0];
rz(0.14106855) q[1];
sx q[1];
rz(-2.8619659) q[1];
sx q[1];
rz(-1.28349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.210189) q[0];
sx q[0];
rz(-2.922466) q[0];
sx q[0];
rz(1.4019843) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0436921) q[2];
sx q[2];
rz(-1.86964) q[2];
sx q[2];
rz(-3.1063307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5775958) q[1];
sx q[1];
rz(-1.3118) q[1];
sx q[1];
rz(1.6902655) q[1];
rz(-0.047824511) q[3];
sx q[3];
rz(-2.3615814) q[3];
sx q[3];
rz(0.081950233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2804602) q[2];
sx q[2];
rz(-0.91117636) q[2];
sx q[2];
rz(-0.99315161) q[2];
rz(1.447621) q[3];
sx q[3];
rz(-1.9341035) q[3];
sx q[3];
rz(1.8705468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.3250378) q[0];
sx q[0];
rz(-3.0386381) q[0];
sx q[0];
rz(-0.8763985) q[0];
rz(1.5986298) q[1];
sx q[1];
rz(-1.8537268) q[1];
sx q[1];
rz(1.1184982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3680822) q[0];
sx q[0];
rz(-0.72282571) q[0];
sx q[0];
rz(2.5834492) q[0];
rz(2.1886979) q[2];
sx q[2];
rz(-0.94091533) q[2];
sx q[2];
rz(0.057531051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7174445) q[1];
sx q[1];
rz(-2.433648) q[1];
sx q[1];
rz(1.397754) q[1];
rz(-pi) q[2];
rz(-2.4966024) q[3];
sx q[3];
rz(-1.3324311) q[3];
sx q[3];
rz(1.1568943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2529651) q[2];
sx q[2];
rz(-0.94433633) q[2];
sx q[2];
rz(1.2523119) q[2];
rz(2.0400932) q[3];
sx q[3];
rz(-1.7606198) q[3];
sx q[3];
rz(-1.291905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0970704) q[0];
sx q[0];
rz(-3.0945393) q[0];
sx q[0];
rz(-2.4142921) q[0];
rz(0.76668382) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(-1.9511706) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.577548) q[0];
sx q[0];
rz(-2.2909031) q[0];
sx q[0];
rz(-1.4029361) q[0];
x q[1];
rz(0.095049871) q[2];
sx q[2];
rz(-2.2288786) q[2];
sx q[2];
rz(-2.6320145) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0990844) q[1];
sx q[1];
rz(-2.7471424) q[1];
sx q[1];
rz(2.45155) q[1];
rz(-pi) q[2];
rz(-2.6383357) q[3];
sx q[3];
rz(-1.0531593) q[3];
sx q[3];
rz(0.6235756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7850354) q[2];
sx q[2];
rz(-1.0871474) q[2];
sx q[2];
rz(1.1615151) q[2];
rz(-1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(-2.13818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7575191) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(-1.9772498) q[2];
sx q[2];
rz(-2.5377574) q[2];
sx q[2];
rz(-0.24518235) q[2];
rz(1.499621) q[3];
sx q[3];
rz(-2.945032) q[3];
sx q[3];
rz(-1.3155027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
