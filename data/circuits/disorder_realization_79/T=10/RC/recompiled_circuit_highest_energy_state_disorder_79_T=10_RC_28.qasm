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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(2.9377687) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(1.6066983) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8730072) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(-2.068822) q[0];
rz(-pi) q[1];
rz(0.34687545) q[2];
sx q[2];
rz(-2.2636608) q[2];
sx q[2];
rz(1.833806) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.038329934) q[1];
sx q[1];
rz(-0.94417773) q[1];
sx q[1];
rz(2.7019397) q[1];
rz(-pi) q[2];
rz(-1.7059709) q[3];
sx q[3];
rz(-1.7723284) q[3];
sx q[3];
rz(-2.7592748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(0.99754769) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7928829) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(1.8541699) q[0];
rz(-2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.2189254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588803) q[0];
sx q[0];
rz(-0.25280372) q[0];
sx q[0];
rz(-0.88835277) q[0];
x q[1];
rz(-2.4697815) q[2];
sx q[2];
rz(-1.3992568) q[2];
sx q[2];
rz(1.951527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7658577) q[1];
sx q[1];
rz(-1.2925783) q[1];
sx q[1];
rz(0.15032676) q[1];
x q[2];
rz(0.35393012) q[3];
sx q[3];
rz(-2.2316389) q[3];
sx q[3];
rz(-2.5315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37811849) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(0.62180579) q[2];
rz(0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(-2.9252692) q[0];
rz(-0.59858876) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(2.6187706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7824088) q[0];
sx q[0];
rz(-2.6893927) q[0];
sx q[0];
rz(-0.81069273) q[0];
rz(-0.36575138) q[2];
sx q[2];
rz(-0.41588441) q[2];
sx q[2];
rz(-1.5378086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93980689) q[1];
sx q[1];
rz(-0.56957376) q[1];
sx q[1];
rz(2.7199183) q[1];
rz(-0.82562311) q[3];
sx q[3];
rz(-2.7047727) q[3];
sx q[3];
rz(-0.88479155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(-1.5175021) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9972123) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(-1.602518) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(-1.296952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3281477) q[0];
sx q[0];
rz(-1.7792276) q[0];
sx q[0];
rz(0.84020241) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0791666) q[2];
sx q[2];
rz(-0.87597825) q[2];
sx q[2];
rz(2.7165987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8140111) q[1];
sx q[1];
rz(-1.9821229) q[1];
sx q[1];
rz(2.7279305) q[1];
rz(-pi) q[2];
rz(-1.2554712) q[3];
sx q[3];
rz(-2.1026169) q[3];
sx q[3];
rz(1.1602311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(-1.0456592) q[2];
rz(-2.9271434) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8056718) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(0.51934284) q[0];
rz(-2.1995811) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(1.5325783) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518884) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(-0.69464442) q[0];
x q[1];
rz(-2.9608594) q[2];
sx q[2];
rz(-1.7109979) q[2];
sx q[2];
rz(0.069958036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6841507) q[1];
sx q[1];
rz(-1.5655715) q[1];
sx q[1];
rz(1.0702151) q[1];
rz(0.71096731) q[3];
sx q[3];
rz(-2.9971854) q[3];
sx q[3];
rz(2.292423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64451009) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(-2.8572148) q[2];
rz(-2.583368) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(2.5010342) q[0];
rz(1.5156281) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(2.9277149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77737264) q[0];
sx q[0];
rz(-1.7988482) q[0];
sx q[0];
rz(1.26012) q[0];
rz(-2.7340545) q[2];
sx q[2];
rz(-0.78754163) q[2];
sx q[2];
rz(2.119273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0391024) q[1];
sx q[1];
rz(-0.63345816) q[1];
sx q[1];
rz(-1.6879115) q[1];
rz(-pi) q[2];
rz(2.0908666) q[3];
sx q[3];
rz(-1.5087391) q[3];
sx q[3];
rz(-2.7620535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(0.17769979) q[2];
rz(-1.7443582) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419256) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(-1.1055111) q[0];
rz(-2.8583177) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(1.4615321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4257129) q[0];
sx q[0];
rz(-2.113803) q[0];
sx q[0];
rz(-2.3433861) q[0];
rz(-pi) q[1];
rz(-0.92723691) q[2];
sx q[2];
rz(-2.3935648) q[2];
sx q[2];
rz(-3.1035556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6532306) q[1];
sx q[1];
rz(-1.686859) q[1];
sx q[1];
rz(-2.577072) q[1];
x q[2];
rz(-2.5652065) q[3];
sx q[3];
rz(-1.295305) q[3];
sx q[3];
rz(-0.34289962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(2.936787) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(-1.9974476) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56267101) q[0];
sx q[0];
rz(-1.5961338) q[0];
sx q[0];
rz(0.48669215) q[0];
rz(-1.2828243) q[2];
sx q[2];
rz(-0.33703732) q[2];
sx q[2];
rz(2.2971643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28539666) q[1];
sx q[1];
rz(-0.90803972) q[1];
sx q[1];
rz(2.0539356) q[1];
rz(-pi) q[2];
rz(-1.2253292) q[3];
sx q[3];
rz(-1.8356712) q[3];
sx q[3];
rz(-0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-2.3365946) q[2];
sx q[2];
rz(-1.8899274) q[2];
rz(-1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(-0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.0018472483) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(1.2563323) q[0];
rz(-3.046335) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(2.6108066) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7408001) q[0];
sx q[0];
rz(-0.2966899) q[0];
sx q[0];
rz(0.99389561) q[0];
rz(-pi) q[1];
rz(1.6094535) q[2];
sx q[2];
rz(-1.8886023) q[2];
sx q[2];
rz(-3.0927049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8618882) q[1];
sx q[1];
rz(-1.8100396) q[1];
sx q[1];
rz(-1.7025392) q[1];
rz(-pi) q[2];
rz(2.6952031) q[3];
sx q[3];
rz(-1.4671031) q[3];
sx q[3];
rz(0.54561347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8119767) q[2];
sx q[2];
rz(-1.5584385) q[2];
sx q[2];
rz(2.6007593) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(0.3987819) q[0];
rz(-1.3840236) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(3.0922281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3545586) q[0];
sx q[0];
rz(-1.462121) q[0];
sx q[0];
rz(-3.0804481) q[0];
rz(-1.1902383) q[2];
sx q[2];
rz(-0.28186381) q[2];
sx q[2];
rz(-2.6388002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80222622) q[1];
sx q[1];
rz(-0.58103937) q[1];
sx q[1];
rz(-1.4895579) q[1];
x q[2];
rz(1.7743763) q[3];
sx q[3];
rz(-2.1927811) q[3];
sx q[3];
rz(2.3443215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(1.9063037) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8265726) q[0];
sx q[0];
rz(-0.44354225) q[0];
sx q[0];
rz(-1.0916239) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(0.51382463) q[2];
sx q[2];
rz(-1.1305625) q[2];
sx q[2];
rz(2.5827575) q[2];
rz(2.1178653) q[3];
sx q[3];
rz(-1.8454875) q[3];
sx q[3];
rz(-1.9133205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
