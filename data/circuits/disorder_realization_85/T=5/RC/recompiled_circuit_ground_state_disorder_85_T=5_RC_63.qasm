OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10211927) q[0];
sx q[0];
rz(-1.6685628) q[0];
sx q[0];
rz(3.0102475) q[0];
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(1.8324469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0378758) q[0];
sx q[0];
rz(-1.193422) q[0];
sx q[0];
rz(2.3233633) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88850682) q[2];
sx q[2];
rz(-0.66239385) q[2];
sx q[2];
rz(-2.2531525) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.943191) q[1];
sx q[1];
rz(-1.3326338) q[1];
sx q[1];
rz(-2.8203301) q[1];
x q[2];
rz(-0.18351002) q[3];
sx q[3];
rz(-1.402463) q[3];
sx q[3];
rz(0.77372293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(2.479539) q[2];
rz(-0.74696294) q[3];
sx q[3];
rz(-0.91619879) q[3];
sx q[3];
rz(1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58716431) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(-1.9594877) q[0];
rz(-0.74554044) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(-0.10339698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1348159) q[0];
sx q[0];
rz(-2.6511483) q[0];
sx q[0];
rz(-0.51014002) q[0];
rz(-pi) q[1];
x q[1];
rz(2.75704) q[2];
sx q[2];
rz(-2.1332624) q[2];
sx q[2];
rz(0.99316521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34161257) q[1];
sx q[1];
rz(-2.3368521) q[1];
sx q[1];
rz(-2.0384203) q[1];
x q[2];
rz(2.313226) q[3];
sx q[3];
rz(-1.6992555) q[3];
sx q[3];
rz(2.1118233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5891002) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(-2.6972771) q[2];
rz(2.0878504) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0788197) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-2.479082) q[0];
rz(-1.484681) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(0.2167162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22610006) q[0];
sx q[0];
rz(-2.2558172) q[0];
sx q[0];
rz(2.1846376) q[0];
x q[1];
rz(-2.4091085) q[2];
sx q[2];
rz(-1.2353131) q[2];
sx q[2];
rz(-1.6215289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89284929) q[1];
sx q[1];
rz(-2.4676552) q[1];
sx q[1];
rz(-1.5065864) q[1];
x q[2];
rz(-0.29189887) q[3];
sx q[3];
rz(-0.51288285) q[3];
sx q[3];
rz(1.8512902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66389877) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(-0.23507512) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821871) q[0];
sx q[0];
rz(-3.0538054) q[0];
sx q[0];
rz(-2.1413595) q[0];
rz(-1.9384711) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(-2.5968754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31123589) q[0];
sx q[0];
rz(-0.85842997) q[0];
sx q[0];
rz(-1.6300549) q[0];
rz(2.4248872) q[2];
sx q[2];
rz(-0.96895987) q[2];
sx q[2];
rz(2.3284137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5627553) q[1];
sx q[1];
rz(-1.3816427) q[1];
sx q[1];
rz(-0.12054969) q[1];
x q[2];
rz(-0.29981837) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(-0.064379582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41877052) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(-2.344632) q[2];
rz(0.289251) q[3];
sx q[3];
rz(-0.97024337) q[3];
sx q[3];
rz(0.42118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61409426) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.2689137) q[0];
rz(0.96619636) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(2.6752245) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9556676) q[0];
sx q[0];
rz(-2.4644682) q[0];
sx q[0];
rz(-2.4311275) q[0];
rz(-pi) q[1];
rz(-0.63299243) q[2];
sx q[2];
rz(-2.5883753) q[2];
sx q[2];
rz(-0.67942441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6154229) q[1];
sx q[1];
rz(-0.96994441) q[1];
sx q[1];
rz(0.38819617) q[1];
rz(-pi) q[2];
rz(2.0556695) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7352778) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(0.80424133) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-1.0420957) q[3];
sx q[3];
rz(3.13412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2020579) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(0.95788389) q[0];
rz(1.6288039) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(1.5527976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2561172) q[0];
sx q[0];
rz(-1.4909571) q[0];
sx q[0];
rz(0.022932963) q[0];
x q[1];
rz(0.27811173) q[2];
sx q[2];
rz(-1.2956923) q[2];
sx q[2];
rz(0.29488568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54974906) q[1];
sx q[1];
rz(-1.7477682) q[1];
sx q[1];
rz(-2.4787477) q[1];
x q[2];
rz(-0.53258606) q[3];
sx q[3];
rz(-0.69952337) q[3];
sx q[3];
rz(-0.12040779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-2.4064257) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(2.2775876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76148024) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(0.6066221) q[0];
rz(0.78423777) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(-1.7971136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356268) q[0];
sx q[0];
rz(-0.84442645) q[0];
sx q[0];
rz(-0.65610049) q[0];
rz(0.44353087) q[2];
sx q[2];
rz(-2.6012528) q[2];
sx q[2];
rz(-2.4176554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91025464) q[1];
sx q[1];
rz(-1.7926072) q[1];
sx q[1];
rz(-2.7342058) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3913888) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(-0.52667945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6091696) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(-2.9621647) q[3];
sx q[3];
rz(-2.1556985) q[3];
sx q[3];
rz(0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0824025) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(2.3865336) q[0];
rz(0.45626196) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(2.0645352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74769831) q[0];
sx q[0];
rz(-1.0193045) q[0];
sx q[0];
rz(-1.8389788) q[0];
rz(-1.6918534) q[2];
sx q[2];
rz(-0.79059764) q[2];
sx q[2];
rz(1.7771174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7491319) q[1];
sx q[1];
rz(-1.2866255) q[1];
sx q[1];
rz(0.47537132) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8313238) q[3];
sx q[3];
rz(-2.5857877) q[3];
sx q[3];
rz(-0.99930857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(-0.73823482) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.8176707) q[3];
sx q[3];
rz(1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(2.451918) q[1];
sx q[1];
rz(-2.4318305) q[1];
sx q[1];
rz(-2.7241657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8807424) q[0];
sx q[0];
rz(-2.1584903) q[0];
sx q[0];
rz(1.8729314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8418188) q[2];
sx q[2];
rz(-1.6437141) q[2];
sx q[2];
rz(2.656183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.058674) q[1];
sx q[1];
rz(-1.6750458) q[1];
sx q[1];
rz(0.41091316) q[1];
x q[2];
rz(-1.1754009) q[3];
sx q[3];
rz(-1.1987276) q[3];
sx q[3];
rz(-1.7898066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4735585) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(3.0958946) q[2];
rz(-2.8081196) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(-1.217655) q[0];
rz(0.24670163) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(0.47952476) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3772904) q[0];
sx q[0];
rz(-2.3545697) q[0];
sx q[0];
rz(2.2020409) q[0];
rz(0.87534753) q[2];
sx q[2];
rz(-1.1743059) q[2];
sx q[2];
rz(-0.04405313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3563167) q[1];
sx q[1];
rz(-2.0802167) q[1];
sx q[1];
rz(1.9059577) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3334951) q[3];
sx q[3];
rz(-1.1166683) q[3];
sx q[3];
rz(-0.25179201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5852927) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(1.7362107) q[2];
rz(-0.65226883) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(-2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(2.3975092) q[1];
sx q[1];
rz(-0.71744812) q[1];
sx q[1];
rz(1.7070028) q[1];
rz(-1.4135398) q[2];
sx q[2];
rz(-0.99307151) q[2];
sx q[2];
rz(2.0057783) q[2];
rz(0.35318315) q[3];
sx q[3];
rz(-1.0241057) q[3];
sx q[3];
rz(3.1288341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
