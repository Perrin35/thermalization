OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(0.13134512) q[0];
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(1.8324469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86517559) q[0];
sx q[0];
rz(-0.88224471) q[0];
sx q[0];
rz(0.49746969) q[0];
x q[1];
rz(-2.2530858) q[2];
sx q[2];
rz(-0.66239385) q[2];
sx q[2];
rz(2.2531525) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6908474) q[1];
sx q[1];
rz(-1.8826797) q[1];
sx q[1];
rz(-1.8212832) q[1];
x q[2];
rz(1.3996436) q[3];
sx q[3];
rz(-1.7516836) q[3];
sx q[3];
rz(2.3756042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(-2.479539) q[2];
rz(2.3946297) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(2.1257187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58716431) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(1.9594877) q[0];
rz(2.3960522) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(3.0381957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4415084) q[0];
sx q[0];
rz(-1.9943976) q[0];
sx q[0];
rz(1.3157428) q[0];
rz(-1.0339884) q[2];
sx q[2];
rz(-0.66945449) q[2];
sx q[2];
rz(2.7975999) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34161257) q[1];
sx q[1];
rz(-0.80474058) q[1];
sx q[1];
rz(-1.1031723) q[1];
x q[2];
rz(-1.7595791) q[3];
sx q[3];
rz(-0.75136614) q[3];
sx q[3];
rz(0.40230678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5891002) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(2.6972771) q[2];
rz(-1.0537423) q[3];
sx q[3];
rz(-3.0709303) q[3];
sx q[3];
rz(-1.9452555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0788197) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-0.66251063) q[0];
rz(1.484681) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(2.9248765) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9154926) q[0];
sx q[0];
rz(-0.88577548) q[0];
sx q[0];
rz(0.95695509) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1323053) q[2];
sx q[2];
rz(-0.88738933) q[2];
sx q[2];
rz(-0.33861288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2487434) q[1];
sx q[1];
rz(-0.67393747) q[1];
sx q[1];
rz(1.6350063) q[1];
x q[2];
rz(2.8496938) q[3];
sx q[3];
rz(-0.51288285) q[3];
sx q[3];
rz(-1.2903024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66389877) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(-0.23507512) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(-0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821871) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-1.0002332) q[0];
rz(1.9384711) q[1];
sx q[1];
rz(-1.6327881) q[1];
sx q[1];
rz(0.54471725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207912) q[0];
sx q[0];
rz(-1.5259597) q[0];
sx q[0];
rz(2.4283571) q[0];
x q[1];
rz(-2.3096309) q[2];
sx q[2];
rz(-2.1429981) q[2];
sx q[2];
rz(-1.9257279) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57883731) q[1];
sx q[1];
rz(-1.75995) q[1];
sx q[1];
rz(3.021043) q[1];
rz(-pi) q[2];
rz(-2.1837213) q[3];
sx q[3];
rz(-1.3231233) q[3];
sx q[3];
rz(1.4644458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41877052) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(-0.79696068) q[2];
rz(-0.289251) q[3];
sx q[3];
rz(-0.97024337) q[3];
sx q[3];
rz(2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61409426) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.8726789) q[0];
rz(-2.1753963) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(-0.46636811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.185925) q[0];
sx q[0];
rz(-2.4644682) q[0];
sx q[0];
rz(2.4311275) q[0];
x q[1];
rz(-2.5086002) q[2];
sx q[2];
rz(-0.55321732) q[2];
sx q[2];
rz(2.4621682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1524017) q[1];
sx q[1];
rz(-2.4394803) q[1];
sx q[1];
rz(-2.0753808) q[1];
rz(-pi) q[2];
rz(-1.0859231) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7352778) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(-2.3373513) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-1.0420957) q[3];
sx q[3];
rz(-0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2020579) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(-2.1837088) q[0];
rz(-1.6288039) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(1.588795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88547546) q[0];
sx q[0];
rz(-1.4909571) q[0];
sx q[0];
rz(3.1186597) q[0];
x q[1];
rz(-2.3423296) q[2];
sx q[2];
rz(-2.7529319) q[2];
sx q[2];
rz(-2.6262019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9839638) q[1];
sx q[1];
rz(-0.92008725) q[1];
sx q[1];
rz(-1.3476861) q[1];
rz(-pi) q[2];
rz(-2.6090066) q[3];
sx q[3];
rz(-2.4420693) q[3];
sx q[3];
rz(-0.12040779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9396886) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(0.73516694) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(2.2775876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801124) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(2.5349706) q[0];
rz(-2.3573549) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(1.7971136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356268) q[0];
sx q[0];
rz(-0.84442645) q[0];
sx q[0];
rz(-0.65610049) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6450301) q[2];
sx q[2];
rz(-1.3482058) q[2];
sx q[2];
rz(2.6816161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1320051) q[1];
sx q[1];
rz(-0.46087206) q[1];
sx q[1];
rz(-0.51746093) q[1];
rz(-1.3913888) q[3];
sx q[3];
rz(-0.47166892) q[3];
sx q[3];
rz(-2.6149132) q[3];
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
rz(2.9621647) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(-2.7348886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0591902) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(0.75505906) q[0];
rz(2.6853307) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(-1.0770575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614968) q[0];
sx q[0];
rz(-1.7984522) q[0];
sx q[0];
rz(0.56786768) q[0];
x q[1];
rz(3.0201689) q[2];
sx q[2];
rz(-2.3540263) q[2];
sx q[2];
rz(1.9483669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4619948) q[1];
sx q[1];
rz(-0.5481998) q[1];
sx q[1];
rz(2.5735709) q[1];
rz(-pi) q[2];
rz(2.1113339) q[3];
sx q[3];
rz(-1.7071299) q[3];
sx q[3];
rz(2.7928074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(-0.73823482) q[2];
rz(-1.632656) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.9807321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974834) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(2.4926376) q[0];
rz(-2.451918) q[1];
sx q[1];
rz(-2.4318305) q[1];
sx q[1];
rz(-0.41742691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605292) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(2.5323243) q[0];
x q[1];
rz(-2.8990977) q[2];
sx q[2];
rz(-2.8333377) q[2];
sx q[2];
rz(0.85390845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.058674) q[1];
sx q[1];
rz(-1.6750458) q[1];
sx q[1];
rz(2.7306795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1754009) q[3];
sx q[3];
rz(-1.1987276) q[3];
sx q[3];
rz(1.7898066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(-3.0958946) q[2];
rz(2.8081196) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(0.12588178) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(1.217655) q[0];
rz(0.24670163) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(-0.47952476) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2829153) q[0];
sx q[0];
rz(-2.0020131) q[0];
sx q[0];
rz(-0.89003508) q[0];
x q[1];
rz(2.2662451) q[2];
sx q[2];
rz(-1.1743059) q[2];
sx q[2];
rz(-3.0975395) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5243234) q[1];
sx q[1];
rz(-1.2795382) q[1];
sx q[1];
rz(-0.53415438) q[1];
rz(0.4654309) q[3];
sx q[3];
rz(-1.7836535) q[3];
sx q[3];
rz(-1.928291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(1.4053819) q[2];
rz(2.4893238) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(-2.4238267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(0.74408342) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-1.7280528) q[2];
sx q[2];
rz(-2.1485211) q[2];
sx q[2];
rz(-1.1358144) q[2];
rz(0.99540972) q[3];
sx q[3];
rz(-1.2708455) q[3];
sx q[3];
rz(1.3686913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
