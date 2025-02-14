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
rz(-3.0120612) q[0];
sx q[0];
rz(-3.01053) q[0];
sx q[0];
rz(0.60453209) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(-0.25783917) q[1];
sx q[1];
rz(-1.2300307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8800921) q[0];
sx q[0];
rz(-1.8207491) q[0];
sx q[0];
rz(0.34323741) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8240439) q[2];
sx q[2];
rz(-2.5032836) q[2];
sx q[2];
rz(0.93250634) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9046272) q[1];
sx q[1];
rz(-0.89348999) q[1];
sx q[1];
rz(-0.64579247) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0139066) q[3];
sx q[3];
rz(-2.0632072) q[3];
sx q[3];
rz(-3.1000053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0872385) q[2];
sx q[2];
rz(-1.9638991) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(-0.44998351) q[3];
sx q[3];
rz(-0.68322244) q[3];
sx q[3];
rz(-2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0403274) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(0.26327565) q[0];
rz(1.0385665) q[1];
sx q[1];
rz(-2.7949605) q[1];
sx q[1];
rz(-0.77478772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446774) q[0];
sx q[0];
rz(-1.7763486) q[0];
sx q[0];
rz(0.62384042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3900142) q[2];
sx q[2];
rz(-0.84779352) q[2];
sx q[2];
rz(2.4931049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0475743) q[1];
sx q[1];
rz(-1.1915922) q[1];
sx q[1];
rz(-3.1017041) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0618013) q[3];
sx q[3];
rz(-0.92872075) q[3];
sx q[3];
rz(1.1428558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7210377) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(-2.5431385) q[2];
rz(-1.7789486) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(0.27073282) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(-1.5077952) q[0];
rz(2.9871121) q[1];
sx q[1];
rz(-1.4198317) q[1];
sx q[1];
rz(-0.78937626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75784439) q[0];
sx q[0];
rz(-2.0998635) q[0];
sx q[0];
rz(1.4704513) q[0];
x q[1];
rz(-1.4285226) q[2];
sx q[2];
rz(-1.8492438) q[2];
sx q[2];
rz(1.5253893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37395135) q[1];
sx q[1];
rz(-2.1726296) q[1];
sx q[1];
rz(-2.7695719) q[1];
rz(-2.0752729) q[3];
sx q[3];
rz(-2.3036727) q[3];
sx q[3];
rz(-1.63597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(-1.6645128) q[2];
rz(-1.9278256) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(1.7267797) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690935) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(-2.5308727) q[0];
rz(-2.4702813) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(-1.3549365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31025654) q[0];
sx q[0];
rz(-0.68643565) q[0];
sx q[0];
rz(-0.75365922) q[0];
rz(-1.2569095) q[2];
sx q[2];
rz(-0.38849026) q[2];
sx q[2];
rz(-1.1780648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6411533) q[1];
sx q[1];
rz(-1.3257003) q[1];
sx q[1];
rz(0.064747253) q[1];
x q[2];
rz(1.4282966) q[3];
sx q[3];
rz(-2.550052) q[3];
sx q[3];
rz(-1.9672036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9352202) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(2.8727403) q[2];
rz(-0.74784652) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(1.6929172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95626962) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(-2.4638033) q[0];
rz(-0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.6544624) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9155898) q[0];
sx q[0];
rz(-1.5446413) q[0];
sx q[0];
rz(1.4348861) q[0];
rz(-0.057737902) q[2];
sx q[2];
rz(-1.1694093) q[2];
sx q[2];
rz(-0.58338469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8890052) q[1];
sx q[1];
rz(-2.41373) q[1];
sx q[1];
rz(0.83830203) q[1];
rz(-3.0934382) q[3];
sx q[3];
rz(-1.9461402) q[3];
sx q[3];
rz(-1.0182235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80374485) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(-2.7860876) q[2];
rz(-2.096094) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(-2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97583714) q[0];
sx q[0];
rz(-2.5582357) q[0];
sx q[0];
rz(3.052886) q[0];
rz(-1.6070131) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-0.075693695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82268084) q[0];
sx q[0];
rz(-1.3782129) q[0];
sx q[0];
rz(2.4680424) q[0];
rz(-pi) q[1];
rz(-3.0906601) q[2];
sx q[2];
rz(-2.3017852) q[2];
sx q[2];
rz(2.5185846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9348889) q[1];
sx q[1];
rz(-0.57414251) q[1];
sx q[1];
rz(-1.8341792) q[1];
x q[2];
rz(-0.57689835) q[3];
sx q[3];
rz(-1.8460974) q[3];
sx q[3];
rz(-2.4979765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82464108) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(1.0392044) q[2];
rz(-0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067326389) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(1.0443895) q[0];
rz(-0.77230612) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-2.4078802) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45424592) q[0];
sx q[0];
rz(-0.92763072) q[0];
sx q[0];
rz(-2.1976794) q[0];
rz(-1.9416642) q[2];
sx q[2];
rz(-1.7108016) q[2];
sx q[2];
rz(-2.3429371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69960475) q[1];
sx q[1];
rz(-1.8612222) q[1];
sx q[1];
rz(1.4493577) q[1];
rz(-pi) q[2];
rz(0.10802631) q[3];
sx q[3];
rz(-0.83893925) q[3];
sx q[3];
rz(-0.15209231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61651984) q[2];
sx q[2];
rz(-2.6313582) q[2];
sx q[2];
rz(-2.5329242) q[2];
rz(-1.0673374) q[3];
sx q[3];
rz(-0.87456861) q[3];
sx q[3];
rz(-0.55060351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50255018) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(-1.6453561) q[0];
rz(-1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(0.38633698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423415) q[0];
sx q[0];
rz(-0.73495293) q[0];
sx q[0];
rz(2.7432261) q[0];
rz(0.069938439) q[2];
sx q[2];
rz(-1.3856263) q[2];
sx q[2];
rz(-0.99157809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2145591) q[1];
sx q[1];
rz(-0.3541358) q[1];
sx q[1];
rz(2.781809) q[1];
x q[2];
rz(0.59054116) q[3];
sx q[3];
rz(-1.5381375) q[3];
sx q[3];
rz(-1.6047359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6284457) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(0.15303843) q[2];
rz(0.89093527) q[3];
sx q[3];
rz(-2.1864083) q[3];
sx q[3];
rz(0.74131596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1560169) q[0];
sx q[0];
rz(-0.22933904) q[0];
sx q[0];
rz(2.7942221) q[0];
rz(0.037847606) q[1];
sx q[1];
rz(-1.9719351) q[1];
sx q[1];
rz(-2.6191424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453953) q[0];
sx q[0];
rz(-1.660083) q[0];
sx q[0];
rz(-1.9849298) q[0];
rz(-0.36549296) q[2];
sx q[2];
rz(-1.5704944) q[2];
sx q[2];
rz(0.049916849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47679893) q[1];
sx q[1];
rz(-1.069624) q[1];
sx q[1];
rz(-1.0887926) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0112004) q[3];
sx q[3];
rz(-1.921145) q[3];
sx q[3];
rz(-0.94061414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6309506) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(1.2003215) q[2];
rz(-2.5837894) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8484304) q[0];
sx q[0];
rz(-2.1577305) q[0];
sx q[0];
rz(1.4137319) q[0];
rz(0.14104715) q[1];
sx q[1];
rz(-1.7660564) q[1];
sx q[1];
rz(2.486855) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4439693) q[0];
sx q[0];
rz(-1.3338517) q[0];
sx q[0];
rz(2.9276642) q[0];
x q[1];
rz(-2.5467196) q[2];
sx q[2];
rz(-2.0030177) q[2];
sx q[2];
rz(3.1405666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84147553) q[1];
sx q[1];
rz(-1.3285331) q[1];
sx q[1];
rz(2.9650142) q[1];
rz(-0.61956866) q[3];
sx q[3];
rz(-1.1244785) q[3];
sx q[3];
rz(-0.60947641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7901788) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(0.72511017) q[2];
rz(-2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-2.5543673) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643322) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(2.019885) q[1];
sx q[1];
rz(-1.5678761) q[1];
sx q[1];
rz(-0.97490464) q[1];
rz(-0.66365343) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(-2.9573351) q[3];
sx q[3];
rz(-0.72132106) q[3];
sx q[3];
rz(-3.0527243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
