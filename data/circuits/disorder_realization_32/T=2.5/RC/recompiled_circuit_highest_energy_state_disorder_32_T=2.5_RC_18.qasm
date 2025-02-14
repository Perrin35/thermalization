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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9204694) q[0];
sx q[0];
rz(-1.2386432) q[0];
sx q[0];
rz(1.8355379) q[0];
rz(-0.31754874) q[2];
sx q[2];
rz(-2.5032836) q[2];
sx q[2];
rz(-2.2090863) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3598944) q[1];
sx q[1];
rz(-0.89907032) q[1];
sx q[1];
rz(-2.2132378) q[1];
rz(-pi) q[2];
rz(1.3377519) q[3];
sx q[3];
rz(-0.50737689) q[3];
sx q[3];
rz(2.9180067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0872385) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(-0.13620201) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-2.4583702) q[3];
sx q[3];
rz(2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1012652) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-2.878317) q[0];
rz(2.1030262) q[1];
sx q[1];
rz(-0.34663215) q[1];
sx q[1];
rz(-0.77478772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8717125) q[0];
sx q[0];
rz(-2.1795666) q[0];
sx q[0];
rz(1.3193498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2293533) q[2];
sx q[2];
rz(-0.99102913) q[2];
sx q[2];
rz(-1.6032246) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9866219) q[1];
sx q[1];
rz(-0.38119527) q[1];
sx q[1];
rz(-1.6705369) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0618013) q[3];
sx q[3];
rz(-0.92872075) q[3];
sx q[3];
rz(-1.1428558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7210377) q[2];
sx q[2];
rz(-1.5169531) q[2];
sx q[2];
rz(0.59845412) q[2];
rz(1.7789486) q[3];
sx q[3];
rz(-0.88038954) q[3];
sx q[3];
rz(2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(-1.5077952) q[0];
rz(-2.9871121) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(2.3522164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95473991) q[0];
sx q[0];
rz(-0.53760872) q[0];
sx q[0];
rz(-0.16969271) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46044431) q[2];
sx q[2];
rz(-0.31184648) q[2];
sx q[2];
rz(2.0058035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.912108) q[1];
sx q[1];
rz(-0.69522017) q[1];
sx q[1];
rz(-1.0840313) q[1];
rz(-pi) q[2];
rz(-0.79929914) q[3];
sx q[3];
rz(-1.2033312) q[3];
sx q[3];
rz(2.7225843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35884759) q[2];
sx q[2];
rz(-0.91155702) q[2];
sx q[2];
rz(1.6645128) q[2];
rz(1.2137671) q[3];
sx q[3];
rz(-1.341308) q[3];
sx q[3];
rz(1.7267797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724991) q[0];
sx q[0];
rz(-1.0558244) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(-0.67131132) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(-1.7866561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63263921) q[0];
sx q[0];
rz(-1.1221948) q[0];
sx q[0];
rz(2.6030356) q[0];
x q[1];
rz(-0.12570555) q[2];
sx q[2];
rz(-1.9393688) q[2];
sx q[2];
rz(1.5154132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5004394) q[1];
sx q[1];
rz(-1.8158923) q[1];
sx q[1];
rz(-3.0768454) q[1];
rz(-pi) q[2];
rz(1.4282966) q[3];
sx q[3];
rz(-2.550052) q[3];
sx q[3];
rz(-1.9672036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9352202) q[2];
sx q[2];
rz(-0.78038961) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(2.3937461) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(-1.4486754) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.185323) q[0];
sx q[0];
rz(-1.3893501) q[0];
sx q[0];
rz(-2.4638033) q[0];
rz(0.14450821) q[1];
sx q[1];
rz(-2.3542207) q[1];
sx q[1];
rz(1.4871303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5337473) q[0];
sx q[0];
rz(-0.13838875) q[0];
sx q[0];
rz(1.7615303) q[0];
x q[1];
rz(1.4356639) q[2];
sx q[2];
rz(-0.40529521) q[2];
sx q[2];
rz(0.43650337) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9095535) q[1];
sx q[1];
rz(-2.0318446) q[1];
sx q[1];
rz(2.1559245) q[1];
rz(-pi) q[2];
x q[2];
rz(0.048154496) q[3];
sx q[3];
rz(-1.1954525) q[3];
sx q[3];
rz(-2.1233692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80374485) q[2];
sx q[2];
rz(-0.2636815) q[2];
sx q[2];
rz(2.7860876) q[2];
rz(2.096094) q[3];
sx q[3];
rz(-1.9020566) q[3];
sx q[3];
rz(-2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97583714) q[0];
sx q[0];
rz(-0.58335692) q[0];
sx q[0];
rz(-3.052886) q[0];
rz(-1.5345796) q[1];
sx q[1];
rz(-1.8314223) q[1];
sx q[1];
rz(-0.075693695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450485) q[0];
sx q[0];
rz(-2.2296784) q[0];
sx q[0];
rz(1.8152899) q[0];
rz(-pi) q[1];
rz(1.5140818) q[2];
sx q[2];
rz(-0.73243388) q[2];
sx q[2];
rz(-2.5947941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5867501) q[1];
sx q[1];
rz(-1.7126709) q[1];
sx q[1];
rz(1.0125404) q[1];
rz(-pi) q[2];
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
rz(-0.82464108) q[2];
sx q[2];
rz(-1.3836766) q[2];
sx q[2];
rz(-1.0392044) q[2];
rz(2.9292987) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(-1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742663) q[0];
sx q[0];
rz(-0.75263158) q[0];
sx q[0];
rz(1.0443895) q[0];
rz(2.3692865) q[1];
sx q[1];
rz(-0.65985313) q[1];
sx q[1];
rz(-0.73371249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334856) q[0];
sx q[0];
rz(-2.2758188) q[0];
sx q[0];
rz(0.66410983) q[0];
x q[1];
rz(-0.15007054) q[2];
sx q[2];
rz(-1.9378621) q[2];
sx q[2];
rz(-2.31524) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1024688) q[1];
sx q[1];
rz(-2.8274635) q[1];
sx q[1];
rz(2.7565095) q[1];
rz(0.10802631) q[3];
sx q[3];
rz(-0.83893925) q[3];
sx q[3];
rz(2.9895003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61651984) q[2];
sx q[2];
rz(-0.51023444) q[2];
sx q[2];
rz(0.60866848) q[2];
rz(-1.0673374) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(-2.5909891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6390425) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(1.6453561) q[0];
rz(-1.3642338) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-0.38633698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3583806) q[0];
sx q[0];
rz(-0.9045426) q[0];
sx q[0];
rz(-1.233564) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.069938439) q[2];
sx q[2];
rz(-1.3856263) q[2];
sx q[2];
rz(-2.1500146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.54562927) q[1];
sx q[1];
rz(-1.2402322) q[1];
sx q[1];
rz(-1.7002374) q[1];
rz(-pi) q[2];
rz(-1.6101077) q[3];
sx q[3];
rz(-2.1609801) q[3];
sx q[3];
rz(-3.085768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6284457) q[2];
sx q[2];
rz(-1.7690423) q[2];
sx q[2];
rz(-2.9885542) q[2];
rz(-0.89093527) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(-2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1560169) q[0];
sx q[0];
rz(-2.9122536) q[0];
sx q[0];
rz(-0.34737059) q[0];
rz(-3.103745) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(-0.52245021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1570826) q[0];
sx q[0];
rz(-1.9831796) q[0];
sx q[0];
rz(-0.097481485) q[0];
rz(-pi) q[1];
rz(-2.7760997) q[2];
sx q[2];
rz(-1.5704944) q[2];
sx q[2];
rz(-0.049916849) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3047239) q[1];
sx q[1];
rz(-2.4608399) q[1];
sx q[1];
rz(-2.4393243) q[1];
rz(-1.0112004) q[3];
sx q[3];
rz(-1.2204477) q[3];
sx q[3];
rz(0.94061414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5106421) q[2];
sx q[2];
rz(-0.46322552) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(-0.55780324) q[3];
sx q[3];
rz(-1.5188981) q[3];
sx q[3];
rz(-0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2931622) q[0];
sx q[0];
rz(-2.1577305) q[0];
sx q[0];
rz(-1.4137319) q[0];
rz(-0.14104715) q[1];
sx q[1];
rz(-1.7660564) q[1];
sx q[1];
rz(-2.486855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6976234) q[0];
sx q[0];
rz(-1.3338517) q[0];
sx q[0];
rz(-0.21392845) q[0];
rz(-pi) q[1];
rz(-0.68871542) q[2];
sx q[2];
rz(-0.71962269) q[2];
sx q[2];
rz(-2.1261499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68653983) q[1];
sx q[1];
rz(-1.3994263) q[1];
sx q[1];
rz(-1.3248578) q[1];
rz(-pi) q[2];
rz(0.68924381) q[3];
sx q[3];
rz(-0.74604496) q[3];
sx q[3];
rz(-0.41710284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(2.4164825) q[2];
rz(2.2509947) q[3];
sx q[3];
rz(-0.85024873) q[3];
sx q[3];
rz(-2.5543673) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17726041) q[0];
sx q[0];
rz(-1.6163419) q[0];
sx q[0];
rz(1.1905715) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(-0.66365343) q[2];
sx q[2];
rz(-1.2812231) q[2];
sx q[2];
rz(-0.44322586) q[2];
rz(-0.18425758) q[3];
sx q[3];
rz(-2.4202716) q[3];
sx q[3];
rz(0.088868401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
