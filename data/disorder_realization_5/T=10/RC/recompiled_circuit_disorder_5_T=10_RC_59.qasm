OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(4.7441626) q[0];
sx q[0];
rz(6.8466853) q[0];
sx q[0];
rz(9.8817649) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5029992) q[0];
sx q[0];
rz(-0.25829691) q[0];
sx q[0];
rz(0.55708416) q[0];
rz(-pi) q[1];
rz(-1.4957501) q[2];
sx q[2];
rz(-1.6610166) q[2];
sx q[2];
rz(-0.98352945) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16986019) q[1];
sx q[1];
rz(-1.5105376) q[1];
sx q[1];
rz(-1.3678642) q[1];
rz(2.071322) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(0.69858944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-2.9764552) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-2.813208) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-2.3388458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734056) q[0];
sx q[0];
rz(-1.7965839) q[0];
sx q[0];
rz(-1.9682103) q[0];
x q[1];
rz(1.0674092) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(-0.78546055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3992577) q[1];
sx q[1];
rz(-0.62082532) q[1];
sx q[1];
rz(-1.4820815) q[1];
rz(-pi) q[2];
rz(-1.4071776) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(2.8676652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-0.96174812) q[3];
sx q[3];
rz(0.081993016) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(-0.28513518) q[0];
rz(-0.51672283) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8782608) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(1.0017298) q[0];
rz(2.2604733) q[2];
sx q[2];
rz(-1.3291877) q[2];
sx q[2];
rz(1.9067665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9078319) q[1];
sx q[1];
rz(-1.7572174) q[1];
sx q[1];
rz(2.9790661) q[1];
x q[2];
rz(0.47311802) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(0.21425042) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-0.66657153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.727107) q[0];
sx q[0];
rz(-1.3654728) q[0];
sx q[0];
rz(1.3798825) q[0];
x q[1];
rz(0.71422691) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(-0.54745882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2506634) q[1];
sx q[1];
rz(-2.2729985) q[1];
sx q[1];
rz(2.3198747) q[1];
rz(-2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(0.028545054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-2.4812223) q[2];
rz(1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3876622) q[0];
sx q[0];
rz(-1.8844863) q[0];
sx q[0];
rz(0.73648209) q[0];
x q[1];
rz(2.1690364) q[2];
sx q[2];
rz(-1.4417357) q[2];
sx q[2];
rz(-1.8051403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0598462) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(2.5043143) q[3];
sx q[3];
rz(-2.0162458) q[3];
sx q[3];
rz(0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(-0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(3.040722) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4911762) q[0];
sx q[0];
rz(-0.6498148) q[0];
sx q[0];
rz(-0.25214904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(3.0234408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(-0.15921758) q[1];
rz(-pi) q[2];
rz(2.7363339) q[3];
sx q[3];
rz(-1.7856981) q[3];
sx q[3];
rz(-0.98291558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4042525) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(-2.960079) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.293175) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(-0.68774736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15935414) q[0];
sx q[0];
rz(-1.9040477) q[0];
sx q[0];
rz(2.9100145) q[0];
x q[1];
rz(-1.8887927) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(-0.064934864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3878165) q[1];
sx q[1];
rz(-3.0101335) q[1];
sx q[1];
rz(-0.17863518) q[1];
x q[2];
rz(1.0024768) q[3];
sx q[3];
rz(-1.5573313) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3380276) q[2];
sx q[2];
rz(-1.4146718) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872902) q[0];
sx q[0];
rz(-1.1371005) q[0];
sx q[0];
rz(-1.4951697) q[0];
rz(1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(1.5944634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63003507) q[1];
sx q[1];
rz(-2.6361598) q[1];
sx q[1];
rz(-3.0384484) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46320398) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(-2.6686058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(2.3525227) q[2];
rz(2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-0.58832204) q[0];
rz(-3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(-2.2081597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(-1.6305001) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2869296) q[2];
sx q[2];
rz(-0.48644201) q[2];
sx q[2];
rz(1.2235299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.583608) q[1];
sx q[1];
rz(-2.4568111) q[1];
sx q[1];
rz(-0.058137356) q[1];
x q[2];
rz(-1.2133573) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(0.77222094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.6588541) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-2.488234) q[0];
rz(-1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-0.18383372) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375658) q[0];
sx q[0];
rz(-1.6871534) q[0];
sx q[0];
rz(-3.0202306) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3751971) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(-0.11608427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2059584) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(0.88840719) q[1];
rz(-pi) q[2];
rz(1.5552181) q[3];
sx q[3];
rz(-0.9874978) q[3];
sx q[3];
rz(2.616236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.2349293) q[3];
sx q[3];
rz(1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-0.72485483) q[2];
sx q[2];
rz(-2.5584494) q[2];
sx q[2];
rz(0.7128788) q[2];
rz(-0.17584569) q[3];
sx q[3];
rz(-0.98777117) q[3];
sx q[3];
rz(-2.5143757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];