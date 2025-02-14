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
rz(0.83346811) q[0];
sx q[0];
rz(-0.8136189) q[0];
sx q[0];
rz(1.5746434) q[0];
rz(-1.1235224) q[1];
sx q[1];
rz(-0.82668537) q[1];
sx q[1];
rz(-2.6348662) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5187162) q[0];
sx q[0];
rz(-2.5883247) q[0];
sx q[0];
rz(3.1068463) q[0];
rz(-pi) q[1];
rz(1.1415403) q[2];
sx q[2];
rz(-0.55396336) q[2];
sx q[2];
rz(-2.9316154) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4611257) q[1];
sx q[1];
rz(-2.8790747) q[1];
sx q[1];
rz(0.21204388) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1024974) q[3];
sx q[3];
rz(-2.2347898) q[3];
sx q[3];
rz(-1.853626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7099521) q[2];
sx q[2];
rz(-1.9441354) q[2];
sx q[2];
rz(0.54090995) q[2];
rz(-1.8139808) q[3];
sx q[3];
rz(-1.4541459) q[3];
sx q[3];
rz(-2.4732164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79880181) q[0];
sx q[0];
rz(-1.7690161) q[0];
sx q[0];
rz(-2.0846682) q[0];
rz(-2.7150555) q[1];
sx q[1];
rz(-1.5142454) q[1];
sx q[1];
rz(3.1372517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2452702) q[0];
sx q[0];
rz(-1.9548804) q[0];
sx q[0];
rz(-0.86309582) q[0];
rz(-pi) q[1];
rz(-0.40965457) q[2];
sx q[2];
rz(-2.1623297) q[2];
sx q[2];
rz(-0.17342527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2935644) q[1];
sx q[1];
rz(-1.5604158) q[1];
sx q[1];
rz(1.4917828) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0375157) q[3];
sx q[3];
rz(-0.53734578) q[3];
sx q[3];
rz(1.6598824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9428375) q[2];
sx q[2];
rz(-1.7386856) q[2];
sx q[2];
rz(-0.37503654) q[2];
rz(0.047920553) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(0.16938773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.38309836) q[0];
sx q[0];
rz(-2.5653745) q[0];
sx q[0];
rz(1.6297485) q[0];
rz(-1.0668628) q[1];
sx q[1];
rz(-1.2991512) q[1];
sx q[1];
rz(-0.80924928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3602995) q[0];
sx q[0];
rz(-1.5341833) q[0];
sx q[0];
rz(0.039860241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.079945076) q[2];
sx q[2];
rz(-2.4221651) q[2];
sx q[2];
rz(1.5647174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7748562) q[1];
sx q[1];
rz(-1.3778169) q[1];
sx q[1];
rz(-2.3425413) q[1];
rz(-0.0076313) q[3];
sx q[3];
rz(-2.7797709) q[3];
sx q[3];
rz(-1.0135038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2201436) q[2];
sx q[2];
rz(-2.4704832) q[2];
sx q[2];
rz(-1.5044093) q[2];
rz(0.97357059) q[3];
sx q[3];
rz(-0.74876553) q[3];
sx q[3];
rz(-2.0481295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022738986) q[0];
sx q[0];
rz(-1.5596507) q[0];
sx q[0];
rz(2.7795025) q[0];
rz(-1.3478966) q[1];
sx q[1];
rz(-1.7575512) q[1];
sx q[1];
rz(1.0628343) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14236951) q[0];
sx q[0];
rz(-2.9175715) q[0];
sx q[0];
rz(-2.260452) q[0];
x q[1];
rz(1.2591909) q[2];
sx q[2];
rz(-0.62436337) q[2];
sx q[2];
rz(-1.7758689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0697067) q[1];
sx q[1];
rz(-0.45821689) q[1];
sx q[1];
rz(0.86104844) q[1];
rz(-pi) q[2];
rz(-2.8873467) q[3];
sx q[3];
rz(-1.3602453) q[3];
sx q[3];
rz(-2.5894534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7245543) q[2];
sx q[2];
rz(-0.78196708) q[2];
sx q[2];
rz(-3.0843201) q[2];
rz(-0.96632424) q[3];
sx q[3];
rz(-2.1939907) q[3];
sx q[3];
rz(1.2053325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.22369497) q[0];
sx q[0];
rz(-0.88858336) q[0];
sx q[0];
rz(2.5433871) q[0];
rz(1.573645) q[1];
sx q[1];
rz(-1.6697829) q[1];
sx q[1];
rz(-0.15984687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1956045) q[0];
sx q[0];
rz(-0.78778446) q[0];
sx q[0];
rz(-2.8986918) q[0];
rz(-0.21128006) q[2];
sx q[2];
rz(-1.5108372) q[2];
sx q[2];
rz(-2.243239) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7084078) q[1];
sx q[1];
rz(-1.9497279) q[1];
sx q[1];
rz(-0.66431944) q[1];
rz(-3.065291) q[3];
sx q[3];
rz(-1.6152744) q[3];
sx q[3];
rz(-0.22212576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9827031) q[2];
sx q[2];
rz(-0.44124678) q[2];
sx q[2];
rz(0.59054792) q[2];
rz(-2.1741137) q[3];
sx q[3];
rz(-1.5524105) q[3];
sx q[3];
rz(-1.9578741) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972873) q[0];
sx q[0];
rz(-2.4794281) q[0];
sx q[0];
rz(-2.1024607) q[0];
rz(0.93152535) q[1];
sx q[1];
rz(-2.4577699) q[1];
sx q[1];
rz(-2.0328124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174579) q[0];
sx q[0];
rz(-2.0724247) q[0];
sx q[0];
rz(-2.5162016) q[0];
rz(-pi) q[1];
rz(0.80486561) q[2];
sx q[2];
rz(-2.1535843) q[2];
sx q[2];
rz(1.031176) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0045864) q[1];
sx q[1];
rz(-0.69939628) q[1];
sx q[1];
rz(-0.57199332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1734725) q[3];
sx q[3];
rz(-0.58486551) q[3];
sx q[3];
rz(-2.7871451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80591566) q[2];
sx q[2];
rz(-0.47199619) q[2];
sx q[2];
rz(-2.8373888) q[2];
rz(0.16600569) q[3];
sx q[3];
rz(-1.3622354) q[3];
sx q[3];
rz(-1.5267052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6840927) q[0];
sx q[0];
rz(-2.9170211) q[0];
sx q[0];
rz(1.7042879) q[0];
rz(-2.8999088) q[1];
sx q[1];
rz(-1.6432523) q[1];
sx q[1];
rz(2.2041352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3110549) q[0];
sx q[0];
rz(-1.3321345) q[0];
sx q[0];
rz(-0.72913362) q[0];
rz(-pi) q[1];
x q[1];
rz(0.030539837) q[2];
sx q[2];
rz(-1.2903155) q[2];
sx q[2];
rz(-2.7860093) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2087349) q[1];
sx q[1];
rz(-0.94122771) q[1];
sx q[1];
rz(0.97811326) q[1];
x q[2];
rz(-0.44643965) q[3];
sx q[3];
rz(-1.0390128) q[3];
sx q[3];
rz(1.4908294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23718701) q[2];
sx q[2];
rz(-2.0439456) q[2];
sx q[2];
rz(-0.22906765) q[2];
rz(1.1711052) q[3];
sx q[3];
rz(-1.3874715) q[3];
sx q[3];
rz(-0.68172014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218629) q[0];
sx q[0];
rz(-0.64058146) q[0];
sx q[0];
rz(1.4135452) q[0];
rz(1.9238663) q[1];
sx q[1];
rz(-1.7037337) q[1];
sx q[1];
rz(-2.6103643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7290658) q[0];
sx q[0];
rz(-2.7587037) q[0];
sx q[0];
rz(-1.0273496) q[0];
x q[1];
rz(0.81905535) q[2];
sx q[2];
rz(-1.2651288) q[2];
sx q[2];
rz(-1.7421987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8063691) q[1];
sx q[1];
rz(-2.1861736) q[1];
sx q[1];
rz(-2.1187617) q[1];
rz(-0.30398627) q[3];
sx q[3];
rz(-1.8517947) q[3];
sx q[3];
rz(-1.3490647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5765877) q[2];
sx q[2];
rz(-1.7030623) q[2];
sx q[2];
rz(0.79565221) q[2];
rz(2.3203881) q[3];
sx q[3];
rz(-2.932565) q[3];
sx q[3];
rz(0.30278277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32383701) q[0];
sx q[0];
rz(-2.9471286) q[0];
sx q[0];
rz(-1.4870462) q[0];
rz(1.6200292) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(1.0618658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855536) q[0];
sx q[0];
rz(-1.0526219) q[0];
sx q[0];
rz(-1.1382301) q[0];
rz(-pi) q[1];
rz(-1.7369607) q[2];
sx q[2];
rz(-0.29815692) q[2];
sx q[2];
rz(0.96337897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3848906) q[1];
sx q[1];
rz(-1.1497918) q[1];
sx q[1];
rz(0.16709354) q[1];
rz(-pi) q[2];
rz(1.1561772) q[3];
sx q[3];
rz(-1.3870169) q[3];
sx q[3];
rz(1.6971677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1294641) q[2];
sx q[2];
rz(-2.0120554) q[2];
sx q[2];
rz(3.0322292) q[2];
rz(0.28956413) q[3];
sx q[3];
rz(-1.3907631) q[3];
sx q[3];
rz(-0.65350986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7107723) q[0];
sx q[0];
rz(-1.1752952) q[0];
sx q[0];
rz(-1.9868504) q[0];
rz(2.7180502) q[1];
sx q[1];
rz(-1.4332899) q[1];
sx q[1];
rz(2.4636041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38289878) q[0];
sx q[0];
rz(-0.80307559) q[0];
sx q[0];
rz(2.4944958) q[0];
rz(-pi) q[1];
rz(-2.8805947) q[2];
sx q[2];
rz(-1.3512865) q[2];
sx q[2];
rz(-1.0936119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0748234) q[1];
sx q[1];
rz(-0.69157234) q[1];
sx q[1];
rz(2.2645386) q[1];
rz(-0.5918233) q[3];
sx q[3];
rz(-1.4524959) q[3];
sx q[3];
rz(0.59676701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8573711) q[2];
sx q[2];
rz(-2.325433) q[2];
sx q[2];
rz(-0.59874272) q[2];
rz(1.5395928) q[3];
sx q[3];
rz(-1.9339823) q[3];
sx q[3];
rz(-2.1046751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4620062) q[0];
sx q[0];
rz(-2.7988667) q[0];
sx q[0];
rz(-1.3890247) q[0];
rz(-3.1287843) q[1];
sx q[1];
rz(-2.7638331) q[1];
sx q[1];
rz(1.4889181) q[1];
rz(2.2118582) q[2];
sx q[2];
rz(-1.6597181) q[2];
sx q[2];
rz(1.5372288) q[2];
rz(-1.011856) q[3];
sx q[3];
rz(-1.5717117) q[3];
sx q[3];
rz(-1.6757884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
