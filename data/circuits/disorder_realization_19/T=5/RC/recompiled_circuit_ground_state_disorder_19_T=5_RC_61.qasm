OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(-0.22766222) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(-1.8546606) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628107) q[0];
sx q[0];
rz(-2.037604) q[0];
sx q[0];
rz(-3.0808671) q[0];
rz(-pi) q[1];
rz(-1.1212767) q[2];
sx q[2];
rz(-1.1444905) q[2];
sx q[2];
rz(2.9232077) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70482774) q[1];
sx q[1];
rz(-2.4840762) q[1];
sx q[1];
rz(-2.2566811) q[1];
rz(-pi) q[2];
rz(1.8836796) q[3];
sx q[3];
rz(-2.1863424) q[3];
sx q[3];
rz(2.2039977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2490354) q[2];
sx q[2];
rz(-2.1302569) q[2];
sx q[2];
rz(-2.2300143) q[2];
rz(-0.75561953) q[3];
sx q[3];
rz(-2.8204462) q[3];
sx q[3];
rz(-0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78254533) q[0];
sx q[0];
rz(-1.3107212) q[0];
sx q[0];
rz(-1.1388592) q[0];
rz(-2.5149939) q[1];
sx q[1];
rz(-0.36830026) q[1];
sx q[1];
rz(-2.0638594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8846466) q[0];
sx q[0];
rz(-0.68832371) q[0];
sx q[0];
rz(-2.5907645) q[0];
rz(-pi) q[1];
rz(2.2347514) q[2];
sx q[2];
rz(-2.5996947) q[2];
sx q[2];
rz(-0.52846891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5087532) q[1];
sx q[1];
rz(-1.5037854) q[1];
sx q[1];
rz(-2.8772023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84799453) q[3];
sx q[3];
rz(-1.902033) q[3];
sx q[3];
rz(3.1183463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44800147) q[2];
sx q[2];
rz(-0.75746626) q[2];
sx q[2];
rz(-0.21270154) q[2];
rz(1.5564144) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(-1.0630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146516) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(-0.84685999) q[0];
rz(-0.2521387) q[1];
sx q[1];
rz(-2.3138901) q[1];
sx q[1];
rz(-0.65021461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668117) q[0];
sx q[0];
rz(-1.2322786) q[0];
sx q[0];
rz(1.6174497) q[0];
x q[1];
rz(-0.6619307) q[2];
sx q[2];
rz(-2.227894) q[2];
sx q[2];
rz(1.3160694) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47874988) q[1];
sx q[1];
rz(-1.2744941) q[1];
sx q[1];
rz(3.0257312) q[1];
rz(-pi) q[2];
rz(1.1232722) q[3];
sx q[3];
rz(-2.3254804) q[3];
sx q[3];
rz(-2.1838674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3113159) q[2];
sx q[2];
rz(-0.69861424) q[2];
sx q[2];
rz(-2.1777731) q[2];
rz(1.7802995) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(-0.042958766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394102) q[0];
sx q[0];
rz(-0.22265156) q[0];
sx q[0];
rz(1.0182925) q[0];
rz(0.88515431) q[1];
sx q[1];
rz(-0.2489018) q[1];
sx q[1];
rz(-0.28908602) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873661) q[0];
sx q[0];
rz(-1.2441846) q[0];
sx q[0];
rz(-1.8892889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9924329) q[2];
sx q[2];
rz(-0.9767864) q[2];
sx q[2];
rz(2.8624812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4319181) q[1];
sx q[1];
rz(-2.5920715) q[1];
sx q[1];
rz(-1.9220244) q[1];
rz(-pi) q[2];
rz(-1.698367) q[3];
sx q[3];
rz(-1.3535548) q[3];
sx q[3];
rz(-2.9590817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79124147) q[2];
sx q[2];
rz(-1.6356607) q[2];
sx q[2];
rz(1.852847) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-2.6668187) q[3];
sx q[3];
rz(0.60540664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8346005) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(2.8173764) q[0];
rz(0.038837198) q[1];
sx q[1];
rz(-0.7380929) q[1];
sx q[1];
rz(1.0968346) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457458) q[0];
sx q[0];
rz(-2.7660094) q[0];
sx q[0];
rz(-0.71806851) q[0];
rz(-1.7261502) q[2];
sx q[2];
rz(-0.58621472) q[2];
sx q[2];
rz(1.4582576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13584863) q[1];
sx q[1];
rz(-1.6127819) q[1];
sx q[1];
rz(-1.2688925) q[1];
rz(1.7750793) q[3];
sx q[3];
rz(-2.7133803) q[3];
sx q[3];
rz(-0.34303676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0754806) q[2];
sx q[2];
rz(-0.26618633) q[2];
sx q[2];
rz(-0.038662635) q[2];
rz(-1.7330811) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(-2.8601638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6283145) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(-2.1224838) q[0];
rz(-2.6844773) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(1.7105182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5605049) q[0];
sx q[0];
rz(-1.3700587) q[0];
sx q[0];
rz(1.5439347) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8879755) q[2];
sx q[2];
rz(-2.1448958) q[2];
sx q[2];
rz(-0.92437896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44164691) q[1];
sx q[1];
rz(-1.6586379) q[1];
sx q[1];
rz(1.2097174) q[1];
rz(-pi) q[2];
rz(-2.2912628) q[3];
sx q[3];
rz(-1.8763233) q[3];
sx q[3];
rz(0.014643365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21923253) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(2.5617981) q[2];
rz(2.842105) q[3];
sx q[3];
rz(-2.6087285) q[3];
sx q[3];
rz(-1.2286435) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86094588) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(-0.10375599) q[0];
rz(2.5170028) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(1.7594899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.493295) q[0];
sx q[0];
rz(-1.7837423) q[0];
sx q[0];
rz(0.84342028) q[0];
rz(2.8417009) q[2];
sx q[2];
rz(-1.6768528) q[2];
sx q[2];
rz(-0.40115717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4098072) q[1];
sx q[1];
rz(-2.5587132) q[1];
sx q[1];
rz(1.0686841) q[1];
x q[2];
rz(1.2861038) q[3];
sx q[3];
rz(-2.2119388) q[3];
sx q[3];
rz(-2.5452328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78911191) q[2];
sx q[2];
rz(-1.666297) q[2];
sx q[2];
rz(2.2769807) q[2];
rz(0.23884808) q[3];
sx q[3];
rz(-2.3819203) q[3];
sx q[3];
rz(0.71340942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637941) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(0.23502769) q[0];
rz(-2.4759953) q[1];
sx q[1];
rz(-2.0033629) q[1];
sx q[1];
rz(-1.6682909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.145133) q[0];
sx q[0];
rz(-1.0211103) q[0];
sx q[0];
rz(0.24497801) q[0];
rz(-pi) q[1];
rz(-2.6742879) q[2];
sx q[2];
rz(-1.4479785) q[2];
sx q[2];
rz(-2.2707224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43117796) q[1];
sx q[1];
rz(-1.3873552) q[1];
sx q[1];
rz(1.4967493) q[1];
rz(-0.45110945) q[3];
sx q[3];
rz(-1.6112956) q[3];
sx q[3];
rz(-0.19865741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9792446) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(2.6806504) q[2];
rz(-1.0848684) q[3];
sx q[3];
rz(-2.3960787) q[3];
sx q[3];
rz(2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0775065) q[0];
sx q[0];
rz(-0.53090799) q[0];
sx q[0];
rz(2.532646) q[0];
rz(-0.56600904) q[1];
sx q[1];
rz(-1.9809664) q[1];
sx q[1];
rz(1.9401248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.214819) q[0];
sx q[0];
rz(-1.578971) q[0];
sx q[0];
rz(-1.7063441) q[0];
rz(-1.1913774) q[2];
sx q[2];
rz(-2.622329) q[2];
sx q[2];
rz(-1.0644703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7215828) q[1];
sx q[1];
rz(-2.3370565) q[1];
sx q[1];
rz(-1.700042) q[1];
x q[2];
rz(-2.5646832) q[3];
sx q[3];
rz(-0.52191496) q[3];
sx q[3];
rz(-1.0229223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0386049) q[2];
sx q[2];
rz(-2.0709585) q[2];
sx q[2];
rz(-1.1919682) q[2];
rz(-1.0209171) q[3];
sx q[3];
rz(-0.20191419) q[3];
sx q[3];
rz(1.6010223) q[3];
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
rz(-2.9182619) q[0];
sx q[0];
rz(-0.20063618) q[0];
sx q[0];
rz(2.1740792) q[0];
rz(-1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(2.7593625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8865693) q[0];
sx q[0];
rz(-0.60211997) q[0];
sx q[0];
rz(1.894879) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76127865) q[2];
sx q[2];
rz(-0.57028162) q[2];
sx q[2];
rz(-2.7340661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8763242) q[1];
sx q[1];
rz(-2.3493715) q[1];
sx q[1];
rz(-1.7382311) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5955673) q[3];
sx q[3];
rz(-1.5848397) q[3];
sx q[3];
rz(-1.394681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.050015673) q[2];
sx q[2];
rz(-2.2735368) q[2];
sx q[2];
rz(-0.78563219) q[2];
rz(3.1146289) q[3];
sx q[3];
rz(-1.6227159) q[3];
sx q[3];
rz(-0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80326573) q[0];
sx q[0];
rz(-1.4563518) q[0];
sx q[0];
rz(-0.8932054) q[0];
rz(-2.4772353) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(-1.8705838) q[2];
sx q[2];
rz(-0.86465093) q[2];
sx q[2];
rz(1.6563889) q[2];
rz(-1.5803278) q[3];
sx q[3];
rz(-0.96612488) q[3];
sx q[3];
rz(0.92074367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
