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
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(6.3732014) q[1];
sx q[1];
rz(6.7572588) q[1];
sx q[1];
rz(13.401539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157566) q[0];
sx q[0];
rz(-0.8716363) q[0];
sx q[0];
rz(1.9833724) q[0];
rz(-pi) q[1];
rz(-2.7273601) q[2];
sx q[2];
rz(-2.4576839) q[2];
sx q[2];
rz(1.6773083) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7427828) q[1];
sx q[1];
rz(-1.2025857) q[1];
sx q[1];
rz(1.8333992) q[1];
rz(-pi) q[2];
rz(-0.99906875) q[3];
sx q[3];
rz(-2.1474827) q[3];
sx q[3];
rz(-0.028923361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0498672) q[2];
sx q[2];
rz(-2.2455402) q[2];
sx q[2];
rz(0.33207616) q[2];
rz(-0.24886985) q[3];
sx q[3];
rz(-1.1955465) q[3];
sx q[3];
rz(2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2700972) q[0];
sx q[0];
rz(-2.8506554) q[0];
sx q[0];
rz(0.53263295) q[0];
rz(-3.0337785) q[1];
sx q[1];
rz(-1.1153406) q[1];
sx q[1];
rz(2.8210988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924393) q[0];
sx q[0];
rz(-2.4863613) q[0];
sx q[0];
rz(-1.4661319) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8024496) q[2];
sx q[2];
rz(-1.4024874) q[2];
sx q[2];
rz(-1.2141808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35831603) q[1];
sx q[1];
rz(-0.22138295) q[1];
sx q[1];
rz(2.3795147) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3555834) q[3];
sx q[3];
rz(-1.8930515) q[3];
sx q[3];
rz(-0.93860579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8104441) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(1.6395052) q[2];
rz(-0.33809996) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(-2.6147208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51760393) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(0.35743085) q[0];
rz(-0.39237157) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(1.9270814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314963) q[0];
sx q[0];
rz(-1.7123245) q[0];
sx q[0];
rz(3.1241199) q[0];
x q[1];
rz(2.8896595) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(2.6643945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4068661) q[1];
sx q[1];
rz(-2.7541783) q[1];
sx q[1];
rz(-1.5785567) q[1];
rz(-1.6271934) q[3];
sx q[3];
rz(-2.5531451) q[3];
sx q[3];
rz(1.9533039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3452722) q[2];
sx q[2];
rz(-2.1638963) q[2];
sx q[2];
rz(-2.7447682) q[2];
rz(-1.2567629) q[3];
sx q[3];
rz(-1.7233012) q[3];
sx q[3];
rz(0.61029148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20763718) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(1.9130094) q[0];
rz(1.928891) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(2.520715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29781326) q[0];
sx q[0];
rz(-0.50619805) q[0];
sx q[0];
rz(-1.9202581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3880896) q[2];
sx q[2];
rz(-0.22543365) q[2];
sx q[2];
rz(1.4788675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.32737) q[1];
sx q[1];
rz(-1.3124221) q[1];
sx q[1];
rz(1.8381005) q[1];
rz(1.4520763) q[3];
sx q[3];
rz(-1.190589) q[3];
sx q[3];
rz(0.79352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2584201) q[2];
sx q[2];
rz(-2.0032538) q[2];
sx q[2];
rz(-2.9849198) q[2];
rz(0.91642085) q[3];
sx q[3];
rz(-1.0333002) q[3];
sx q[3];
rz(2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(0.39749417) q[0];
rz(-0.022857895) q[1];
sx q[1];
rz(-2.6409179) q[1];
sx q[1];
rz(2.4668677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0894273) q[0];
sx q[0];
rz(-0.31712714) q[0];
sx q[0];
rz(0.84256147) q[0];
x q[1];
rz(1.6653453) q[2];
sx q[2];
rz(-1.6942021) q[2];
sx q[2];
rz(-2.9837554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6276363) q[1];
sx q[1];
rz(-1.4003203) q[1];
sx q[1];
rz(-1.7254616) q[1];
x q[2];
rz(-0.27366112) q[3];
sx q[3];
rz(-1.869264) q[3];
sx q[3];
rz(-0.65200114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(0.69925365) q[2];
rz(1.3462542) q[3];
sx q[3];
rz(-0.56763879) q[3];
sx q[3];
rz(-2.4301372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(0.97459546) q[1];
sx q[1];
rz(-1.3037126) q[1];
sx q[1];
rz(1.7030565) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33773003) q[0];
sx q[0];
rz(-0.32050214) q[0];
sx q[0];
rz(2.8021373) q[0];
x q[1];
rz(1.3149777) q[2];
sx q[2];
rz(-1.5808357) q[2];
sx q[2];
rz(-2.5217587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12425646) q[1];
sx q[1];
rz(-1.8030232) q[1];
sx q[1];
rz(-0.88084765) q[1];
rz(-0.4532279) q[3];
sx q[3];
rz(-2.6010644) q[3];
sx q[3];
rz(1.5858397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67821104) q[2];
sx q[2];
rz(-1.7794098) q[2];
sx q[2];
rz(-1.8360651) q[2];
rz(-1.8259004) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(2.2731884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797794) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(1.6449991) q[0];
rz(0.98622259) q[1];
sx q[1];
rz(-1.58135) q[1];
sx q[1];
rz(0.49759069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80579306) q[0];
sx q[0];
rz(-1.8247274) q[0];
sx q[0];
rz(0.098761244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1946218) q[2];
sx q[2];
rz(-1.1188036) q[2];
sx q[2];
rz(2.3561881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.033173) q[1];
sx q[1];
rz(-2.7153176) q[1];
sx q[1];
rz(1.1854965) q[1];
rz(-pi) q[2];
rz(1.6646181) q[3];
sx q[3];
rz(-1.854343) q[3];
sx q[3];
rz(0.16167262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30279532) q[2];
sx q[2];
rz(-1.5668198) q[2];
sx q[2];
rz(-0.29067579) q[2];
rz(2.931328) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(2.1073585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57598376) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(1.1822816) q[0];
rz(-0.15444175) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(2.2363037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8715031) q[0];
sx q[0];
rz(-0.2729899) q[0];
sx q[0];
rz(-1.6156625) q[0];
rz(-0.45640517) q[2];
sx q[2];
rz(-1.3564566) q[2];
sx q[2];
rz(0.53949196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0539848) q[1];
sx q[1];
rz(-1.4248669) q[1];
sx q[1];
rz(-2.3838359) q[1];
rz(-pi) q[2];
rz(-1.0175627) q[3];
sx q[3];
rz(-2.746547) q[3];
sx q[3];
rz(1.0329122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0514544) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(-0.95019379) q[2];
rz(-0.072362445) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(0.70934057) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224834) q[0];
sx q[0];
rz(-0.95471946) q[0];
sx q[0];
rz(2.0461653) q[0];
rz(1.0911881) q[1];
sx q[1];
rz(-0.90091101) q[1];
sx q[1];
rz(-0.99517623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54992095) q[0];
sx q[0];
rz(-2.4732145) q[0];
sx q[0];
rz(-0.10381283) q[0];
x q[1];
rz(2.2500317) q[2];
sx q[2];
rz(-2.7195647) q[2];
sx q[2];
rz(2.0419378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0867445) q[1];
sx q[1];
rz(-1.5252171) q[1];
sx q[1];
rz(0.084910141) q[1];
rz(-pi) q[2];
rz(0.82251151) q[3];
sx q[3];
rz(-2.872618) q[3];
sx q[3];
rz(-1.9654056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5114078) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(2.4616145) q[2];
rz(-1.6857111) q[3];
sx q[3];
rz(-2.2472491) q[3];
sx q[3];
rz(1.9623914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.7162914) q[0];
sx q[0];
rz(-0.93451262) q[0];
sx q[0];
rz(-2.5262078) q[0];
rz(-1.3353434) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(1.3478442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6123209) q[0];
sx q[0];
rz(-2.7833496) q[0];
sx q[0];
rz(-0.1310346) q[0];
rz(-pi) q[1];
rz(-0.87839076) q[2];
sx q[2];
rz(-2.2579402) q[2];
sx q[2];
rz(2.2855482) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7536698) q[1];
sx q[1];
rz(-1.7311454) q[1];
sx q[1];
rz(2.760375) q[1];
rz(-1.4487292) q[3];
sx q[3];
rz(-2.3231835) q[3];
sx q[3];
rz(-0.19644745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24796692) q[2];
sx q[2];
rz(-1.8411571) q[2];
sx q[2];
rz(0.51353961) q[2];
rz(-0.14704554) q[3];
sx q[3];
rz(-0.5439609) q[3];
sx q[3];
rz(-1.2646382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2697987) q[0];
sx q[0];
rz(-0.88903058) q[0];
sx q[0];
rz(-0.24118184) q[0];
rz(-1.3006032) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(1.945449) q[2];
sx q[2];
rz(-1.4383264) q[2];
sx q[2];
rz(-1.2551398) q[2];
rz(1.1250238) q[3];
sx q[3];
rz(-0.99109886) q[3];
sx q[3];
rz(-2.9352321) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
