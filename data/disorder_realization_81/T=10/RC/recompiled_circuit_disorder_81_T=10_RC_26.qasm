OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(-1.7629495) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9076951) q[0];
sx q[0];
rz(-1.6595708) q[0];
sx q[0];
rz(1.5855473) q[0];
x q[1];
rz(-1.7832463) q[2];
sx q[2];
rz(-2.6531086) q[2];
sx q[2];
rz(0.24220315) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3528459) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(2.7098141) q[1];
rz(-pi) q[2];
rz(2.2545635) q[3];
sx q[3];
rz(-1.6449882) q[3];
sx q[3];
rz(2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(-2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79016722) q[0];
sx q[0];
rz(-1.5895491) q[0];
sx q[0];
rz(0.055939527) q[0];
rz(0.29859121) q[2];
sx q[2];
rz(-2.9332187) q[2];
sx q[2];
rz(-1.6076455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8801404) q[1];
sx q[1];
rz(-1.6958691) q[1];
sx q[1];
rz(-2.2373881) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69497377) q[3];
sx q[3];
rz(-1.4222099) q[3];
sx q[3];
rz(2.9050764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(-2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3804669) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-0.82021964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41528156) q[0];
sx q[0];
rz(-1.6831241) q[0];
sx q[0];
rz(0.88322722) q[0];
x q[1];
rz(2.064346) q[2];
sx q[2];
rz(-1.764467) q[2];
sx q[2];
rz(0.94656241) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0606544) q[1];
sx q[1];
rz(-1.2186236) q[1];
sx q[1];
rz(-1.1206131) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0355755) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(2.9860949) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-2.2303175) q[0];
rz(0.43831929) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.320425) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94538044) q[0];
sx q[0];
rz(-1.5811265) q[0];
sx q[0];
rz(1.1611847) q[0];
x q[1];
rz(-0.40789149) q[2];
sx q[2];
rz(-2.6558999) q[2];
sx q[2];
rz(2.8913468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3418386) q[1];
sx q[1];
rz(-0.8952039) q[1];
sx q[1];
rz(-1.9104596) q[1];
x q[2];
rz(-2.5370595) q[3];
sx q[3];
rz(-2.2570838) q[3];
sx q[3];
rz(-1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(1.140973) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300366) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.226549) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.3006166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039935) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(-3.0689737) q[0];
x q[1];
rz(1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-0.14771151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(0.57002108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8793094) q[3];
sx q[3];
rz(-1.7119006) q[3];
sx q[3];
rz(1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(2.208476) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(1.5531497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827305) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(-1.9431252) q[0];
x q[1];
rz(-1.2049098) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(-2.6423955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99332422) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(0.12970129) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2736069) q[3];
sx q[3];
rz(-1.6276974) q[3];
sx q[3];
rz(-1.1662607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(2.4643331) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(0.97704926) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1099694) q[0];
sx q[0];
rz(-0.97951802) q[0];
sx q[0];
rz(1.6615608) q[0];
x q[1];
rz(-0.00074978272) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-0.11242871) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6923185) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(-1.7983789) q[1];
rz(-pi) q[2];
rz(-2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(2.6543806) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(0.69865984) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(1.2493856) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65803618) q[0];
sx q[0];
rz(-1.5244966) q[0];
sx q[0];
rz(2.9306843) q[0];
x q[1];
rz(2.2970389) q[2];
sx q[2];
rz(-2.484998) q[2];
sx q[2];
rz(3.055228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82259761) q[1];
sx q[1];
rz(-1.5486451) q[1];
sx q[1];
rz(-1.5763361) q[1];
rz(-pi) q[2];
rz(1.9551679) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(-2.3068908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1406527) q[0];
sx q[0];
rz(-1.5322598) q[0];
sx q[0];
rz(1.7059822) q[0];
rz(-pi) q[1];
rz(0.75002807) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(-2.7761369) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7087047) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(-1.9087285) q[1];
rz(-pi) q[2];
rz(-2.5556373) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(-0.5639329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(-2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(3.1179324) q[0];
rz(0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-0.67217174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-2.0477717) q[0];
x q[1];
rz(-2.6300738) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(2.9715003) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9384267) q[1];
sx q[1];
rz(-2.2397579) q[1];
sx q[1];
rz(2.0069564) q[1];
rz(0.43283312) q[3];
sx q[3];
rz(-1.7107309) q[3];
sx q[3];
rz(0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(0.36995861) q[2];
rz(-1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.9406208) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5621915) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.205668) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(-2.8430812) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
