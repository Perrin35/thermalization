OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(-2.9773657) q[0];
rz(1.7832463) q[2];
sx q[2];
rz(-2.6531086) q[2];
sx q[2];
rz(2.8993895) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2921819) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.8241747) q[1];
rz(-2.2545635) q[3];
sx q[3];
rz(-1.6449882) q[3];
sx q[3];
rz(-2.2825953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59250295) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(-0.94353765) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(-1.1862322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378368) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-0.32365139) q[0];
rz(-1.6329174) q[2];
sx q[2];
rz(-1.3717692) q[2];
sx q[2];
rz(-1.8387427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1520878) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.771404) q[1];
rz(-2.4466189) q[3];
sx q[3];
rz(-1.7193828) q[3];
sx q[3];
rz(-0.23651628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3804669) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911644) q[0];
sx q[0];
rz(-0.69520742) q[0];
sx q[0];
rz(1.3948963) q[0];
rz(1.1782896) q[2];
sx q[2];
rz(-2.6143392) q[2];
sx q[2];
rz(0.96780992) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.324675) q[1];
sx q[1];
rz(-1.9915238) q[1];
sx q[1];
rz(-0.38751985) q[1];
rz(-pi) q[2];
rz(3.0235602) q[3];
sx q[3];
rz(-1.1088088) q[3];
sx q[3];
rz(3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.320425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1962122) q[0];
sx q[0];
rz(-1.5604661) q[0];
sx q[0];
rz(-1.9804079) q[0];
x q[1];
rz(1.3643866) q[2];
sx q[2];
rz(-1.1279391) q[2];
sx q[2];
rz(2.4368311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79975407) q[1];
sx q[1];
rz(-0.8952039) q[1];
sx q[1];
rz(1.231133) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96418013) q[3];
sx q[3];
rz(-2.2607431) q[3];
sx q[3];
rz(2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0358255) q[2];
sx q[2];
rz(-0.92869174) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(-1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.8409761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39869719) q[0];
sx q[0];
rz(-1.5455751) q[0];
sx q[0];
rz(2.7874649) q[0];
rz(-pi) q[1];
rz(-1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(0.14771151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3805441) q[1];
sx q[1];
rz(-0.64861464) q[1];
sx q[1];
rz(-0.56306871) q[1];
rz(-pi) q[2];
rz(-2.8793094) q[3];
sx q[3];
rz(-1.7119006) q[3];
sx q[3];
rz(1.2348246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79648298) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(2.5911962) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4295411) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(2.3343711) q[0];
x q[1];
rz(1.6278218) q[2];
sx q[2];
rz(-0.36643039) q[2];
sx q[2];
rz(-1.0183522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-2.4273708) q[1];
sx q[1];
rz(0.12970129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4828959) q[3];
sx q[3];
rz(-2.4368736) q[3];
sx q[3];
rz(-2.8040915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-0.67725956) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1099694) q[0];
sx q[0];
rz(-0.97951802) q[0];
sx q[0];
rz(1.4800319) q[0];
rz(-pi) q[1];
rz(1.5706967) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(0.11167234) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3358826) q[1];
sx q[1];
rz(-1.6478331) q[1];
sx q[1];
rz(1.9105934) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.8922071) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65803618) q[0];
sx q[0];
rz(-1.617096) q[0];
sx q[0];
rz(2.9306843) q[0];
x q[1];
rz(-1.0480568) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(-2.0975031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5640806) q[1];
sx q[1];
rz(-0.022833303) q[1];
sx q[1];
rz(2.8965685) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8532393) q[3];
sx q[3];
rz(-2.1835612) q[3];
sx q[3];
rz(1.8307277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(-1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1406527) q[0];
sx q[0];
rz(-1.6093328) q[0];
sx q[0];
rz(1.7059822) q[0];
rz(0.55677982) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(2.5533887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4328879) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(1.9087285) q[1];
x q[2];
rz(2.5556373) q[3];
sx q[3];
rz(-1.5065985) q[3];
sx q[3];
rz(2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7260872) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(-2.0477717) q[0];
rz(3.1250481) q[2];
sx q[2];
rz(-0.51157727) q[2];
sx q[2];
rz(1.3862773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9384267) q[1];
sx q[1];
rz(-2.2397579) q[1];
sx q[1];
rz(2.0069564) q[1];
rz(-pi) q[2];
rz(-2.8176114) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(-1.4676263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(-2.771634) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5621915) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(1.9696708) q[2];
sx q[2];
rz(-1.7181859) q[2];
sx q[2];
rz(-2.0167375) q[2];
rz(0.29851144) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
