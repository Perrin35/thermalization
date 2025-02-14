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
rz(-2.3075624) q[0];
sx q[0];
rz(2.7355255) q[0];
sx q[0];
rz(14.526215) q[0];
rz(1.740068) q[1];
sx q[1];
rz(5.8086173) q[1];
sx q[1];
rz(9.557815) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8048153) q[0];
sx q[0];
rz(-0.93685564) q[0];
sx q[0];
rz(0.32685771) q[0];
rz(-pi) q[1];
rz(-0.51271397) q[2];
sx q[2];
rz(-1.1774592) q[2];
sx q[2];
rz(-2.8296997) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.052357167) q[1];
sx q[1];
rz(-1.7947788) q[1];
sx q[1];
rz(1.8941982) q[1];
rz(-pi) q[2];
rz(-0.36433371) q[3];
sx q[3];
rz(-1.0291463) q[3];
sx q[3];
rz(0.9389824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0066321) q[2];
sx q[2];
rz(-1.0198318) q[2];
sx q[2];
rz(-2.4572065) q[2];
rz(1.1413752) q[3];
sx q[3];
rz(-2.8751774) q[3];
sx q[3];
rz(-0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77966493) q[0];
sx q[0];
rz(-0.68988887) q[0];
sx q[0];
rz(-0.46098125) q[0];
rz(1.1191248) q[1];
sx q[1];
rz(-0.42367595) q[1];
sx q[1];
rz(-2.8295595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0756931) q[0];
sx q[0];
rz(-1.5585494) q[0];
sx q[0];
rz(3.1135983) q[0];
rz(-pi) q[1];
rz(0.80993311) q[2];
sx q[2];
rz(-1.6259127) q[2];
sx q[2];
rz(0.3050366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4496135) q[1];
sx q[1];
rz(-0.17573729) q[1];
sx q[1];
rz(-1.537561) q[1];
rz(-pi) q[2];
rz(-1.9309405) q[3];
sx q[3];
rz(-0.68406287) q[3];
sx q[3];
rz(1.6863914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1877039) q[2];
sx q[2];
rz(-1.1053332) q[2];
sx q[2];
rz(0.40404955) q[2];
rz(-0.36014253) q[3];
sx q[3];
rz(-1.4934243) q[3];
sx q[3];
rz(-0.68296105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78524154) q[0];
sx q[0];
rz(-1.0026362) q[0];
sx q[0];
rz(-0.32401618) q[0];
rz(1.0600545) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(0.70045984) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4361079) q[0];
sx q[0];
rz(-1.7388845) q[0];
sx q[0];
rz(-1.7343434) q[0];
x q[1];
rz(-2.926658) q[2];
sx q[2];
rz(-1.9086468) q[2];
sx q[2];
rz(2.1825254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8328402) q[1];
sx q[1];
rz(-2.3165417) q[1];
sx q[1];
rz(-1.8925335) q[1];
rz(0.19566124) q[3];
sx q[3];
rz(-1.1563468) q[3];
sx q[3];
rz(2.5959542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6274274) q[2];
sx q[2];
rz(-2.8943987) q[2];
sx q[2];
rz(-0.79666454) q[2];
rz(-1.1967777) q[3];
sx q[3];
rz(-1.5408885) q[3];
sx q[3];
rz(0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45348039) q[0];
sx q[0];
rz(-2.1402833) q[0];
sx q[0];
rz(-1.3230811) q[0];
rz(0.46609136) q[1];
sx q[1];
rz(-2.2345462) q[1];
sx q[1];
rz(-1.9922493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7931343) q[0];
sx q[0];
rz(-2.0597947) q[0];
sx q[0];
rz(-0.49599802) q[0];
rz(2.0987112) q[2];
sx q[2];
rz(-0.68860523) q[2];
sx q[2];
rz(-1.7797888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9117078) q[1];
sx q[1];
rz(-2.9344892) q[1];
sx q[1];
rz(-2.9207499) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2364616) q[3];
sx q[3];
rz(-1.4282246) q[3];
sx q[3];
rz(-0.73745525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65679067) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(2.316324) q[2];
rz(1.2800062) q[3];
sx q[3];
rz(-1.5214336) q[3];
sx q[3];
rz(-0.72117225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44329062) q[0];
sx q[0];
rz(-1.2606786) q[0];
sx q[0];
rz(1.3662421) q[0];
rz(1.0900991) q[1];
sx q[1];
rz(-1.6366199) q[1];
sx q[1];
rz(1.3390138) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96025673) q[0];
sx q[0];
rz(-0.68013817) q[0];
sx q[0];
rz(-2.8614317) q[0];
x q[1];
rz(2.5191804) q[2];
sx q[2];
rz(-2.174985) q[2];
sx q[2];
rz(-2.6780918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4046062) q[1];
sx q[1];
rz(-1.4530208) q[1];
sx q[1];
rz(2.9280258) q[1];
rz(-pi) q[2];
rz(-2.3424266) q[3];
sx q[3];
rz(-0.46119565) q[3];
sx q[3];
rz(-2.9950664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1947386) q[2];
sx q[2];
rz(-1.4118492) q[2];
sx q[2];
rz(2.8080688) q[2];
rz(-2.5528095) q[3];
sx q[3];
rz(-0.47313658) q[3];
sx q[3];
rz(-0.80140448) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81353417) q[0];
sx q[0];
rz(-1.7840339) q[0];
sx q[0];
rz(-1.300977) q[0];
rz(1.0282372) q[1];
sx q[1];
rz(-2.6515549) q[1];
sx q[1];
rz(-2.6992544) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8632354) q[0];
sx q[0];
rz(-1.8148737) q[0];
sx q[0];
rz(-2.6574479) q[0];
rz(-pi) q[1];
rz(-0.84965923) q[2];
sx q[2];
rz(-2.4772212) q[2];
sx q[2];
rz(-2.0400782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3446253) q[1];
sx q[1];
rz(-1.7360002) q[1];
sx q[1];
rz(-2.7112673) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4668999) q[3];
sx q[3];
rz(-1.8891617) q[3];
sx q[3];
rz(-1.1863134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(0.4717353) q[2];
rz(1.6996023) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(0.26427463) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1137375) q[0];
sx q[0];
rz(-1.6738482) q[0];
sx q[0];
rz(-0.16673985) q[0];
rz(-0.79404229) q[1];
sx q[1];
rz(-1.200095) q[1];
sx q[1];
rz(0.099743191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4605527) q[0];
sx q[0];
rz(-1.0175219) q[0];
sx q[0];
rz(2.7727292) q[0];
rz(-pi) q[1];
rz(-2.7664423) q[2];
sx q[2];
rz(-0.98723961) q[2];
sx q[2];
rz(0.75022682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.197021) q[1];
sx q[1];
rz(-0.58496956) q[1];
sx q[1];
rz(2.9119063) q[1];
rz(-pi) q[2];
rz(2.7230303) q[3];
sx q[3];
rz(-2.1893614) q[3];
sx q[3];
rz(1.5582066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86228592) q[2];
sx q[2];
rz(-2.0760832) q[2];
sx q[2];
rz(-0.12969895) q[2];
rz(-1.2636403) q[3];
sx q[3];
rz(-1.940515) q[3];
sx q[3];
rz(2.9960846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97067863) q[0];
sx q[0];
rz(-0.92249528) q[0];
sx q[0];
rz(0.42651919) q[0];
rz(-1.0921987) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(-1.0027286) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0195983) q[0];
sx q[0];
rz(-1.3126123) q[0];
sx q[0];
rz(-0.95457558) q[0];
rz(0.72288798) q[2];
sx q[2];
rz(-1.5010415) q[2];
sx q[2];
rz(-0.43089128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5432524) q[1];
sx q[1];
rz(-1.9045254) q[1];
sx q[1];
rz(2.0327492) q[1];
rz(-pi) q[2];
rz(-0.2360446) q[3];
sx q[3];
rz(-1.1268643) q[3];
sx q[3];
rz(0.027396552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2670474) q[2];
sx q[2];
rz(-0.55570498) q[2];
sx q[2];
rz(3.0181273) q[2];
rz(-0.63001436) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(0.71648487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89710871) q[0];
sx q[0];
rz(-2.4429595) q[0];
sx q[0];
rz(-0.83936349) q[0];
rz(0.15946236) q[1];
sx q[1];
rz(-0.99226743) q[1];
sx q[1];
rz(-0.92963591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146067) q[0];
sx q[0];
rz(-2.2898798) q[0];
sx q[0];
rz(1.9107463) q[0];
rz(-0.18126296) q[2];
sx q[2];
rz(-2.5428704) q[2];
sx q[2];
rz(-1.8229654) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60022053) q[1];
sx q[1];
rz(-1.6959136) q[1];
sx q[1];
rz(2.255463) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93299039) q[3];
sx q[3];
rz(-1.3853042) q[3];
sx q[3];
rz(0.18292566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9440426) q[2];
sx q[2];
rz(-1.1669179) q[2];
sx q[2];
rz(2.3738532) q[2];
rz(0.36457101) q[3];
sx q[3];
rz(-1.7578099) q[3];
sx q[3];
rz(3.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0070873) q[0];
sx q[0];
rz(-2.5275079) q[0];
sx q[0];
rz(-0.8836723) q[0];
rz(-0.20835749) q[1];
sx q[1];
rz(-2.6388984) q[1];
sx q[1];
rz(1.3349894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9656031) q[0];
sx q[0];
rz(-1.635968) q[0];
sx q[0];
rz(-1.6287281) q[0];
rz(-pi) q[1];
rz(-0.31317385) q[2];
sx q[2];
rz(-0.54944456) q[2];
sx q[2];
rz(-0.012481364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22182759) q[1];
sx q[1];
rz(-1.034267) q[1];
sx q[1];
rz(2.6229658) q[1];
rz(-pi) q[2];
rz(0.86746704) q[3];
sx q[3];
rz(-2.2604979) q[3];
sx q[3];
rz(2.3610601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9779382) q[2];
sx q[2];
rz(-1.6341219) q[2];
sx q[2];
rz(3.0955691) q[2];
rz(0.16511354) q[3];
sx q[3];
rz(-2.8486227) q[3];
sx q[3];
rz(-1.9273812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5788427) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(2.9764755) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(-0.44331917) q[2];
sx q[2];
rz(-1.8675315) q[2];
sx q[2];
rz(0.055589614) q[2];
rz(1.0551999) q[3];
sx q[3];
rz(-1.996481) q[3];
sx q[3];
rz(0.79660637) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
