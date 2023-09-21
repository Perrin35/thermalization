OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(-2.2111501) q[0];
rz(-2.5118877) q[1];
sx q[1];
rz(-1.1344818) q[1];
sx q[1];
rz(1.1073444) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168092) q[0];
sx q[0];
rz(-0.9073572) q[0];
sx q[0];
rz(-1.1532564) q[0];
x q[1];
rz(-0.8503352) q[2];
sx q[2];
rz(-1.9413661) q[2];
sx q[2];
rz(2.4239899) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2335637) q[1];
sx q[1];
rz(-1.4587914) q[1];
sx q[1];
rz(-1.8506552) q[1];
rz(-pi) q[2];
rz(-2.2060478) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(-2.2176544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(-1.8135653) q[2];
rz(2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-3.0096171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593453) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(-0.88062084) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-0.41795119) q[1];
sx q[1];
rz(2.3243288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5215223) q[0];
sx q[0];
rz(-1.2261454) q[0];
sx q[0];
rz(0.14981139) q[0];
rz(-1.0802644) q[2];
sx q[2];
rz(-1.763952) q[2];
sx q[2];
rz(0.9888538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9643147) q[1];
sx q[1];
rz(-0.69523584) q[1];
sx q[1];
rz(-2.7362105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2137787) q[3];
sx q[3];
rz(-2.6414053) q[3];
sx q[3];
rz(0.57750765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-0.36402738) q[2];
rz(2.1552127) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(-1.6769489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36689511) q[0];
sx q[0];
rz(-0.8323454) q[0];
sx q[0];
rz(0.96631518) q[0];
rz(-0.19293383) q[1];
sx q[1];
rz(-2.0529592) q[1];
sx q[1];
rz(-1.4470709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149949) q[0];
sx q[0];
rz(-1.5532171) q[0];
sx q[0];
rz(-0.21501644) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(1.9384055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6009532) q[1];
sx q[1];
rz(-0.97266957) q[1];
sx q[1];
rz(-1.1126493) q[1];
rz(-pi) q[2];
rz(-1.0560016) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(1.521829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7029999) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(2.036371) q[2];
rz(-2.3953719) q[3];
sx q[3];
rz(-1.6522224) q[3];
sx q[3];
rz(0.97572774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063909) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(-1.4659457) q[0];
rz(-2.8566467) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(2.7526061) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3272414) q[0];
sx q[0];
rz(-1.2885433) q[0];
sx q[0];
rz(1.3319356) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0573425) q[2];
sx q[2];
rz(-1.3975189) q[2];
sx q[2];
rz(-1.372352) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.092408808) q[1];
sx q[1];
rz(-2.1158943) q[1];
sx q[1];
rz(-1.070302) q[1];
x q[2];
rz(0.22104927) q[3];
sx q[3];
rz(-0.76023686) q[3];
sx q[3];
rz(1.5787214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4185562) q[2];
sx q[2];
rz(-1.1636461) q[2];
sx q[2];
rz(2.6848865) q[2];
rz(1.6263973) q[3];
sx q[3];
rz(-0.94907343) q[3];
sx q[3];
rz(0.28234282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221955) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(2.7815681) q[0];
rz(2.4941764) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(0.49450758) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9279328) q[0];
sx q[0];
rz(-1.4466009) q[0];
sx q[0];
rz(0.28161063) q[0];
rz(-0.039426609) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(-2.6089422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8066986) q[1];
sx q[1];
rz(-0.94167275) q[1];
sx q[1];
rz(-1.4966399) q[1];
x q[2];
rz(2.6974929) q[3];
sx q[3];
rz(-1.5478304) q[3];
sx q[3];
rz(-2.3820153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4313724) q[2];
sx q[2];
rz(-2.4857095) q[2];
sx q[2];
rz(2.8430856) q[2];
rz(2.8295637) q[3];
sx q[3];
rz(-1.8108862) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1453778) q[0];
sx q[0];
rz(-0.72308102) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(3.1205102) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(1.9063937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.9039246) q[0];
sx q[0];
rz(1.8963277) q[0];
rz(-pi) q[1];
rz(2.4371106) q[2];
sx q[2];
rz(-1.8330169) q[2];
sx q[2];
rz(-0.30202497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8630353) q[1];
sx q[1];
rz(-0.99579358) q[1];
sx q[1];
rz(2.0726191) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33552334) q[3];
sx q[3];
rz(-0.24970679) q[3];
sx q[3];
rz(-1.6186796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6482676) q[2];
sx q[2];
rz(-0.92612925) q[2];
sx q[2];
rz(2.7098999) q[2];
rz(1.7290944) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39111185) q[0];
sx q[0];
rz(-1.1688122) q[0];
sx q[0];
rz(-0.11211638) q[0];
rz(-0.21513367) q[1];
sx q[1];
rz(-1.5605749) q[1];
sx q[1];
rz(2.0281866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8002692) q[0];
sx q[0];
rz(-1.1870664) q[0];
sx q[0];
rz(1.8510438) q[0];
rz(1.3117123) q[2];
sx q[2];
rz(-1.8763181) q[2];
sx q[2];
rz(-0.53369001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6326633) q[1];
sx q[1];
rz(-0.76970657) q[1];
sx q[1];
rz(2.6973004) q[1];
rz(1.1235085) q[3];
sx q[3];
rz(-1.7637858) q[3];
sx q[3];
rz(-2.6393294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-2.4089456) q[2];
sx q[2];
rz(3.0155638) q[2];
rz(-2.0942988) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(-0.70820156) q[3];
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
rz(0.66013181) q[0];
sx q[0];
rz(-2.3829057) q[0];
sx q[0];
rz(1.6814167) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-1.0943202) q[1];
sx q[1];
rz(-1.2089027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9019933) q[0];
sx q[0];
rz(-1.308846) q[0];
sx q[0];
rz(0.082017935) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2751341) q[2];
sx q[2];
rz(-0.87477113) q[2];
sx q[2];
rz(-2.8528086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9637451) q[1];
sx q[1];
rz(-2.3235465) q[1];
sx q[1];
rz(2.6820116) q[1];
rz(-pi) q[2];
rz(2.9221411) q[3];
sx q[3];
rz(-2.1900574) q[3];
sx q[3];
rz(-1.0964364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7631491) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(1.8219927) q[2];
rz(-0.59213263) q[3];
sx q[3];
rz(-1.7254646) q[3];
sx q[3];
rz(-3.106451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(1.7804902) q[0];
rz(1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-0.39168721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254588) q[0];
sx q[0];
rz(-0.43457169) q[0];
sx q[0];
rz(0.01390121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7462256) q[2];
sx q[2];
rz(-0.48465109) q[2];
sx q[2];
rz(-1.6389099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3561331) q[1];
sx q[1];
rz(-1.2411989) q[1];
sx q[1];
rz(2.3258924) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12556062) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.7988127) q[2];
sx q[2];
rz(-1.8048145) q[2];
rz(-2.3796066) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(-0.80037642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.8005463) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(-1.4292498) q[1];
sx q[1];
rz(-1.0639023) q[1];
sx q[1];
rz(0.16194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8853332) q[0];
sx q[0];
rz(-0.67978871) q[0];
sx q[0];
rz(-1.1023561) q[0];
x q[1];
rz(1.6185206) q[2];
sx q[2];
rz(-2.952791) q[2];
sx q[2];
rz(1.7392841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.098617741) q[1];
sx q[1];
rz(-2.9264755) q[1];
sx q[1];
rz(2.2153562) q[1];
rz(0.52311388) q[3];
sx q[3];
rz(-1.599708) q[3];
sx q[3];
rz(-1.6541964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95742115) q[2];
sx q[2];
rz(-1.1051757) q[2];
sx q[2];
rz(1.7133678) q[2];
rz(1.1994294) q[3];
sx q[3];
rz(-2.0482443) q[3];
sx q[3];
rz(1.3706346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409055) q[0];
sx q[0];
rz(-0.4754684) q[0];
sx q[0];
rz(-1.1204002) q[0];
rz(1.7715001) q[1];
sx q[1];
rz(-0.94540989) q[1];
sx q[1];
rz(2.170845) q[1];
rz(-0.8926819) q[2];
sx q[2];
rz(-0.18845367) q[2];
sx q[2];
rz(-1.0343196) q[2];
rz(-0.45392848) q[3];
sx q[3];
rz(-2.6700927) q[3];
sx q[3];
rz(2.3910458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];