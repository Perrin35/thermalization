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
rz(0.93267814) q[0];
sx q[0];
rz(-2.7490766) q[0];
sx q[0];
rz(2.0445332) q[0];
rz(1.2708083) q[1];
sx q[1];
rz(-1.1975809) q[1];
sx q[1];
rz(-1.7041915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2559833) q[0];
sx q[0];
rz(-2.3893111) q[0];
sx q[0];
rz(2.6153569) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9095638) q[2];
sx q[2];
rz(-1.4995607) q[2];
sx q[2];
rz(2.247612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7157876) q[1];
sx q[1];
rz(-2.1658848) q[1];
sx q[1];
rz(-2.7462105) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1278856) q[3];
sx q[3];
rz(-2.4936112) q[3];
sx q[3];
rz(1.614384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8561594) q[2];
sx q[2];
rz(-2.9069052) q[2];
sx q[2];
rz(0.50133809) q[2];
rz(1.3242138) q[3];
sx q[3];
rz(-2.5960077) q[3];
sx q[3];
rz(-1.4668303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45601869) q[0];
sx q[0];
rz(-0.40559232) q[0];
sx q[0];
rz(-1.4612041) q[0];
rz(2.9623518) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(0.63757149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6939664) q[0];
sx q[0];
rz(-1.7991095) q[0];
sx q[0];
rz(2.8986321) q[0];
rz(-pi) q[1];
rz(1.7098956) q[2];
sx q[2];
rz(-0.31965986) q[2];
sx q[2];
rz(-0.57397288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7934139) q[1];
sx q[1];
rz(-1.2169588) q[1];
sx q[1];
rz(-0.89576141) q[1];
rz(-pi) q[2];
rz(-2.4781391) q[3];
sx q[3];
rz(-1.2689948) q[3];
sx q[3];
rz(-0.38340195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3011938) q[2];
sx q[2];
rz(-1.0777148) q[2];
sx q[2];
rz(3.0974498) q[2];
rz(2.5282777) q[3];
sx q[3];
rz(-1.8200579) q[3];
sx q[3];
rz(0.80673748) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052385656) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(-1.3432304) q[0];
rz(-1.7228458) q[1];
sx q[1];
rz(-0.90444618) q[1];
sx q[1];
rz(-0.80352965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65261025) q[0];
sx q[0];
rz(-1.5119512) q[0];
sx q[0];
rz(0.68025689) q[0];
rz(-1.1131002) q[2];
sx q[2];
rz(-2.105793) q[2];
sx q[2];
rz(2.0488103) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1413925) q[1];
sx q[1];
rz(-1.8484925) q[1];
sx q[1];
rz(-3.030409) q[1];
rz(-1.9827051) q[3];
sx q[3];
rz(-2.0406699) q[3];
sx q[3];
rz(-1.8793233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.558305) q[2];
sx q[2];
rz(-0.95197695) q[2];
sx q[2];
rz(1.0179016) q[2];
rz(-1.2244276) q[3];
sx q[3];
rz(-1.3320351) q[3];
sx q[3];
rz(1.1289736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99854904) q[0];
sx q[0];
rz(-2.3426549) q[0];
sx q[0];
rz(2.1813188) q[0];
rz(0.1943365) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(-3.1365373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086297) q[0];
sx q[0];
rz(-0.89187183) q[0];
sx q[0];
rz(-0.70569785) q[0];
rz(-pi) q[1];
rz(-1.6089392) q[2];
sx q[2];
rz(-1.7542217) q[2];
sx q[2];
rz(0.18193744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7654842) q[1];
sx q[1];
rz(-1.4579726) q[1];
sx q[1];
rz(2.9214431) q[1];
rz(-pi) q[2];
rz(1.051733) q[3];
sx q[3];
rz(-1.4421652) q[3];
sx q[3];
rz(-2.7277456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2527577) q[2];
sx q[2];
rz(-1.5450666) q[2];
sx q[2];
rz(-1.6952391) q[2];
rz(-1.3085922) q[3];
sx q[3];
rz(-1.7596217) q[3];
sx q[3];
rz(2.5368209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62110353) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(2.5526168) q[0];
rz(-3.1336054) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(-2.0569107) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046291489) q[0];
sx q[0];
rz(-1.3372375) q[0];
sx q[0];
rz(-0.26136847) q[0];
rz(-pi) q[1];
x q[1];
rz(0.029726278) q[2];
sx q[2];
rz(-0.27966248) q[2];
sx q[2];
rz(-0.047911876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1850622) q[1];
sx q[1];
rz(-1.1780199) q[1];
sx q[1];
rz(-2.3529733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067361319) q[3];
sx q[3];
rz(-0.90927059) q[3];
sx q[3];
rz(-2.1811821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40953723) q[2];
sx q[2];
rz(-1.0169225) q[2];
sx q[2];
rz(-1.766905) q[2];
rz(0.095666766) q[3];
sx q[3];
rz(-1.5868264) q[3];
sx q[3];
rz(2.5110551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.10504) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(-1.8010944) q[0];
rz(0.66572491) q[1];
sx q[1];
rz(-1.7981073) q[1];
sx q[1];
rz(2.417477) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3984864) q[0];
sx q[0];
rz(-1.2091276) q[0];
sx q[0];
rz(2.818011) q[0];
rz(2.3409136) q[2];
sx q[2];
rz(-2.4086047) q[2];
sx q[2];
rz(1.6782325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9812079) q[1];
sx q[1];
rz(-1.1120468) q[1];
sx q[1];
rz(-1.7225411) q[1];
rz(-pi) q[2];
rz(-0.0047896623) q[3];
sx q[3];
rz(-1.6671204) q[3];
sx q[3];
rz(2.9504058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.93806997) q[2];
sx q[2];
rz(-1.6911493) q[2];
sx q[2];
rz(-0.098048992) q[2];
rz(2.8010662) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.606474) q[0];
sx q[0];
rz(-2.4381194) q[0];
sx q[0];
rz(-0.055543609) q[0];
rz(-2.7659888) q[1];
sx q[1];
rz(-2.3665078) q[1];
sx q[1];
rz(1.4449878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9268413) q[0];
sx q[0];
rz(-0.80252778) q[0];
sx q[0];
rz(-1.4237798) q[0];
x q[1];
rz(-3.08899) q[2];
sx q[2];
rz(-1.3173283) q[2];
sx q[2];
rz(0.11620493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1174795) q[1];
sx q[1];
rz(-1.584519) q[1];
sx q[1];
rz(0.3776457) q[1];
rz(-2.4395326) q[3];
sx q[3];
rz(-2.5581048) q[3];
sx q[3];
rz(-1.7398011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0982509) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(3.1116389) q[3];
sx q[3];
rz(-1.5104537) q[3];
sx q[3];
rz(2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444262) q[0];
sx q[0];
rz(-0.63969669) q[0];
sx q[0];
rz(-1.2514914) q[0];
rz(-2.331612) q[1];
sx q[1];
rz(-2.4046343) q[1];
sx q[1];
rz(-1.6862148) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43622323) q[0];
sx q[0];
rz(-1.7972662) q[0];
sx q[0];
rz(2.3821077) q[0];
rz(-pi) q[1];
rz(-3.0461629) q[2];
sx q[2];
rz(-2.61728) q[2];
sx q[2];
rz(1.4316259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6941713) q[1];
sx q[1];
rz(-0.93915597) q[1];
sx q[1];
rz(-2.3775435) q[1];
x q[2];
rz(-0.64469962) q[3];
sx q[3];
rz(-2.19826) q[3];
sx q[3];
rz(-0.51034605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.84888419) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.1513618) q[2];
rz(-0.014160841) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(-1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26361156) q[0];
sx q[0];
rz(-2.4284555) q[0];
sx q[0];
rz(-0.94733316) q[0];
rz(0.5961279) q[1];
sx q[1];
rz(-1.7201741) q[1];
sx q[1];
rz(1.5625585) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347711) q[0];
sx q[0];
rz(-1.3993503) q[0];
sx q[0];
rz(-2.417592) q[0];
rz(2.8364592) q[2];
sx q[2];
rz(-2.6465979) q[2];
sx q[2];
rz(2.309805) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5489997) q[1];
sx q[1];
rz(-2.3301417) q[1];
sx q[1];
rz(-2.215272) q[1];
rz(-pi) q[2];
rz(1.3684016) q[3];
sx q[3];
rz(-1.7949275) q[3];
sx q[3];
rz(3.0355767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8922213) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(0.48386827) q[2];
rz(-2.5441235) q[3];
sx q[3];
rz(-1.6941083) q[3];
sx q[3];
rz(2.5634403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396962) q[0];
sx q[0];
rz(-1.551349) q[0];
sx q[0];
rz(-2.7751112) q[0];
rz(-1.81555) q[1];
sx q[1];
rz(-0.96249023) q[1];
sx q[1];
rz(2.0749626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1078306) q[0];
sx q[0];
rz(-2.1559546) q[0];
sx q[0];
rz(2.5385227) q[0];
rz(-pi) q[1];
rz(-1.616912) q[2];
sx q[2];
rz(-0.96509714) q[2];
sx q[2];
rz(-0.92026383) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6910256) q[1];
sx q[1];
rz(-0.98916173) q[1];
sx q[1];
rz(2.8514991) q[1];
rz(1.5068866) q[3];
sx q[3];
rz(-0.47807594) q[3];
sx q[3];
rz(-1.687754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8804973) q[2];
sx q[2];
rz(-0.77793056) q[2];
sx q[2];
rz(1.0658537) q[2];
rz(-0.26782688) q[3];
sx q[3];
rz(-1.4466176) q[3];
sx q[3];
rz(2.6563787) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5814701) q[0];
sx q[0];
rz(-2.1385834) q[0];
sx q[0];
rz(-2.9271097) q[0];
rz(-2.8196234) q[1];
sx q[1];
rz(-1.3170769) q[1];
sx q[1];
rz(-1.9016686) q[1];
rz(-1.5628846) q[2];
sx q[2];
rz(-1.8521761) q[2];
sx q[2];
rz(0.076836486) q[2];
rz(1.1371218) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(1.6655123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
