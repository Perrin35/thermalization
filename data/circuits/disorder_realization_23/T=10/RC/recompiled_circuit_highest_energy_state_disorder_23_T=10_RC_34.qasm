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
rz(-2.5459557) q[0];
sx q[0];
rz(-0.41002265) q[0];
sx q[0];
rz(0.79467839) q[0];
rz(-1.0533286) q[1];
sx q[1];
rz(-0.95212189) q[1];
sx q[1];
rz(-1.9917537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4021989) q[0];
sx q[0];
rz(-1.6065242) q[0];
sx q[0];
rz(1.2767919) q[0];
x q[1];
rz(1.8914521) q[2];
sx q[2];
rz(-2.3806664) q[2];
sx q[2];
rz(-0.93283949) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2268035) q[1];
sx q[1];
rz(-2.7046596) q[1];
sx q[1];
rz(0.41586693) q[1];
rz(-pi) q[2];
rz(0.46709316) q[3];
sx q[3];
rz(-1.8008999) q[3];
sx q[3];
rz(-2.439374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4475693) q[2];
sx q[2];
rz(-0.99267107) q[2];
sx q[2];
rz(-0.64219323) q[2];
rz(2.0726974) q[3];
sx q[3];
rz(-1.5725719) q[3];
sx q[3];
rz(1.4890891) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.883413) q[0];
sx q[0];
rz(-0.0029819948) q[0];
sx q[0];
rz(0.51134837) q[0];
rz(1.6343575) q[1];
sx q[1];
rz(-1.2751445) q[1];
sx q[1];
rz(1.3828329) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16973142) q[0];
sx q[0];
rz(-1.3644553) q[0];
sx q[0];
rz(-1.9623164) q[0];
rz(-pi) q[1];
rz(-0.34592705) q[2];
sx q[2];
rz(-0.9285766) q[2];
sx q[2];
rz(1.3042892) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1369776) q[1];
sx q[1];
rz(-1.8796693) q[1];
sx q[1];
rz(-1.1501794) q[1];
rz(-pi) q[2];
rz(0.17937998) q[3];
sx q[3];
rz(-2.6112643) q[3];
sx q[3];
rz(1.3190003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.6766659) q[2];
sx q[2];
rz(-1.9366465) q[2];
sx q[2];
rz(-3.0386772) q[2];
rz(-2.0788705) q[3];
sx q[3];
rz(-1.9461742) q[3];
sx q[3];
rz(-0.11500558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677143) q[0];
sx q[0];
rz(-1.0391087) q[0];
sx q[0];
rz(0.1097196) q[0];
rz(0.32490718) q[1];
sx q[1];
rz(-1.9903245) q[1];
sx q[1];
rz(2.6085764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4895297) q[0];
sx q[0];
rz(-2.8127413) q[0];
sx q[0];
rz(2.7016599) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.052003666) q[2];
sx q[2];
rz(-1.5947761) q[2];
sx q[2];
rz(-2.5131186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17688454) q[1];
sx q[1];
rz(-2.1358651) q[1];
sx q[1];
rz(-2.3749897) q[1];
rz(-pi) q[2];
rz(0.43709938) q[3];
sx q[3];
rz(-2.282939) q[3];
sx q[3];
rz(-2.8959078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.816232) q[2];
sx q[2];
rz(-1.8133138) q[2];
sx q[2];
rz(1.4556966) q[2];
rz(0.085518941) q[3];
sx q[3];
rz(-2.072008) q[3];
sx q[3];
rz(0.94676179) q[3];
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
rz(-2.7607255) q[0];
sx q[0];
rz(-0.65422288) q[0];
sx q[0];
rz(2.4192659) q[0];
rz(1.5707312) q[1];
sx q[1];
rz(-1.1474362) q[1];
sx q[1];
rz(3.0991203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4054779) q[0];
sx q[0];
rz(-1.7543989) q[0];
sx q[0];
rz(0.45001375) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4804777) q[2];
sx q[2];
rz(-1.5549763) q[2];
sx q[2];
rz(2.9243369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60085697) q[1];
sx q[1];
rz(-1.7607948) q[1];
sx q[1];
rz(0.79278391) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3298395) q[3];
sx q[3];
rz(-0.78380871) q[3];
sx q[3];
rz(-2.2026521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2461569) q[2];
sx q[2];
rz(-1.7652721) q[2];
sx q[2];
rz(-2.8975272) q[2];
rz(0.79143381) q[3];
sx q[3];
rz(-1.3118728) q[3];
sx q[3];
rz(0.70146504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0303845) q[0];
sx q[0];
rz(-3.0272439) q[0];
sx q[0];
rz(2.9682888) q[0];
rz(2.2068842) q[1];
sx q[1];
rz(-0.84988958) q[1];
sx q[1];
rz(0.25397837) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5565709) q[0];
sx q[0];
rz(-2.3818172) q[0];
sx q[0];
rz(1.4069446) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2082269) q[2];
sx q[2];
rz(-2.1340279) q[2];
sx q[2];
rz(-2.8993487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9137302) q[1];
sx q[1];
rz(-0.87764064) q[1];
sx q[1];
rz(-1.5437897) q[1];
rz(-1.1831102) q[3];
sx q[3];
rz(-1.630894) q[3];
sx q[3];
rz(-0.60754402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0737334) q[2];
sx q[2];
rz(-2.8974055) q[2];
sx q[2];
rz(-1.4402639) q[2];
rz(-2.4654147) q[3];
sx q[3];
rz(-1.22236) q[3];
sx q[3];
rz(-3.0618099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.9192231) q[0];
sx q[0];
rz(-1.3105404) q[0];
sx q[0];
rz(0.72810143) q[0];
rz(-1.9758196) q[1];
sx q[1];
rz(-2.246942) q[1];
sx q[1];
rz(2.5808835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8090208) q[0];
sx q[0];
rz(-0.14223465) q[0];
sx q[0];
rz(0.47551544) q[0];
rz(-0.54975551) q[2];
sx q[2];
rz(-1.4696331) q[2];
sx q[2];
rz(1.8515974) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.06212) q[1];
sx q[1];
rz(-1.9626856) q[1];
sx q[1];
rz(0.23418871) q[1];
x q[2];
rz(-1.6339763) q[3];
sx q[3];
rz(-1.6681021) q[3];
sx q[3];
rz(-2.7246876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59549436) q[2];
sx q[2];
rz(-1.4608773) q[2];
sx q[2];
rz(2.0687436) q[2];
rz(-2.8192375) q[3];
sx q[3];
rz(-2.7287546) q[3];
sx q[3];
rz(0.391092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.76674616) q[0];
sx q[0];
rz(-1.640919) q[0];
sx q[0];
rz(0.40819502) q[0];
rz(0.32632581) q[1];
sx q[1];
rz(-2.7944481) q[1];
sx q[1];
rz(-1.4761285) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2583019) q[0];
sx q[0];
rz(-1.6101478) q[0];
sx q[0];
rz(1.395306) q[0];
rz(-pi) q[1];
rz(-1.9397975) q[2];
sx q[2];
rz(-0.59518669) q[2];
sx q[2];
rz(-0.068142224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1211675) q[1];
sx q[1];
rz(-1.5161991) q[1];
sx q[1];
rz(1.4538641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.514634) q[3];
sx q[3];
rz(-2.857465) q[3];
sx q[3];
rz(1.8654902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8866715) q[2];
sx q[2];
rz(-1.6559867) q[2];
sx q[2];
rz(0.86611789) q[2];
rz(2.4050889) q[3];
sx q[3];
rz(-1.085142) q[3];
sx q[3];
rz(0.88732639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0199468) q[0];
sx q[0];
rz(-3.096088) q[0];
sx q[0];
rz(-1.2787) q[0];
rz(2.1389029) q[1];
sx q[1];
rz(-1.8808782) q[1];
sx q[1];
rz(-0.444828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76679517) q[0];
sx q[0];
rz(-2.594769) q[0];
sx q[0];
rz(1.8575791) q[0];
rz(-pi) q[1];
rz(2.629822) q[2];
sx q[2];
rz(-0.79088849) q[2];
sx q[2];
rz(-0.98423401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7463716) q[1];
sx q[1];
rz(-2.3276797) q[1];
sx q[1];
rz(-0.16415747) q[1];
rz(-pi) q[2];
rz(-0.37737198) q[3];
sx q[3];
rz(-1.3283973) q[3];
sx q[3];
rz(0.19211543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43850809) q[2];
sx q[2];
rz(-1.4139516) q[2];
sx q[2];
rz(-0.50602305) q[2];
rz(0.94432008) q[3];
sx q[3];
rz(-3.1157065) q[3];
sx q[3];
rz(-0.73616141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29189062) q[0];
sx q[0];
rz(-0.99221748) q[0];
sx q[0];
rz(-0.010183656) q[0];
rz(-0.61095515) q[1];
sx q[1];
rz(-2.1612031) q[1];
sx q[1];
rz(-0.082286509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1707538) q[0];
sx q[0];
rz(-2.0505095) q[0];
sx q[0];
rz(-1.1520349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9111229) q[2];
sx q[2];
rz(-2.2374212) q[2];
sx q[2];
rz(2.8067592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6895143) q[1];
sx q[1];
rz(-1.6321215) q[1];
sx q[1];
rz(-0.21558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.786834) q[3];
sx q[3];
rz(-2.2240765) q[3];
sx q[3];
rz(0.21524425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33132195) q[2];
sx q[2];
rz(-1.4484582) q[2];
sx q[2];
rz(-0.70915478) q[2];
rz(-1.6592615) q[3];
sx q[3];
rz(-2.5100561) q[3];
sx q[3];
rz(0.46019301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3191147) q[0];
sx q[0];
rz(-2.3111486) q[0];
sx q[0];
rz(2.5122232) q[0];
rz(-1.4156226) q[1];
sx q[1];
rz(-1.4868088) q[1];
sx q[1];
rz(1.9633044) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45635763) q[0];
sx q[0];
rz(-0.97536063) q[0];
sx q[0];
rz(2.6490982) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35188108) q[2];
sx q[2];
rz(-2.0469672) q[2];
sx q[2];
rz(1.1426403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7892742) q[1];
sx q[1];
rz(-1.4613918) q[1];
sx q[1];
rz(2.3477702) q[1];
x q[2];
rz(2.1363821) q[3];
sx q[3];
rz(-0.93200383) q[3];
sx q[3];
rz(1.8106724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87469953) q[2];
sx q[2];
rz(-1.5723672) q[2];
sx q[2];
rz(-2.3229522) q[2];
rz(-1.6308174) q[3];
sx q[3];
rz(-1.2997593) q[3];
sx q[3];
rz(0.40217933) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4928987) q[0];
sx q[0];
rz(-2.0317827) q[0];
sx q[0];
rz(-1.6775525) q[0];
rz(1.2661487) q[1];
sx q[1];
rz(-1.1504953) q[1];
sx q[1];
rz(-1.6082416) q[1];
rz(-0.89546236) q[2];
sx q[2];
rz(-1.3668899) q[2];
sx q[2];
rz(-1.5820506) q[2];
rz(1.6615909) q[3];
sx q[3];
rz(-1.3868854) q[3];
sx q[3];
rz(-1.8793061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
