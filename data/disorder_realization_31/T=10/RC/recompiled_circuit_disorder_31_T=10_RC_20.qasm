OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52656093) q[0];
sx q[0];
rz(-2.5685413) q[0];
sx q[0];
rz(-0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9561477) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(0.74622112) q[0];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(2.8710499) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1229413) q[1];
sx q[1];
rz(-2.7416347) q[1];
sx q[1];
rz(2.8040228) q[1];
rz(-0.59431608) q[3];
sx q[3];
rz(-1.672097) q[3];
sx q[3];
rz(-3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75630674) q[0];
sx q[0];
rz(-3.1161838) q[0];
sx q[0];
rz(2.4029762) q[0];
x q[1];
rz(-0.39312675) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(-1.7413505) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3719912) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(0.12188697) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96315893) q[3];
sx q[3];
rz(-0.31664407) q[3];
sx q[3];
rz(0.39099993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9804304) q[0];
sx q[0];
rz(-2.3083901) q[0];
sx q[0];
rz(-2.3479793) q[0];
rz(1.0513564) q[2];
sx q[2];
rz(-2.4161077) q[2];
sx q[2];
rz(-1.6348334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6669238) q[1];
sx q[1];
rz(-2.3033934) q[1];
sx q[1];
rz(1.0679507) q[1];
rz(0.46488071) q[3];
sx q[3];
rz(-0.14296963) q[3];
sx q[3];
rz(0.79899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.320257) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757641) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(-0.34253828) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4565059) q[2];
sx q[2];
rz(-1.473046) q[2];
sx q[2];
rz(-3.0435261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0879285) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(-0.019836516) q[1];
x q[2];
rz(-1.4471099) q[3];
sx q[3];
rz(-1.6061155) q[3];
sx q[3];
rz(-0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-2.2713984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4207626) q[0];
sx q[0];
rz(-1.2569068) q[0];
sx q[0];
rz(1.3563966) q[0];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-2.4325271) q[2];
sx q[2];
rz(2.2774334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.732547) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(1.9802666) q[1];
rz(-pi) q[2];
rz(-1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(-0.8997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.4917096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0817889) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(-2.9102737) q[0];
rz(-pi) q[1];
rz(-0.50478023) q[2];
sx q[2];
rz(-2.3961888) q[2];
sx q[2];
rz(-0.2573075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41960934) q[1];
sx q[1];
rz(-1.0227385) q[1];
sx q[1];
rz(-2.4649058) q[1];
x q[2];
rz(3.0400279) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(2.9543119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(-0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-0.033989865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8695553) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(1.3433334) q[0];
rz(-pi) q[1];
rz(-2.2529644) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(-1.467848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1696724) q[1];
sx q[1];
rz(-1.6213413) q[1];
sx q[1];
rz(-0.89269841) q[1];
rz(-pi) q[2];
rz(-0.61492413) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59003419) q[0];
sx q[0];
rz(-2.708205) q[0];
sx q[0];
rz(1.5584459) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3169575) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(3.0658714) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48965028) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(0.55649346) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(-1.8898213) q[0];
rz(-0.66954008) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(0.30977419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(-2.8938328) q[0];
rz(1.3938815) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(-0.44405802) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.579638) q[1];
sx q[1];
rz(-1.5127752) q[1];
sx q[1];
rz(-1.7768363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7731032) q[3];
sx q[3];
rz(-0.62928761) q[3];
sx q[3];
rz(3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(2.4859378) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.9357095) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9046281) q[0];
sx q[0];
rz(-1.4693854) q[0];
sx q[0];
rz(1.2786091) q[0];
rz(-pi) q[1];
rz(3.0707804) q[2];
sx q[2];
rz(-0.76528463) q[2];
sx q[2];
rz(-0.27809696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26350281) q[1];
sx q[1];
rz(-0.51551688) q[1];
sx q[1];
rz(2.0104736) q[1];
rz(3.1390926) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(2.1713184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832227) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-2.0680188) q[2];
sx q[2];
rz(-1.1560658) q[2];
sx q[2];
rz(1.8552468) q[2];
rz(1.0740888) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];