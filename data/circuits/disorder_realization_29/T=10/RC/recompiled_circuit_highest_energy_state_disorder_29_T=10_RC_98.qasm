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
rz(0.88275498) q[0];
sx q[0];
rz(-2.9967699) q[0];
sx q[0];
rz(0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(-0.96938649) q[1];
sx q[1];
rz(-0.17925395) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0517901) q[0];
sx q[0];
rz(-1.8500916) q[0];
sx q[0];
rz(-2.0816878) q[0];
x q[1];
rz(-0.72291908) q[2];
sx q[2];
rz(-1.5953879) q[2];
sx q[2];
rz(-1.6484156) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7764531) q[1];
sx q[1];
rz(-1.498933) q[1];
sx q[1];
rz(0.30993575) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93710812) q[3];
sx q[3];
rz(-1.5230012) q[3];
sx q[3];
rz(2.5339047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90193191) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(0.88879746) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773847) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(1.8119716) q[0];
rz(-1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764424) q[0];
sx q[0];
rz(-0.2731384) q[0];
sx q[0];
rz(-1.1821724) q[0];
rz(-0.99807605) q[2];
sx q[2];
rz(-2.6055899) q[2];
sx q[2];
rz(1.052945) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0985579) q[1];
sx q[1];
rz(-1.398067) q[1];
sx q[1];
rz(-0.14202001) q[1];
rz(-2.3584189) q[3];
sx q[3];
rz(-2.0370954) q[3];
sx q[3];
rz(-2.9766022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(0.19134276) q[2];
rz(0.30609104) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(-2.7171296) q[0];
rz(1.3407432) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(2.4901966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6474045) q[0];
sx q[0];
rz(-1.5636824) q[0];
sx q[0];
rz(-1.7344385) q[0];
rz(-0.38369757) q[2];
sx q[2];
rz(-1.7406775) q[2];
sx q[2];
rz(3.1184514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5391158) q[1];
sx q[1];
rz(-2.0664363) q[1];
sx q[1];
rz(1.3631166) q[1];
rz(-pi) q[2];
rz(-2.1742854) q[3];
sx q[3];
rz(-0.27316948) q[3];
sx q[3];
rz(-2.4045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.035630781) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(-2.4448709) q[2];
rz(-2.4524955) q[3];
sx q[3];
rz(-1.9870116) q[3];
sx q[3];
rz(0.2564297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.43831393) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(1.0003723) q[0];
rz(-1.4601624) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.2132852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5942237) q[0];
sx q[0];
rz(-1.4495736) q[0];
sx q[0];
rz(-1.8531043) q[0];
x q[1];
rz(1.103066) q[2];
sx q[2];
rz(-2.2437895) q[2];
sx q[2];
rz(1.8813934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.225093) q[1];
sx q[1];
rz(-1.8986393) q[1];
sx q[1];
rz(1.6986548) q[1];
x q[2];
rz(0.090661006) q[3];
sx q[3];
rz(-1.034015) q[3];
sx q[3];
rz(-1.8859552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(-2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475912) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(1.1119615) q[1];
sx q[1];
rz(-1.0589736) q[1];
sx q[1];
rz(1.8066822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955973) q[0];
sx q[0];
rz(-2.3438518) q[0];
sx q[0];
rz(2.255385) q[0];
x q[1];
rz(-1.6755988) q[2];
sx q[2];
rz(-0.99250472) q[2];
sx q[2];
rz(1.0225909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4547382) q[1];
sx q[1];
rz(-1.6461186) q[1];
sx q[1];
rz(-1.7416864) q[1];
rz(-pi) q[2];
rz(1.2615292) q[3];
sx q[3];
rz(-2.5430508) q[3];
sx q[3];
rz(0.33613294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(-2.000957) q[2];
rz(-3.0349777) q[3];
sx q[3];
rz(-1.6116319) q[3];
sx q[3];
rz(2.831736) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(-3.1383681) q[0];
rz(3.0942753) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.3501732) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4989717) q[0];
sx q[0];
rz(-1.3184217) q[0];
sx q[0];
rz(0.060615505) q[0];
x q[1];
rz(-2.240594) q[2];
sx q[2];
rz(-1.0987079) q[2];
sx q[2];
rz(0.53874082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9474831) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(-0.67621381) q[1];
rz(1.7058701) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(-2.6378353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1515767) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(0.081175096) q[2];
rz(-2.2422527) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154295) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(-2.6370866) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-2.1213396) q[1];
sx q[1];
rz(-3.1212433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5581455) q[0];
sx q[0];
rz(-0.37675315) q[0];
sx q[0];
rz(-1.3135629) q[0];
rz(2.11436) q[2];
sx q[2];
rz(-0.24212101) q[2];
sx q[2];
rz(-2.0106135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0416866) q[1];
sx q[1];
rz(-1.3803687) q[1];
sx q[1];
rz(-1.8407525) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4844839) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-2.8750471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67123479) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.7669558) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90877157) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(-2.0594647) q[0];
rz(-2.1807561) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(0.94295162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7060407) q[0];
sx q[0];
rz(-2.0461975) q[0];
sx q[0];
rz(-0.33669223) q[0];
rz(2.9081557) q[2];
sx q[2];
rz(-2.1899411) q[2];
sx q[2];
rz(2.7987628) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3245736) q[1];
sx q[1];
rz(-2.1789411) q[1];
sx q[1];
rz(-2.2961246) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10412962) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(0.96398523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(-1.1478434) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.743478) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107133) q[0];
sx q[0];
rz(-2.78237) q[0];
sx q[0];
rz(3.0391589) q[0];
rz(0.61406413) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(0.817743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2984268) q[0];
sx q[0];
rz(-0.90230391) q[0];
sx q[0];
rz(-2.2737204) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4253584) q[2];
sx q[2];
rz(-0.67611968) q[2];
sx q[2];
rz(-0.506625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9407258) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(-0.84152842) q[1];
rz(-pi) q[2];
rz(0.93877403) q[3];
sx q[3];
rz(-1.1891439) q[3];
sx q[3];
rz(1.3248688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3488591) q[2];
sx q[2];
rz(-2.0642955) q[2];
sx q[2];
rz(-2.8279772) q[2];
rz(-0.46401986) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(-2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713276) q[0];
sx q[0];
rz(-0.79065228) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(2.927921) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(0.22458354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41193) q[0];
sx q[0];
rz(-2.2673635) q[0];
sx q[0];
rz(2.4921472) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9360673) q[2];
sx q[2];
rz(-2.1238616) q[2];
sx q[2];
rz(-1.5129364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1285291) q[1];
sx q[1];
rz(-0.54058248) q[1];
sx q[1];
rz(-2.9180727) q[1];
rz(-2.7216689) q[3];
sx q[3];
rz(-0.89717275) q[3];
sx q[3];
rz(-2.6633584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97682041) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(-2.9898047) q[2];
rz(0.031489059) q[3];
sx q[3];
rz(-1.3969235) q[3];
sx q[3];
rz(-0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0781773) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(2.0582485) q[1];
sx q[1];
rz(-1.5647519) q[1];
sx q[1];
rz(-1.5819989) q[1];
rz(0.2486817) q[2];
sx q[2];
rz(-2.9338825) q[2];
sx q[2];
rz(2.6354229) q[2];
rz(0.29240378) q[3];
sx q[3];
rz(-2.79861) q[3];
sx q[3];
rz(2.5220148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
