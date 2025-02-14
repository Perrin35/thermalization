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
rz(3.0796666) q[0];
sx q[0];
rz(-2.5653337) q[0];
sx q[0];
rz(1.040333) q[0];
rz(-pi) q[1];
rz(2.4186736) q[2];
sx q[2];
rz(-1.5462048) q[2];
sx q[2];
rz(-1.493177) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36513953) q[1];
sx q[1];
rz(-1.498933) q[1];
sx q[1];
rz(0.30993575) q[1];
rz(0.93710812) q[3];
sx q[3];
rz(-1.5230012) q[3];
sx q[3];
rz(-2.5339047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90193191) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(2.530976) q[2];
rz(2.2527952) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(2.4077267) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773847) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.3296211) q[0];
rz(1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(2.2775473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9784911) q[0];
sx q[0];
rz(-1.8231043) q[0];
sx q[0];
rz(3.0358394) q[0];
x q[1];
rz(-2.8301622) q[2];
sx q[2];
rz(-2.0144147) q[2];
sx q[2];
rz(-1.6964427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.043034779) q[1];
sx q[1];
rz(-1.7435257) q[1];
sx q[1];
rz(-0.14202001) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95324272) q[3];
sx q[3];
rz(-2.2526178) q[3];
sx q[3];
rz(2.1566331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(-2.9502499) q[2];
rz(-2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3092344) q[0];
sx q[0];
rz(-1.2795871) q[0];
sx q[0];
rz(-2.7171296) q[0];
rz(1.8008495) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(0.65139604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021928) q[0];
sx q[0];
rz(-0.16379539) q[0];
sx q[0];
rz(1.6144362) q[0];
rz(-pi) q[1];
rz(2.7578951) q[2];
sx q[2];
rz(-1.4009152) q[2];
sx q[2];
rz(0.023141247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5391158) q[1];
sx q[1];
rz(-2.0664363) q[1];
sx q[1];
rz(1.3631166) q[1];
rz(-2.983909) q[3];
sx q[3];
rz(-1.7947949) q[3];
sx q[3];
rz(-1.3583314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1059619) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(2.4448709) q[2];
rz(2.4524955) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7032787) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(-2.1412204) q[0];
rz(-1.6814303) q[1];
sx q[1];
rz(-1.5337475) q[1];
sx q[1];
rz(1.2132852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1532261) q[0];
sx q[0];
rz(-1.2906162) q[0];
sx q[0];
rz(-3.0154254) q[0];
x q[1];
rz(0.51474173) q[2];
sx q[2];
rz(-0.79823433) q[2];
sx q[2];
rz(1.9412184) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4459232) q[1];
sx q[1];
rz(-1.4497821) q[1];
sx q[1];
rz(-0.33034973) q[1];
rz(-3.0509316) q[3];
sx q[3];
rz(-1.034015) q[3];
sx q[3];
rz(-1.8859552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(0.064112045) q[2];
rz(2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(0.39976111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.0475912) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(0.79175788) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9324786) q[0];
sx q[0];
rz(-0.9830342) q[0];
sx q[0];
rz(-2.5665119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58081268) q[2];
sx q[2];
rz(-1.6585095) q[2];
sx q[2];
rz(2.650819) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68685442) q[1];
sx q[1];
rz(-1.4954741) q[1];
sx q[1];
rz(1.3999062) q[1];
x q[2];
rz(0.20466699) q[3];
sx q[3];
rz(-1.0042666) q[3];
sx q[3];
rz(0.03290225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(2.000957) q[2];
rz(-0.10661495) q[3];
sx q[3];
rz(-1.6116319) q[3];
sx q[3];
rz(-2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8150197) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(0.003224592) q[0];
rz(3.0942753) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.3501732) q[1];
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
rz(-2.0428847) q[2];
sx q[2];
rz(2.6028518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9474831) q[1];
sx q[1];
rz(-2.1305741) q[1];
sx q[1];
rz(-0.67621381) q[1];
rz(1.7914823) q[3];
sx q[3];
rz(-1.5410564) q[3];
sx q[3];
rz(1.9427751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99001592) q[2];
sx q[2];
rz(-2.9559957) q[2];
sx q[2];
rz(-3.0604176) q[2];
rz(-2.2422527) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154295) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(0.02034932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.394144) q[0];
sx q[0];
rz(-1.4770623) q[0];
sx q[0];
rz(-1.205349) q[0];
rz(1.3624914) q[2];
sx q[2];
rz(-1.4464738) q[2];
sx q[2];
rz(0.97036874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0416866) q[1];
sx q[1];
rz(-1.3803687) q[1];
sx q[1];
rz(1.8407525) q[1];
rz(1.4844839) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-0.26654551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.3746369) q[2];
rz(-1.5291519) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(-3.0488739) q[3];
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
rz(-0.90877157) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(-2.198641) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1652007) q[0];
sx q[0];
rz(-1.2726674) q[0];
sx q[0];
rz(1.0714897) q[0];
rz(-pi) q[1];
rz(-1.2569179) q[2];
sx q[2];
rz(-2.4853443) q[2];
sx q[2];
rz(-0.73168025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3228938) q[1];
sx q[1];
rz(-0.9095236) q[1];
sx q[1];
rz(-0.7612919) q[1];
rz(-1.8887159) q[3];
sx q[3];
rz(-2.8189481) q[3];
sx q[3];
rz(-1.8436197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(-1.9937493) q[2];
rz(-0.36192274) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(1.743478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73087937) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(0.10243375) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(2.3238497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9303904) q[0];
sx q[0];
rz(-2.1028535) q[0];
sx q[0];
rz(-0.802687) q[0];
x q[1];
rz(-1.7162343) q[2];
sx q[2];
rz(-2.465473) q[2];
sx q[2];
rz(-0.506625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2008668) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(-2.3000642) q[1];
rz(-pi) q[2];
rz(-0.97408143) q[3];
sx q[3];
rz(-0.72457216) q[3];
sx q[3];
rz(-2.4251079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(-2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(-2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87026507) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(2.927921) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(-2.9170091) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29447039) q[0];
sx q[0];
rz(-2.0531539) q[0];
sx q[0];
rz(2.3806974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5245536) q[2];
sx q[2];
rz(-2.4894425) q[2];
sx q[2];
rz(-0.99936501) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27222363) q[1];
sx q[1];
rz(-2.096513) q[1];
sx q[1];
rz(-1.7030667) q[1];
x q[2];
rz(-2.0430327) q[3];
sx q[3];
rz(-2.3655112) q[3];
sx q[3];
rz(0.14296338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(3.013124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0634154) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(2.892911) q[2];
sx q[2];
rz(-0.20771019) q[2];
sx q[2];
rz(-0.50616979) q[2];
rz(-0.29240378) q[3];
sx q[3];
rz(-0.3429827) q[3];
sx q[3];
rz(-0.61957785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
