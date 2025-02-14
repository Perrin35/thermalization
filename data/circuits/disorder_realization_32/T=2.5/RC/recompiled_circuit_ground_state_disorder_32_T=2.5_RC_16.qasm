OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(5.8869047) q[1];
sx q[1];
rz(6.1904391) q[1];
sx q[1];
rz(9.9775597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39550135) q[0];
sx q[0];
rz(-1.9736909) q[0];
sx q[0];
rz(-2.7159116) q[0];
rz(2.3371614) q[2];
sx q[2];
rz(-1.0519769) q[2];
sx q[2];
rz(-0.08229736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0243135) q[1];
sx q[1];
rz(-1.2543884) q[1];
sx q[1];
rz(0.98543075) q[1];
rz(1.0200653) q[3];
sx q[3];
rz(-2.4006776) q[3];
sx q[3];
rz(-1.4108085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7736241) q[2];
sx q[2];
rz(-1.6398733) q[2];
sx q[2];
rz(3.0527414) q[2];
rz(0.39696524) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001584) q[0];
sx q[0];
rz(-2.6341697) q[0];
sx q[0];
rz(-0.85451025) q[0];
rz(-0.62141934) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(1.052676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7490279) q[0];
sx q[0];
rz(-1.6727007) q[0];
sx q[0];
rz(-3.0046731) q[0];
rz(-pi) q[1];
rz(-0.82727716) q[2];
sx q[2];
rz(-2.5649539) q[2];
sx q[2];
rz(2.1183543) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32089511) q[1];
sx q[1];
rz(-1.7716494) q[1];
sx q[1];
rz(-0.83495514) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0016453513) q[3];
sx q[3];
rz(-0.96685322) q[3];
sx q[3];
rz(1.920948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66775995) q[2];
sx q[2];
rz(-2.2809873) q[2];
sx q[2];
rz(2.1916126) q[2];
rz(2.1330323) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(2.8620201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968593) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(1.7945633) q[0];
rz(1.0579717) q[1];
sx q[1];
rz(-0.24669138) q[1];
sx q[1];
rz(3.0314441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059123813) q[0];
sx q[0];
rz(-2.0284855) q[0];
sx q[0];
rz(-0.35510606) q[0];
rz(3.0579849) q[2];
sx q[2];
rz(-1.4232529) q[2];
sx q[2];
rz(-1.7718441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0855838) q[1];
sx q[1];
rz(-0.97404503) q[1];
sx q[1];
rz(-1.6996154) q[1];
rz(0.83602704) q[3];
sx q[3];
rz(-2.5213833) q[3];
sx q[3];
rz(2.7171752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9614253) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(0.041672826) q[2];
rz(0.20798072) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(-0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623077) q[0];
sx q[0];
rz(-1.3212181) q[0];
sx q[0];
rz(2.1606309) q[0];
rz(0.43308577) q[1];
sx q[1];
rz(-1.4739477) q[1];
sx q[1];
rz(-1.6281737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53674066) q[0];
sx q[0];
rz(-1.3910595) q[0];
sx q[0];
rz(0.40412235) q[0];
x q[1];
rz(3.1279706) q[2];
sx q[2];
rz(-2.2421466) q[2];
sx q[2];
rz(-1.7126132) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42964307) q[1];
sx q[1];
rz(-2.6338661) q[1];
sx q[1];
rz(-1.3583378) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1918729) q[3];
sx q[3];
rz(-2.4746568) q[3];
sx q[3];
rz(-2.2569424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5341586) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(-2.329211) q[2];
rz(-1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(-1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3355376) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(-1.6988276) q[0];
rz(-2.6834148) q[1];
sx q[1];
rz(-1.6280326) q[1];
sx q[1];
rz(-2.3853669) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4837912) q[0];
sx q[0];
rz(-1.4157802) q[0];
sx q[0];
rz(0.94021262) q[0];
rz(-pi) q[1];
rz(0.67839038) q[2];
sx q[2];
rz(-2.3783461) q[2];
sx q[2];
rz(-0.89707546) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9994558) q[1];
sx q[1];
rz(-0.161006) q[1];
sx q[1];
rz(-2.5672037) q[1];
x q[2];
rz(-1.2268158) q[3];
sx q[3];
rz(-1.6046673) q[3];
sx q[3];
rz(0.86465166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68236399) q[2];
sx q[2];
rz(-1.7533147) q[2];
sx q[2];
rz(-0.96561042) q[2];
rz(-0.56973488) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(-1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7392893) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(-1.7560316) q[0];
rz(0.7631453) q[1];
sx q[1];
rz(-1.6940073) q[1];
sx q[1];
rz(-3.0398583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0245314) q[0];
sx q[0];
rz(-0.62705886) q[0];
sx q[0];
rz(1.7282906) q[0];
rz(-0.11631233) q[2];
sx q[2];
rz(-2.7564133) q[2];
sx q[2];
rz(1.0508089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6468874) q[1];
sx q[1];
rz(-2.8637619) q[1];
sx q[1];
rz(1.3459413) q[1];
x q[2];
rz(2.3699573) q[3];
sx q[3];
rz(-2.2766678) q[3];
sx q[3];
rz(-2.7350458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0330641) q[2];
sx q[2];
rz(-1.4427002) q[2];
sx q[2];
rz(-0.70890439) q[2];
rz(-2.3809643) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(2.2061548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0660504) q[0];
sx q[0];
rz(-1.1533371) q[0];
sx q[0];
rz(2.4205038) q[0];
rz(-2.1636294) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.5616547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2943731) q[0];
sx q[0];
rz(-0.10589639) q[0];
sx q[0];
rz(-0.76283331) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.220854) q[2];
sx q[2];
rz(-2.5169543) q[2];
sx q[2];
rz(-0.24764316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60266337) q[1];
sx q[1];
rz(-0.70405761) q[1];
sx q[1];
rz(-0.62967794) q[1];
x q[2];
rz(-2.0098814) q[3];
sx q[3];
rz(-0.71936047) q[3];
sx q[3];
rz(-2.2124825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(-1.52012) q[2];
rz(2.704845) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.082551) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(0.82373291) q[0];
rz(1.1388904) q[1];
sx q[1];
rz(-1.3528119) q[1];
sx q[1];
rz(1.9653856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19768077) q[0];
sx q[0];
rz(-1.0689387) q[0];
sx q[0];
rz(-0.051086257) q[0];
rz(-pi) q[1];
rz(0.55032879) q[2];
sx q[2];
rz(-2.4541353) q[2];
sx q[2];
rz(-2.2000809) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35299435) q[1];
sx q[1];
rz(-1.561015) q[1];
sx q[1];
rz(-0.23810975) q[1];
rz(-pi) q[2];
rz(1.6105936) q[3];
sx q[3];
rz(-2.1615513) q[3];
sx q[3];
rz(-1.9692957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99522432) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(2.5541019) q[2];
rz(2.7557709) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189482) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(2.5318085) q[0];
rz(-1.3439641) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(-0.15300289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6082108) q[0];
sx q[0];
rz(-2.675288) q[0];
sx q[0];
rz(-0.043688579) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010092334) q[2];
sx q[2];
rz(-2.3463425) q[2];
sx q[2];
rz(-0.45999664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3435095) q[1];
sx q[1];
rz(-2.2127418) q[1];
sx q[1];
rz(-0.023333418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0172264) q[3];
sx q[3];
rz(-0.77249902) q[3];
sx q[3];
rz(1.8541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1028221) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(-0.62225303) q[2];
rz(-1.2004987) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(-0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(-2.0390873) q[1];
sx q[1];
rz(-0.016977221) q[1];
sx q[1];
rz(0.65974081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099946215) q[0];
sx q[0];
rz(-0.92586854) q[0];
sx q[0];
rz(-1.3553331) q[0];
rz(-pi) q[1];
rz(-1.1573359) q[2];
sx q[2];
rz(-1.454118) q[2];
sx q[2];
rz(-1.199374) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28478947) q[1];
sx q[1];
rz(-1.2592013) q[1];
sx q[1];
rz(1.9143399) q[1];
rz(-pi) q[2];
rz(-2.2485178) q[3];
sx q[3];
rz(-1.6117981) q[3];
sx q[3];
rz(1.6779016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3416662) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(-2.9426835) q[2];
rz(-2.343446) q[3];
sx q[3];
rz(-1.1856368) q[3];
sx q[3];
rz(-1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220919) q[0];
sx q[0];
rz(-1.6035447) q[0];
sx q[0];
rz(-1.0731687) q[0];
rz(0.78492289) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(1.9368108) q[2];
sx q[2];
rz(-3.0095965) q[2];
sx q[2];
rz(-2.1684564) q[2];
rz(-1.8446696) q[3];
sx q[3];
rz(-2.2215138) q[3];
sx q[3];
rz(0.81946269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
