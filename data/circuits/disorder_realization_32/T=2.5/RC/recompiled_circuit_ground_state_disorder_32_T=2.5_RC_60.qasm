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
rz(2.7453121) q[1];
sx q[1];
rz(-3.0488465) q[1];
sx q[1];
rz(2.588811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39550135) q[0];
sx q[0];
rz(-1.1679018) q[0];
sx q[0];
rz(-2.7159116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80443126) q[2];
sx q[2];
rz(-2.0896157) q[2];
sx q[2];
rz(-0.08229736) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1172792) q[1];
sx q[1];
rz(-1.2543884) q[1];
sx q[1];
rz(-2.1561619) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6951233) q[3];
sx q[3];
rz(-0.95799082) q[3];
sx q[3];
rz(-1.0366576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7736241) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(-3.0527414) q[2];
rz(2.7446274) q[3];
sx q[3];
rz(-2.1019955) q[3];
sx q[3];
rz(1.6840434) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001584) q[0];
sx q[0];
rz(-2.6341697) q[0];
sx q[0];
rz(-2.2870824) q[0];
rz(0.62141934) q[1];
sx q[1];
rz(-2.8050551) q[1];
sx q[1];
rz(-2.0889166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54224211) q[0];
sx q[0];
rz(-0.17048888) q[0];
sx q[0];
rz(-2.4986099) q[0];
x q[1];
rz(2.0172878) q[2];
sx q[2];
rz(-1.948775) q[2];
sx q[2];
rz(-1.93731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0704144) q[1];
sx q[1];
rz(-2.2885808) q[1];
sx q[1];
rz(0.26805265) q[1];
x q[2];
rz(-3.1399473) q[3];
sx q[3];
rz(-0.96685322) q[3];
sx q[3];
rz(1.2206447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66775995) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(2.1916126) q[2];
rz(1.0085603) q[3];
sx q[3];
rz(-2.5086094) q[3];
sx q[3];
rz(0.27957255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(1.3470294) q[0];
rz(2.0836209) q[1];
sx q[1];
rz(-0.24669138) q[1];
sx q[1];
rz(-3.0314441) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[1];
x q[1];
rz(1.0588689) q[2];
sx q[2];
rz(-2.9721568) q[2];
sx q[2];
rz(-1.2537341) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55787841) q[1];
sx q[1];
rz(-1.4643351) q[1];
sx q[1];
rz(0.6006247) q[1];
rz(-pi) q[2];
rz(-2.3055656) q[3];
sx q[3];
rz(-2.5213833) q[3];
sx q[3];
rz(2.7171752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9614253) q[2];
sx q[2];
rz(-1.5632997) q[2];
sx q[2];
rz(3.0999198) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(-0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623077) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(0.98096171) q[0];
rz(0.43308577) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(-1.513419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0312322) q[0];
sx q[0];
rz(-1.1735532) q[0];
sx q[0];
rz(-1.7658986) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89940091) q[2];
sx q[2];
rz(-1.581462) q[2];
sx q[2];
rz(-0.15029066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95483741) q[1];
sx q[1];
rz(-1.4680956) q[1];
sx q[1];
rz(-1.0726733) q[1];
rz(0.93928503) q[3];
sx q[3];
rz(-1.3399249) q[3];
sx q[3];
rz(2.758647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5341586) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(2.329211) q[2];
rz(-1.8937998) q[3];
sx q[3];
rz(-1.2500074) q[3];
sx q[3];
rz(1.8947424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8060551) q[0];
sx q[0];
rz(-1.0220818) q[0];
sx q[0];
rz(-1.442765) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.6280326) q[1];
sx q[1];
rz(-2.3853669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6578015) q[0];
sx q[0];
rz(-1.7258125) q[0];
sx q[0];
rz(0.94021262) q[0];
rz(-2.4632023) q[2];
sx q[2];
rz(-2.3783461) q[2];
sx q[2];
rz(2.2445172) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56179433) q[1];
sx q[1];
rz(-1.4358014) q[1];
sx q[1];
rz(-1.6588103) q[1];
rz(-1.4706572) q[3];
sx q[3];
rz(-0.34557811) q[3];
sx q[3];
rz(-2.341193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68236399) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(-2.1759822) q[2];
rz(-2.5718578) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7392893) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(1.7560316) q[0];
rz(0.7631453) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(3.0398583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0245314) q[0];
sx q[0];
rz(-2.5145338) q[0];
sx q[0];
rz(1.4133021) q[0];
rz(3.0252803) q[2];
sx q[2];
rz(-2.7564133) q[2];
sx q[2];
rz(1.0508089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49470529) q[1];
sx q[1];
rz(-0.27783074) q[1];
sx q[1];
rz(1.3459413) q[1];
x q[2];
rz(2.2564628) q[3];
sx q[3];
rz(-0.99381522) q[3];
sx q[3];
rz(-2.5654405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0330641) q[2];
sx q[2];
rz(-1.4427002) q[2];
sx q[2];
rz(0.70890439) q[2];
rz(2.3809643) q[3];
sx q[3];
rz(-1.1516738) q[3];
sx q[3];
rz(-0.93543783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075542299) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(-0.72108889) q[0];
rz(2.1636294) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.579938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52873351) q[0];
sx q[0];
rz(-1.4943143) q[0];
sx q[0];
rz(1.4974817) q[0];
rz(-pi) q[1];
rz(1.220854) q[2];
sx q[2];
rz(-2.5169543) q[2];
sx q[2];
rz(0.24764316) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5389293) q[1];
sx q[1];
rz(-0.70405761) q[1];
sx q[1];
rz(-2.5119147) q[1];
rz(2.2411602) q[3];
sx q[3];
rz(-1.2868902) q[3];
sx q[3];
rz(0.30208029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-2.6140116) q[2];
sx q[2];
rz(1.52012) q[2];
rz(0.4367477) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(-2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590416) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(-2.3178597) q[0];
rz(1.1388904) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(-1.9653856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9439119) q[0];
sx q[0];
rz(-2.072654) q[0];
sx q[0];
rz(-0.051086257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9763821) q[2];
sx q[2];
rz(-2.1422804) q[2];
sx q[2];
rz(-1.6125082) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8835134) q[1];
sx q[1];
rz(-2.9032859) q[1];
sx q[1];
rz(3.1001453) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0823364) q[3];
sx q[3];
rz(-2.5496581) q[3];
sx q[3];
rz(-1.2436641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1463683) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(-2.5541019) q[2];
rz(-2.7557709) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-0.90252701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189482) q[0];
sx q[0];
rz(-2.1653439) q[0];
sx q[0];
rz(0.60978419) q[0];
rz(-1.7976286) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(-2.9885898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5333819) q[0];
sx q[0];
rz(-2.675288) q[0];
sx q[0];
rz(3.0979041) q[0];
rz(-pi) q[1];
rz(-3.1315003) q[2];
sx q[2];
rz(-2.3463425) q[2];
sx q[2];
rz(-0.45999664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3045522) q[1];
sx q[1];
rz(-2.4992832) q[1];
sx q[1];
rz(-1.539597) q[1];
x q[2];
rz(-2.6681247) q[3];
sx q[3];
rz(-0.93515474) q[3];
sx q[3];
rz(-1.1421695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0387705) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(0.62225303) q[2];
rz(-1.9410939) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95933652) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(-1.1025053) q[1];
sx q[1];
rz(-0.016977221) q[1];
sx q[1];
rz(2.4818518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3400354) q[0];
sx q[0];
rz(-1.3990972) q[0];
sx q[0];
rz(0.6562018) q[0];
x q[1];
rz(-0.12730403) q[2];
sx q[2];
rz(-1.1603171) q[2];
sx q[2];
rz(0.32038996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3952418) q[1];
sx q[1];
rz(-1.8971656) q[1];
sx q[1];
rz(-2.8119948) q[1];
rz(-pi) q[2];
rz(-0.052614502) q[3];
sx q[3];
rz(-0.89375118) q[3];
sx q[3];
rz(0.074124215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79992646) q[2];
sx q[2];
rz(-3.0604) q[2];
sx q[2];
rz(2.9426835) q[2];
rz(-2.343446) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(1.2822436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9195008) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(-2.3566698) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(-1.6941407) q[2];
sx q[2];
rz(-1.5236749) q[2];
sx q[2];
rz(2.9070367) q[2];
rz(-1.296923) q[3];
sx q[3];
rz(-0.92007888) q[3];
sx q[3];
rz(-2.32213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
