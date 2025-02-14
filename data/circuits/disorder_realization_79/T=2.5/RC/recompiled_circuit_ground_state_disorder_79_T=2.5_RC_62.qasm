OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5918936) q[0];
sx q[0];
rz(-2.0285719) q[0];
sx q[0];
rz(2.9099303) q[0];
rz(-2.7910233) q[1];
sx q[1];
rz(-2.9313593) q[1];
sx q[1];
rz(-2.3793013) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0526844) q[0];
sx q[0];
rz(-1.1529973) q[0];
sx q[0];
rz(-2.4521928) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76114391) q[2];
sx q[2];
rz(-2.5147438) q[2];
sx q[2];
rz(-0.055920211) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.9673358) q[1];
sx q[1];
rz(-2.9474181) q[1];
sx q[1];
rz(2.8513191) q[1];
rz(-pi) q[2];
rz(-2.2985994) q[3];
sx q[3];
rz(-2.5268433) q[3];
sx q[3];
rz(-0.45045567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98213696) q[2];
sx q[2];
rz(-2.5348713) q[2];
sx q[2];
rz(-0.43178001) q[2];
rz(0.86709705) q[3];
sx q[3];
rz(-0.95621395) q[3];
sx q[3];
rz(1.6352656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302439) q[0];
sx q[0];
rz(-0.83006492) q[0];
sx q[0];
rz(-1.9507971) q[0];
rz(-1.2442773) q[1];
sx q[1];
rz(-2.0121274) q[1];
sx q[1];
rz(-0.58070374) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8180926) q[0];
sx q[0];
rz(-0.43002263) q[0];
sx q[0];
rz(-1.9097206) q[0];
x q[1];
rz(-0.9208619) q[2];
sx q[2];
rz(-2.5061786) q[2];
sx q[2];
rz(-0.17104761) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.001284842) q[1];
sx q[1];
rz(-1.2540037) q[1];
sx q[1];
rz(0.66441925) q[1];
x q[2];
rz(-1.034581) q[3];
sx q[3];
rz(-2.381122) q[3];
sx q[3];
rz(1.5602668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0018491) q[2];
sx q[2];
rz(-1.7364343) q[2];
sx q[2];
rz(1.639036) q[2];
rz(2.9338845) q[3];
sx q[3];
rz(-1.3108871) q[3];
sx q[3];
rz(1.0072072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19596066) q[0];
sx q[0];
rz(-2.0646844) q[0];
sx q[0];
rz(2.9929602) q[0];
rz(-1.0736939) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(-0.78280848) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4769094) q[0];
sx q[0];
rz(-0.17897716) q[0];
sx q[0];
rz(-1.8838691) q[0];
rz(-pi) q[1];
rz(-1.5284848) q[2];
sx q[2];
rz(-1.0382922) q[2];
sx q[2];
rz(1.8489727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69082123) q[1];
sx q[1];
rz(-1.3714002) q[1];
sx q[1];
rz(-2.6197207) q[1];
x q[2];
rz(1.820676) q[3];
sx q[3];
rz(-0.80554038) q[3];
sx q[3];
rz(-0.66242927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9624376) q[2];
sx q[2];
rz(-2.16733) q[2];
sx q[2];
rz(-1.1243189) q[2];
rz(3.1149241) q[3];
sx q[3];
rz(-2.4494438) q[3];
sx q[3];
rz(-0.28287014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7836595) q[0];
sx q[0];
rz(-1.1555576) q[0];
sx q[0];
rz(1.577277) q[0];
rz(2.0951001) q[1];
sx q[1];
rz(-2.107403) q[1];
sx q[1];
rz(1.1276721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8148147) q[0];
sx q[0];
rz(-0.59140697) q[0];
sx q[0];
rz(0.018494292) q[0];
x q[1];
rz(-1.3443742) q[2];
sx q[2];
rz(-2.7993996) q[2];
sx q[2];
rz(2.7777517) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31364031) q[1];
sx q[1];
rz(-1.8871739) q[1];
sx q[1];
rz(0.72319855) q[1];
rz(3.042759) q[3];
sx q[3];
rz(-1.4526794) q[3];
sx q[3];
rz(1.57847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6005738) q[2];
sx q[2];
rz(-1.8709196) q[2];
sx q[2];
rz(-0.18826558) q[2];
rz(2.6514734) q[3];
sx q[3];
rz(-1.6907588) q[3];
sx q[3];
rz(-1.8671487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.984943) q[0];
sx q[0];
rz(-1.5644512) q[0];
sx q[0];
rz(-1.4666602) q[0];
rz(2.6690392) q[1];
sx q[1];
rz(-2.2652389) q[1];
sx q[1];
rz(2.1579425) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60118942) q[0];
sx q[0];
rz(-1.9944257) q[0];
sx q[0];
rz(-1.3715959) q[0];
rz(-1.3784658) q[2];
sx q[2];
rz(-1.1149613) q[2];
sx q[2];
rz(-2.6050732) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2347348) q[1];
sx q[1];
rz(-0.51381293) q[1];
sx q[1];
rz(-1.6870901) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8365588) q[3];
sx q[3];
rz(-1.6739707) q[3];
sx q[3];
rz(2.3931695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5094362) q[2];
sx q[2];
rz(-1.7551273) q[2];
sx q[2];
rz(-0.40718386) q[2];
rz(1.8699649) q[3];
sx q[3];
rz(-2.7224702) q[3];
sx q[3];
rz(-0.65740383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42095175) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(-1.2366914) q[0];
rz(-2.0115133) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(2.0225661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73858628) q[0];
sx q[0];
rz(-1.6386014) q[0];
sx q[0];
rz(2.9878776) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5594312) q[2];
sx q[2];
rz(-0.88554875) q[2];
sx q[2];
rz(0.2174046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0862217) q[1];
sx q[1];
rz(-1.9920336) q[1];
sx q[1];
rz(1.4122541) q[1];
rz(-1.9672333) q[3];
sx q[3];
rz(-1.9810976) q[3];
sx q[3];
rz(-1.0956216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5469024) q[2];
sx q[2];
rz(-1.12744) q[2];
sx q[2];
rz(-2.8315869) q[2];
rz(1.4340495) q[3];
sx q[3];
rz(-1.5177582) q[3];
sx q[3];
rz(-1.2448509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2731648) q[0];
sx q[0];
rz(-2.7719066) q[0];
sx q[0];
rz(-2.0694859) q[0];
rz(-1.3605236) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(1.2881813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38561197) q[0];
sx q[0];
rz(-2.1588584) q[0];
sx q[0];
rz(2.5137492) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4425051) q[2];
sx q[2];
rz(-1.0396084) q[2];
sx q[2];
rz(-0.99260222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9065422) q[1];
sx q[1];
rz(-1.8120011) q[1];
sx q[1];
rz(1.3969087) q[1];
rz(-pi) q[2];
rz(2.0615408) q[3];
sx q[3];
rz(-2.2778058) q[3];
sx q[3];
rz(2.3243927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0897022) q[2];
sx q[2];
rz(-1.5081729) q[2];
sx q[2];
rz(-1.4349597) q[2];
rz(2.6658304) q[3];
sx q[3];
rz(-0.69403726) q[3];
sx q[3];
rz(-2.7580875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.39110228) q[0];
sx q[0];
rz(-0.35677156) q[0];
sx q[0];
rz(1.179689) q[0];
rz(0.36243311) q[1];
sx q[1];
rz(-1.06203) q[1];
sx q[1];
rz(1.5865883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060136) q[0];
sx q[0];
rz(-1.7227253) q[0];
sx q[0];
rz(-0.32468995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4612975) q[2];
sx q[2];
rz(-1.1914502) q[2];
sx q[2];
rz(2.2915886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4175222) q[1];
sx q[1];
rz(-2.2322025) q[1];
sx q[1];
rz(2.577642) q[1];
rz(0.15896564) q[3];
sx q[3];
rz(-1.0086802) q[3];
sx q[3];
rz(-1.5288439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72593752) q[2];
sx q[2];
rz(-2.026181) q[2];
sx q[2];
rz(-1.9647145) q[2];
rz(2.7365909) q[3];
sx q[3];
rz(-2.160725) q[3];
sx q[3];
rz(2.7485031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554351) q[0];
sx q[0];
rz(-0.09849184) q[0];
sx q[0];
rz(-1.3354907) q[0];
rz(-2.2291741) q[1];
sx q[1];
rz(-1.4211979) q[1];
sx q[1];
rz(0.040249912) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089913) q[0];
sx q[0];
rz(-0.17290056) q[0];
sx q[0];
rz(1.7937051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77004899) q[2];
sx q[2];
rz(-0.75919861) q[2];
sx q[2];
rz(-0.46261132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98513888) q[1];
sx q[1];
rz(-0.6985526) q[1];
sx q[1];
rz(0.88628873) q[1];
x q[2];
rz(-2.9222708) q[3];
sx q[3];
rz(-2.3613033) q[3];
sx q[3];
rz(1.8530451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5361629) q[2];
sx q[2];
rz(-1.9998113) q[2];
sx q[2];
rz(2.6119168) q[2];
rz(-1.4179432) q[3];
sx q[3];
rz(-0.46897408) q[3];
sx q[3];
rz(2.6403707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.0025075992) q[0];
sx q[0];
rz(-1.3676099) q[0];
sx q[0];
rz(0.031877192) q[0];
rz(1.9335951) q[1];
sx q[1];
rz(-0.54787689) q[1];
sx q[1];
rz(2.8338103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8610749) q[0];
sx q[0];
rz(-1.8372179) q[0];
sx q[0];
rz(-2.0841588) q[0];
x q[1];
rz(1.304302) q[2];
sx q[2];
rz(-1.6937758) q[2];
sx q[2];
rz(-0.83713712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3584328) q[1];
sx q[1];
rz(-0.98020411) q[1];
sx q[1];
rz(0.62548754) q[1];
rz(-pi) q[2];
rz(-0.14282222) q[3];
sx q[3];
rz(-1.4429314) q[3];
sx q[3];
rz(-0.82500237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.53177437) q[2];
sx q[2];
rz(-2.4419624) q[2];
sx q[2];
rz(-0.043665234) q[2];
rz(1.6627436) q[3];
sx q[3];
rz(-2.6378529) q[3];
sx q[3];
rz(1.3879205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8762348) q[0];
sx q[0];
rz(-1.4690514) q[0];
sx q[0];
rz(1.5734191) q[0];
rz(0.028989446) q[1];
sx q[1];
rz(-2.034076) q[1];
sx q[1];
rz(-2.1709002) q[1];
rz(1.7912279) q[2];
sx q[2];
rz(-1.6904808) q[2];
sx q[2];
rz(-1.3338989) q[2];
rz(-2.0652931) q[3];
sx q[3];
rz(-0.87135472) q[3];
sx q[3];
rz(3.0490124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
