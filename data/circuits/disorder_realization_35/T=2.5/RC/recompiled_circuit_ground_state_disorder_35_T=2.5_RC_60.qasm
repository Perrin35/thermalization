OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(0.068109186) q[0];
sx q[0];
rz(13.343233) q[0];
rz(1.8344954) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(-1.3219272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711219) q[0];
sx q[0];
rz(-1.4323455) q[0];
sx q[0];
rz(-2.2884503) q[0];
x q[1];
rz(-2.6451557) q[2];
sx q[2];
rz(-1.8775216) q[2];
sx q[2];
rz(-2.6076743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3275431) q[1];
sx q[1];
rz(-0.40990489) q[1];
sx q[1];
rz(1.4973212) q[1];
rz(-pi) q[2];
rz(2.1573651) q[3];
sx q[3];
rz(-0.34474686) q[3];
sx q[3];
rz(2.1380827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(2.784101) q[2];
rz(-0.56719559) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(-0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95328632) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(1.8081007) q[0];
rz(2.7045344) q[1];
sx q[1];
rz(-2.5407365) q[1];
sx q[1];
rz(2.3016047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8020426) q[0];
sx q[0];
rz(-1.7239445) q[0];
sx q[0];
rz(1.0569554) q[0];
rz(-pi) q[1];
rz(-2.8905021) q[2];
sx q[2];
rz(-0.60949627) q[2];
sx q[2];
rz(-1.239038) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4639791) q[1];
sx q[1];
rz(-2.2871454) q[1];
sx q[1];
rz(-2.3542627) q[1];
rz(-pi) q[2];
rz(1.0209818) q[3];
sx q[3];
rz(-1.5375759) q[3];
sx q[3];
rz(-1.1266421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4240894) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(1.7355512) q[2];
rz(-2.9794335) q[3];
sx q[3];
rz(-1.5806961) q[3];
sx q[3];
rz(-2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(0.42174569) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(0.27580321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037791) q[0];
sx q[0];
rz(-0.40291726) q[0];
sx q[0];
rz(-0.62312868) q[0];
rz(-pi) q[1];
rz(-0.47068542) q[2];
sx q[2];
rz(-1.7690725) q[2];
sx q[2];
rz(-0.40775611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.083407) q[1];
sx q[1];
rz(-0.32489932) q[1];
sx q[1];
rz(2.9778381) q[1];
rz(-pi) q[2];
rz(-1.6869322) q[3];
sx q[3];
rz(-0.7051055) q[3];
sx q[3];
rz(0.0034611387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3383823) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(2.4513643) q[2];
rz(2.7817182) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(0.38351044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91367078) q[0];
sx q[0];
rz(-1.0424732) q[0];
sx q[0];
rz(-0.87164718) q[0];
rz(-0.40920416) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86996468) q[0];
sx q[0];
rz(-1.6883435) q[0];
sx q[0];
rz(-2.6948476) q[0];
rz(-1.4710452) q[2];
sx q[2];
rz(-1.776223) q[2];
sx q[2];
rz(-0.86459407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4620381) q[1];
sx q[1];
rz(-1.7741388) q[1];
sx q[1];
rz(2.773079) q[1];
rz(-pi) q[2];
rz(1.7658224) q[3];
sx q[3];
rz(-0.93149501) q[3];
sx q[3];
rz(-2.2488058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(1.4247591) q[3];
sx q[3];
rz(-0.32130876) q[3];
sx q[3];
rz(-0.43535522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7216126) q[0];
sx q[0];
rz(-1.551832) q[0];
sx q[0];
rz(-0.87442526) q[0];
rz(-0.32886109) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(1.0985451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4916353) q[0];
sx q[0];
rz(-1.945561) q[0];
sx q[0];
rz(-1.6894421) q[0];
x q[1];
rz(-0.83580612) q[2];
sx q[2];
rz(-1.5975738) q[2];
sx q[2];
rz(1.4876897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2862773) q[1];
sx q[1];
rz(-0.5868656) q[1];
sx q[1];
rz(-2.1333958) q[1];
x q[2];
rz(0.27733488) q[3];
sx q[3];
rz(-2.1430052) q[3];
sx q[3];
rz(-0.25857833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99981368) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-0.67406526) q[2];
rz(-1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(0.28356799) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0312626) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(-1.7562235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0273756) q[0];
sx q[0];
rz(-1.3189335) q[0];
sx q[0];
rz(-0.92433874) q[0];
rz(-pi) q[1];
rz(0.49853168) q[2];
sx q[2];
rz(-1.4953002) q[2];
sx q[2];
rz(2.1488862) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1277571) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(0.54658039) q[1];
x q[2];
rz(2.9870944) q[3];
sx q[3];
rz(-1.5348001) q[3];
sx q[3];
rz(-2.846623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(0.37609491) q[2];
rz(-2.0070576) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019808708) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(0.10251775) q[0];
rz(-0.58492297) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.9580511) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2823806) q[0];
sx q[0];
rz(-1.6995653) q[0];
sx q[0];
rz(2.9984401) q[0];
x q[1];
rz(-1.2361707) q[2];
sx q[2];
rz(-0.93572223) q[2];
sx q[2];
rz(-2.2872567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2940604) q[1];
sx q[1];
rz(-0.69662635) q[1];
sx q[1];
rz(1.9737592) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46681301) q[3];
sx q[3];
rz(-1.3315505) q[3];
sx q[3];
rz(-1.6166663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3465603) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(1.0605158) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.878433) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(2.8286381) q[0];
rz(-2.1921659) q[1];
sx q[1];
rz(-1.7890309) q[1];
sx q[1];
rz(-1.4535646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4623186) q[0];
sx q[0];
rz(-1.3841713) q[0];
sx q[0];
rz(-1.4786722) q[0];
x q[1];
rz(-0.61954478) q[2];
sx q[2];
rz(-1.6352374) q[2];
sx q[2];
rz(-3.0321414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6387299) q[1];
sx q[1];
rz(-2.4045378) q[1];
sx q[1];
rz(2.8046737) q[1];
x q[2];
rz(-1.9096776) q[3];
sx q[3];
rz(-1.1383071) q[3];
sx q[3];
rz(0.56604715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7877385) q[2];
sx q[2];
rz(-0.9757897) q[2];
sx q[2];
rz(-1.806951) q[2];
rz(1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036309328) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(-0.69806725) q[0];
rz(-0.87567323) q[1];
sx q[1];
rz(-1.061729) q[1];
sx q[1];
rz(-2.7811513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35931906) q[0];
sx q[0];
rz(-0.8018612) q[0];
sx q[0];
rz(-0.62863056) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7588604) q[2];
sx q[2];
rz(-1.6305411) q[2];
sx q[2];
rz(2.003423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0086509) q[1];
sx q[1];
rz(-1.1906149) q[1];
sx q[1];
rz(-0.93312414) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8713851) q[3];
sx q[3];
rz(-2.5761309) q[3];
sx q[3];
rz(0.61494857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91161072) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(-2.3243375) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(2.9221007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(3.1095374) q[0];
rz(-1.8219148) q[1];
sx q[1];
rz(-1.7664884) q[1];
sx q[1];
rz(0.65418902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830156) q[0];
sx q[0];
rz(-2.9111283) q[0];
sx q[0];
rz(2.5925267) q[0];
x q[1];
rz(2.3532969) q[2];
sx q[2];
rz(-1.239778) q[2];
sx q[2];
rz(2.2904879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75958222) q[1];
sx q[1];
rz(-1.1076265) q[1];
sx q[1];
rz(0.5447556) q[1];
rz(-2.8579162) q[3];
sx q[3];
rz(-0.95796889) q[3];
sx q[3];
rz(3.0949288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0111982) q[2];
sx q[2];
rz(-2.0989213) q[2];
sx q[2];
rz(2.1350258) q[2];
rz(0.69342962) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(1.8600195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(2.2609932) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(-2.6237189) q[2];
sx q[2];
rz(-0.93954115) q[2];
sx q[2];
rz(2.1672135) q[2];
rz(2.4453658) q[3];
sx q[3];
rz(-1.8693557) q[3];
sx q[3];
rz(-2.1341677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
