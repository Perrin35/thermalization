OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(-2.0071964) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463319) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(-2.6803826) q[0];
rz(1.2167591) q[2];
sx q[2];
rz(-2.3166222) q[2];
sx q[2];
rz(1.1302469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2156148) q[1];
sx q[1];
rz(-1.8324592) q[1];
sx q[1];
rz(1.2222626) q[1];
rz(-pi) q[2];
rz(-1.8144242) q[3];
sx q[3];
rz(-1.4654136) q[3];
sx q[3];
rz(1.6107314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2873043) q[0];
sx q[0];
rz(-2.0137631) q[0];
sx q[0];
rz(-1.7914377) q[0];
rz(-0.77483564) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(1.6006084) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5920168) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(-1.7906584) q[1];
rz(0.0094718178) q[3];
sx q[3];
rz(-2.1433899) q[3];
sx q[3];
rz(-0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817292) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-2.1121315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345966) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(0.2683123) q[0];
x q[1];
rz(2.4918409) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(-2.8284555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0440812) q[1];
sx q[1];
rz(-1.4315726) q[1];
sx q[1];
rz(2.984725) q[1];
rz(-pi) q[2];
rz(-3.0089278) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(-2.0987299) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.5159336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(2.9177833) q[0];
x q[1];
rz(3.0558673) q[2];
sx q[2];
rz(-1.013373) q[2];
sx q[2];
rz(-2.2103708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3486466) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(1.2299041) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3514148) q[3];
sx q[3];
rz(-2.1483768) q[3];
sx q[3];
rz(-0.95336174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1631854) q[0];
sx q[0];
rz(-2.1911231) q[0];
sx q[0];
rz(-2.083367) q[0];
x q[1];
rz(-2.4159555) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(-2.5051136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85256165) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(-2.186071) q[1];
rz(-pi) q[2];
rz(-0.026167913) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.3775795) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-2.6766052) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842448) q[0];
sx q[0];
rz(-1.2372969) q[0];
sx q[0];
rz(-1.1112945) q[0];
rz(-pi) q[1];
rz(2.8000185) q[2];
sx q[2];
rz(-1.5014868) q[2];
sx q[2];
rz(-2.9937033) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8919233) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(-3.0628672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(0.15596095) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0322745) q[0];
sx q[0];
rz(-2.8335857) q[0];
sx q[0];
rz(0.59662915) q[0];
rz(-pi) q[1];
rz(-2.6290226) q[2];
sx q[2];
rz(-1.9480431) q[2];
sx q[2];
rz(0.23114983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29282001) q[1];
sx q[1];
rz(-2.9584868) q[1];
sx q[1];
rz(1.8715026) q[1];
rz(-1.0714674) q[3];
sx q[3];
rz(-2.0519749) q[3];
sx q[3];
rz(-2.1035946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(-1.0726661) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-2.8930194) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8817431) q[0];
sx q[0];
rz(-1.3215085) q[0];
sx q[0];
rz(2.6134406) q[0];
rz(-pi) q[1];
rz(1.0609264) q[2];
sx q[2];
rz(-1.9976127) q[2];
sx q[2];
rz(-0.61444297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.12411815) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(-0.89213051) q[1];
rz(1.400984) q[3];
sx q[3];
rz(-1.6465934) q[3];
sx q[3];
rz(2.4236987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(2.8588262) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.3185906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193664) q[0];
sx q[0];
rz(-2.0113693) q[0];
sx q[0];
rz(0.088935436) q[0];
rz(-0.96005) q[2];
sx q[2];
rz(-2.1337482) q[2];
sx q[2];
rz(0.27573953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9642155) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(-1.7978976) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9570458) q[3];
sx q[3];
rz(-0.85856122) q[3];
sx q[3];
rz(-1.4854747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(-1.261196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1051837) q[0];
sx q[0];
rz(-1.9793708) q[0];
sx q[0];
rz(-1.6856134) q[0];
rz(-pi) q[1];
rz(0.59826675) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(-2.604904) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.48325) q[1];
sx q[1];
rz(-1.5066506) q[1];
sx q[1];
rz(-2.642753) q[1];
rz(2.946978) q[3];
sx q[3];
rz(-0.6932887) q[3];
sx q[3];
rz(0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-0.99909487) q[2];
sx q[2];
rz(-1.5487164) q[2];
sx q[2];
rz(-1.6185417) q[2];
rz(1.431987) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];