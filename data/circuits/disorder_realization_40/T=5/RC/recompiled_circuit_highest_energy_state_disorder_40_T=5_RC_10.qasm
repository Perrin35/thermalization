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
rz(-1.4200014) q[0];
sx q[0];
rz(5.3237259) q[0];
sx q[0];
rz(11.525679) q[0];
rz(2.8031082) q[1];
sx q[1];
rz(-1.6522633) q[1];
sx q[1];
rz(-1.7117865) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734935) q[0];
sx q[0];
rz(-1.8419319) q[0];
sx q[0];
rz(-1.447729) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0645039) q[2];
sx q[2];
rz(-0.66254751) q[2];
sx q[2];
rz(-2.8665989) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80624455) q[1];
sx q[1];
rz(-1.8523289) q[1];
sx q[1];
rz(2.4068969) q[1];
rz(-0.5654556) q[3];
sx q[3];
rz(-1.3420336) q[3];
sx q[3];
rz(-0.33293399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1374984) q[2];
sx q[2];
rz(-1.4278922) q[2];
sx q[2];
rz(-3.0924228) q[2];
rz(-0.54685012) q[3];
sx q[3];
rz(-0.85421908) q[3];
sx q[3];
rz(1.2949519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8787254) q[0];
sx q[0];
rz(-2.290949) q[0];
sx q[0];
rz(3.1021297) q[0];
rz(-2.4354758) q[1];
sx q[1];
rz(-1.1742274) q[1];
sx q[1];
rz(-1.2605234) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2606756) q[0];
sx q[0];
rz(-1.6291143) q[0];
sx q[0];
rz(-1.1526945) q[0];
rz(-1.0320382) q[2];
sx q[2];
rz(-0.97698394) q[2];
sx q[2];
rz(-2.3597882) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15490048) q[1];
sx q[1];
rz(-3.0873484) q[1];
sx q[1];
rz(2.2510347) q[1];
x q[2];
rz(-0.32716635) q[3];
sx q[3];
rz(-2.9370902) q[3];
sx q[3];
rz(-0.12302264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6588821) q[2];
sx q[2];
rz(-1.005268) q[2];
sx q[2];
rz(-0.91747326) q[2];
rz(0.42258036) q[3];
sx q[3];
rz(-1.8924507) q[3];
sx q[3];
rz(2.1036072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.68040401) q[0];
sx q[0];
rz(-0.55484158) q[0];
sx q[0];
rz(0.94938266) q[0];
rz(-1.8934911) q[1];
sx q[1];
rz(-2.5578942) q[1];
sx q[1];
rz(1.5138352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7010403) q[0];
sx q[0];
rz(-1.3258385) q[0];
sx q[0];
rz(-0.37491684) q[0];
rz(-pi) q[1];
rz(-1.0870088) q[2];
sx q[2];
rz(-1.3788333) q[2];
sx q[2];
rz(-1.0707945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6520034) q[1];
sx q[1];
rz(-1.6449274) q[1];
sx q[1];
rz(-0.28848047) q[1];
rz(-pi) q[2];
rz(-3.0120609) q[3];
sx q[3];
rz(-1.8567652) q[3];
sx q[3];
rz(0.27559552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44024399) q[2];
sx q[2];
rz(-1.1260208) q[2];
sx q[2];
rz(-1.6181642) q[2];
rz(-2.8273888) q[3];
sx q[3];
rz(-1.0256297) q[3];
sx q[3];
rz(1.5020802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1057798) q[0];
sx q[0];
rz(-0.19136763) q[0];
sx q[0];
rz(-2.9845003) q[0];
rz(3.0063903) q[1];
sx q[1];
rz(-2.356485) q[1];
sx q[1];
rz(-0.33321998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858118) q[0];
sx q[0];
rz(-0.69792047) q[0];
sx q[0];
rz(-2.3960953) q[0];
rz(-pi) q[1];
rz(1.7945788) q[2];
sx q[2];
rz(-2.4007041) q[2];
sx q[2];
rz(-0.2303094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1830755) q[1];
sx q[1];
rz(-2.6974899) q[1];
sx q[1];
rz(-2.2217795) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2122194) q[3];
sx q[3];
rz(-0.74508706) q[3];
sx q[3];
rz(2.2555971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80804431) q[2];
sx q[2];
rz(-0.6490038) q[2];
sx q[2];
rz(-1.6928847) q[2];
rz(-0.73271218) q[3];
sx q[3];
rz(-1.2196187) q[3];
sx q[3];
rz(1.609751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919959) q[0];
sx q[0];
rz(-1.4866328) q[0];
sx q[0];
rz(2.9100371) q[0];
rz(-0.72136503) q[1];
sx q[1];
rz(-2.2504579) q[1];
sx q[1];
rz(-0.84842938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68616406) q[0];
sx q[0];
rz(-1.9047184) q[0];
sx q[0];
rz(1.0592106) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6171378) q[2];
sx q[2];
rz(-2.5046621) q[2];
sx q[2];
rz(0.54383155) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49461708) q[1];
sx q[1];
rz(-1.500795) q[1];
sx q[1];
rz(0.38826298) q[1];
x q[2];
rz(-1.6759962) q[3];
sx q[3];
rz(-1.8064515) q[3];
sx q[3];
rz(2.8812129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.007175) q[2];
sx q[2];
rz(-0.60322064) q[2];
sx q[2];
rz(-1.2009386) q[2];
rz(2.4522771) q[3];
sx q[3];
rz(-1.9715693) q[3];
sx q[3];
rz(0.56572604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597964) q[0];
sx q[0];
rz(-1.5317651) q[0];
sx q[0];
rz(1.7933886) q[0];
rz(0.7181522) q[1];
sx q[1];
rz(-0.57128692) q[1];
sx q[1];
rz(-1.8498373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367158) q[0];
sx q[0];
rz(-1.3463536) q[0];
sx q[0];
rz(2.7925452) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0935502) q[2];
sx q[2];
rz(-2.4746918) q[2];
sx q[2];
rz(2.1490522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3499432) q[1];
sx q[1];
rz(-1.7393438) q[1];
sx q[1];
rz(1.4838292) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060543493) q[3];
sx q[3];
rz(-2.745564) q[3];
sx q[3];
rz(2.372449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1527839) q[2];
sx q[2];
rz(-2.1765716) q[2];
sx q[2];
rz(2.7394845) q[2];
rz(1.0058588) q[3];
sx q[3];
rz(-2.918225) q[3];
sx q[3];
rz(1.4373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.724029) q[0];
sx q[0];
rz(-0.14823866) q[0];
sx q[0];
rz(1.8120026) q[0];
rz(-1.8136464) q[1];
sx q[1];
rz(-1.4168394) q[1];
sx q[1];
rz(-0.45164576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76021323) q[0];
sx q[0];
rz(-0.45256786) q[0];
sx q[0];
rz(-0.20968881) q[0];
rz(0.42518227) q[2];
sx q[2];
rz(-1.6984792) q[2];
sx q[2];
rz(-0.38056254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51367769) q[1];
sx q[1];
rz(-1.6213525) q[1];
sx q[1];
rz(1.3804269) q[1];
rz(-pi) q[2];
rz(2.3450646) q[3];
sx q[3];
rz(-0.8332656) q[3];
sx q[3];
rz(-3.0298373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20246501) q[2];
sx q[2];
rz(-1.7697325) q[2];
sx q[2];
rz(-0.28787127) q[2];
rz(1.6599844) q[3];
sx q[3];
rz(-2.2899254) q[3];
sx q[3];
rz(1.6847346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.562029) q[0];
sx q[0];
rz(-0.96051878) q[0];
sx q[0];
rz(-0.55291837) q[0];
rz(1.8909854) q[1];
sx q[1];
rz(-2.7313373) q[1];
sx q[1];
rz(1.7611354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4794248) q[0];
sx q[0];
rz(-1.933455) q[0];
sx q[0];
rz(-0.5525241) q[0];
rz(-pi) q[1];
rz(-3.0199261) q[2];
sx q[2];
rz(-0.98699283) q[2];
sx q[2];
rz(2.9505299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1565754) q[1];
sx q[1];
rz(-2.9010198) q[1];
sx q[1];
rz(2.9678718) q[1];
x q[2];
rz(-0.80069009) q[3];
sx q[3];
rz(-1.4501867) q[3];
sx q[3];
rz(-2.1757954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55014253) q[2];
sx q[2];
rz(-1.3478792) q[2];
sx q[2];
rz(2.4264917) q[2];
rz(-0.044655785) q[3];
sx q[3];
rz(-0.5831334) q[3];
sx q[3];
rz(-2.513212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(1.968068) q[0];
sx q[0];
rz(-2.3723497) q[0];
sx q[0];
rz(-0.64088696) q[0];
rz(1.6454654) q[1];
sx q[1];
rz(-2.75664) q[1];
sx q[1];
rz(2.4900751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2007569) q[0];
sx q[0];
rz(-0.88229942) q[0];
sx q[0];
rz(-2.8066638) q[0];
rz(-pi) q[1];
rz(-1.3369433) q[2];
sx q[2];
rz(-1.0923947) q[2];
sx q[2];
rz(-0.054577915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99371893) q[1];
sx q[1];
rz(-0.9269971) q[1];
sx q[1];
rz(-2.6753475) q[1];
rz(-pi) q[2];
rz(2.9504029) q[3];
sx q[3];
rz(-1.228473) q[3];
sx q[3];
rz(0.99016193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7351825) q[2];
sx q[2];
rz(-2.0969022) q[2];
sx q[2];
rz(-2.6835105) q[2];
rz(3.1297019) q[3];
sx q[3];
rz(-1.0606822) q[3];
sx q[3];
rz(0.021765821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7387725) q[0];
sx q[0];
rz(-0.194855) q[0];
sx q[0];
rz(-3.1126157) q[0];
rz(2.2150691) q[1];
sx q[1];
rz(-1.69311) q[1];
sx q[1];
rz(-0.40625939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0615912) q[0];
sx q[0];
rz(-1.5621981) q[0];
sx q[0];
rz(-0.28010578) q[0];
rz(1.1154956) q[2];
sx q[2];
rz(-0.69820729) q[2];
sx q[2];
rz(-0.9309665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0613289) q[1];
sx q[1];
rz(-0.29326639) q[1];
sx q[1];
rz(1.6786511) q[1];
rz(-pi) q[2];
rz(2.3417937) q[3];
sx q[3];
rz(-2.1635028) q[3];
sx q[3];
rz(2.0519902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77068344) q[2];
sx q[2];
rz(-1.6434881) q[2];
sx q[2];
rz(2.5468199) q[2];
rz(0.70991436) q[3];
sx q[3];
rz(-0.9570595) q[3];
sx q[3];
rz(2.1962568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5668673) q[0];
sx q[0];
rz(-0.95483268) q[0];
sx q[0];
rz(0.86190295) q[0];
rz(2.8361539) q[1];
sx q[1];
rz(-1.3258691) q[1];
sx q[1];
rz(0.23205145) q[1];
rz(-1.1535063) q[2];
sx q[2];
rz(-1.8883033) q[2];
sx q[2];
rz(0.89486833) q[2];
rz(-2.6628982) q[3];
sx q[3];
rz(-2.6810418) q[3];
sx q[3];
rz(-2.2543805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
