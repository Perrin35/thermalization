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
rz(1.3067955) q[0];
sx q[0];
rz(-0.84796325) q[0];
sx q[0];
rz(-0.037394878) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(1.4451292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045528) q[0];
sx q[0];
rz(-1.7524961) q[0];
sx q[0];
rz(-0.43614014) q[0];
x q[1];
rz(3.0710241) q[2];
sx q[2];
rz(-1.6733992) q[2];
sx q[2];
rz(3.0377394) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8822669) q[1];
sx q[1];
rz(-1.0333459) q[1];
sx q[1];
rz(0.56973704) q[1];
rz(-3.0274452) q[3];
sx q[3];
rz(-2.1329013) q[3];
sx q[3];
rz(-1.8207681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41363132) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(1.5067345) q[2];
rz(-0.19351752) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(-2.1957652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0993318) q[0];
sx q[0];
rz(-1.3657382) q[0];
sx q[0];
rz(-0.17901626) q[0];
rz(1.6049113) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(1.590033) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5478599) q[0];
sx q[0];
rz(-1.6681801) q[0];
sx q[0];
rz(-1.2800526) q[0];
rz(-0.29921542) q[2];
sx q[2];
rz(-0.73198527) q[2];
sx q[2];
rz(0.74886403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.034160567) q[1];
sx q[1];
rz(-1.9692268) q[1];
sx q[1];
rz(0.28657367) q[1];
x q[2];
rz(-0.47081866) q[3];
sx q[3];
rz(-1.075346) q[3];
sx q[3];
rz(1.9268056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4701074) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(3.1191077) q[2];
rz(-2.2454028) q[3];
sx q[3];
rz(-2.360207) q[3];
sx q[3];
rz(-1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.80980587) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(1.608954) q[0];
rz(-1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(-0.87361139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7067817) q[0];
sx q[0];
rz(-1.0651089) q[0];
sx q[0];
rz(-2.3825453) q[0];
rz(-0.64675764) q[2];
sx q[2];
rz(-2.0271747) q[2];
sx q[2];
rz(2.6025085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70097476) q[1];
sx q[1];
rz(-1.7243885) q[1];
sx q[1];
rz(0.48753341) q[1];
x q[2];
rz(-2.8779712) q[3];
sx q[3];
rz(-2.859451) q[3];
sx q[3];
rz(-0.96708114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5736299) q[2];
sx q[2];
rz(-0.45845389) q[2];
sx q[2];
rz(2.6652794) q[2];
rz(1.025398) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(-2.8477113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697701) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(2.5452132) q[0];
rz(2.3884933) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(-2.7900043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.387334) q[0];
sx q[0];
rz(-1.046858) q[0];
sx q[0];
rz(1.7225527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2542721) q[2];
sx q[2];
rz(-1.9565689) q[2];
sx q[2];
rz(-1.9172457) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6476485) q[1];
sx q[1];
rz(-0.46174368) q[1];
sx q[1];
rz(-1.0346725) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91686317) q[3];
sx q[3];
rz(-1.0427949) q[3];
sx q[3];
rz(2.7478086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4190462) q[2];
sx q[2];
rz(-1.6435577) q[2];
sx q[2];
rz(-0.92918116) q[2];
rz(1.2292713) q[3];
sx q[3];
rz(-3.1049187) q[3];
sx q[3];
rz(1.8721972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1790328) q[0];
sx q[0];
rz(-0.61229175) q[0];
sx q[0];
rz(-2.6137733) q[0];
rz(2.5714696) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(-2.747587) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0150945) q[0];
sx q[0];
rz(-1.0929035) q[0];
sx q[0];
rz(-1.7897578) q[0];
rz(-1.9377685) q[2];
sx q[2];
rz(-1.0093371) q[2];
sx q[2];
rz(-2.9740564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24007639) q[1];
sx q[1];
rz(-1.6887159) q[1];
sx q[1];
rz(1.7941956) q[1];
rz(-2.3218003) q[3];
sx q[3];
rz(-0.91320437) q[3];
sx q[3];
rz(-0.70943225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7492619) q[2];
sx q[2];
rz(-1.2886084) q[2];
sx q[2];
rz(-2.0070845) q[2];
rz(0.025731651) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(0.64190763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57539097) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(-0.45912418) q[0];
rz(-2.3205914) q[1];
sx q[1];
rz(-0.88290015) q[1];
sx q[1];
rz(0.51148907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59508649) q[0];
sx q[0];
rz(-1.0985785) q[0];
sx q[0];
rz(-0.32924633) q[0];
rz(-0.29634997) q[2];
sx q[2];
rz(-0.68542333) q[2];
sx q[2];
rz(2.1330331) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59907179) q[1];
sx q[1];
rz(-0.59172809) q[1];
sx q[1];
rz(-0.27246957) q[1];
x q[2];
rz(-1.7396163) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5795827) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(-0.75135922) q[2];
rz(0.56685081) q[3];
sx q[3];
rz(-1.6040498) q[3];
sx q[3];
rz(-1.3273299) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.421627) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(-0.0075465329) q[0];
rz(-2.201572) q[1];
sx q[1];
rz(-2.4766141) q[1];
sx q[1];
rz(3.0090581) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7791393) q[0];
sx q[0];
rz(-0.7374239) q[0];
sx q[0];
rz(2.2064184) q[0];
rz(-pi) q[1];
rz(0.36302249) q[2];
sx q[2];
rz(-1.6953812) q[2];
sx q[2];
rz(1.4113791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79515565) q[1];
sx q[1];
rz(-0.43722414) q[1];
sx q[1];
rz(-0.45262419) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9136732) q[3];
sx q[3];
rz(-1.7352805) q[3];
sx q[3];
rz(0.13492385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14564766) q[2];
sx q[2];
rz(-1.4952679) q[2];
sx q[2];
rz(0.10124595) q[2];
rz(2.4920987) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(-1.2671965) q[0];
rz(-1.2665117) q[1];
sx q[1];
rz(-0.15788618) q[1];
sx q[1];
rz(0.74323851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38108724) q[0];
sx q[0];
rz(-0.60956565) q[0];
sx q[0];
rz(0.35778348) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8536123) q[2];
sx q[2];
rz(-2.3415945) q[2];
sx q[2];
rz(1.204042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.083588138) q[1];
sx q[1];
rz(-2.8561855) q[1];
sx q[1];
rz(2.1607375) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.322213) q[3];
sx q[3];
rz(-1.0257105) q[3];
sx q[3];
rz(-1.2631288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1242975) q[2];
sx q[2];
rz(-0.20169078) q[2];
sx q[2];
rz(1.7151493) q[2];
rz(-1.0157478) q[3];
sx q[3];
rz(-1.7044715) q[3];
sx q[3];
rz(-0.82971853) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5947241) q[0];
sx q[0];
rz(-0.65383738) q[0];
sx q[0];
rz(0.015901707) q[0];
rz(1.5025899) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(-2.1471088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8581945) q[0];
sx q[0];
rz(-2.1680729) q[0];
sx q[0];
rz(-2.6676763) q[0];
x q[1];
rz(2.7632668) q[2];
sx q[2];
rz(-1.6294453) q[2];
sx q[2];
rz(0.24967061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82518903) q[1];
sx q[1];
rz(-1.087552) q[1];
sx q[1];
rz(-2.2903328) q[1];
rz(-pi) q[2];
x q[2];
rz(1.307147) q[3];
sx q[3];
rz(-0.52894512) q[3];
sx q[3];
rz(-0.45252132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.7113643) q[2];
sx q[2];
rz(-2.963781) q[2];
rz(-1.4583679) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(1.873707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.6630702) q[0];
sx q[0];
rz(-1.8926184) q[0];
sx q[0];
rz(1.2836237) q[0];
rz(-0.78370699) q[1];
sx q[1];
rz(-2.3678534) q[1];
sx q[1];
rz(1.3073889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2268902) q[0];
sx q[0];
rz(-0.44825867) q[0];
sx q[0];
rz(2.9284555) q[0];
rz(-1.9587014) q[2];
sx q[2];
rz(-1.5989306) q[2];
sx q[2];
rz(-1.494734) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5518903) q[1];
sx q[1];
rz(-1.3968727) q[1];
sx q[1];
rz(0.059373418) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0211759) q[3];
sx q[3];
rz(-1.2116287) q[3];
sx q[3];
rz(-3.1001665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6002097) q[2];
sx q[2];
rz(-1.7317438) q[2];
sx q[2];
rz(-2.6516338) q[2];
rz(-1.1146891) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212696) q[0];
sx q[0];
rz(-1.2280432) q[0];
sx q[0];
rz(-1.8228774) q[0];
rz(1.2976788) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(2.8823356) q[2];
sx q[2];
rz(-1.3217864) q[2];
sx q[2];
rz(-1.9733236) q[2];
rz(-2.277247) q[3];
sx q[3];
rz(-1.723524) q[3];
sx q[3];
rz(1.9836457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
