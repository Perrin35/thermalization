OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(2.6020738) q[1];
sx q[1];
rz(10.625216) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1737328) q[0];
sx q[0];
rz(-1.1636486) q[0];
sx q[0];
rz(-1.834154) q[0];
rz(1.0317694) q[2];
sx q[2];
rz(-1.4570149) q[2];
sx q[2];
rz(0.10345085) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8564954) q[1];
sx q[1];
rz(-0.76438099) q[1];
sx q[1];
rz(0.38715036) q[1];
rz(-pi) q[2];
rz(-0.79521631) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(-0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(0.14848407) q[0];
x q[1];
rz(2.4321062) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(-0.62765861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1670926) q[1];
sx q[1];
rz(-1.316861) q[1];
sx q[1];
rz(2.7212935) q[1];
rz(-pi) q[2];
rz(-1.4278533) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(1.8247719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41985837) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255071) q[0];
sx q[0];
rz(-1.320991) q[0];
sx q[0];
rz(-0.020629701) q[0];
rz(-pi) q[1];
rz(1.8823207) q[2];
sx q[2];
rz(-1.4853962) q[2];
sx q[2];
rz(2.4740919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64297134) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(0.69795124) q[1];
rz(-pi) q[2];
rz(1.8886391) q[3];
sx q[3];
rz(-2.7405973) q[3];
sx q[3];
rz(-0.81438118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(0.88090849) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8857408) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(-2.6141502) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23668004) q[2];
sx q[2];
rz(-1.5152021) q[2];
sx q[2];
rz(2.6829164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.199898) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(-0.56337507) q[1];
x q[2];
rz(-0.36066182) q[3];
sx q[3];
rz(-0.70913991) q[3];
sx q[3];
rz(3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2252786) q[0];
sx q[0];
rz(-1.7597223) q[0];
sx q[0];
rz(-0.96154763) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0448858) q[2];
sx q[2];
rz(-0.57069639) q[2];
sx q[2];
rz(1.8622461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67152126) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(2.3684711) q[1];
rz(-pi) q[2];
rz(0.14857265) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(-3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(-1.0587143) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(-0.44021846) q[0];
rz(2.5541359) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(2.6641252) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7212972) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(-0.23423127) q[1];
rz(0.70763564) q[3];
sx q[3];
rz(-1.4228627) q[3];
sx q[3];
rz(-1.1843475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.6301427) q[0];
rz(1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232988) q[0];
sx q[0];
rz(-2.5677486) q[0];
sx q[0];
rz(-2.7159575) q[0];
x q[1];
rz(0.62692554) q[2];
sx q[2];
rz(-1.727384) q[2];
sx q[2];
rz(-0.92948929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48982606) q[1];
sx q[1];
rz(-0.090432743) q[1];
sx q[1];
rz(-2.138278) q[1];
x q[2];
rz(0.82960415) q[3];
sx q[3];
rz(-1.7332352) q[3];
sx q[3];
rz(0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4293489) q[0];
sx q[0];
rz(-1.1375543) q[0];
sx q[0];
rz(0.54746763) q[0];
rz(-pi) q[1];
x q[1];
rz(1.278644) q[2];
sx q[2];
rz(-0.38049618) q[2];
sx q[2];
rz(-0.7064864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(0.559505) q[1];
rz(0.29406677) q[3];
sx q[3];
rz(-1.2122279) q[3];
sx q[3];
rz(2.8420574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8326571) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(2.1946857) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614721) q[0];
sx q[0];
rz(-1.2456018) q[0];
sx q[0];
rz(1.8878216) q[0];
rz(-pi) q[1];
rz(1.5309179) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(3.1181042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0526035) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(-0.25119541) q[1];
rz(-pi) q[2];
rz(-0.19091786) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-2.4972829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(-2.9677662) q[0];
rz(-pi) q[1];
rz(2.9115453) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(-2.7237797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.555968) q[1];
sx q[1];
rz(-1.2871736) q[1];
sx q[1];
rz(2.7908299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.321407) q[3];
sx q[3];
rz(-0.60884848) q[3];
sx q[3];
rz(-1.604515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(2.1951998) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(-0.89623981) q[3];
sx q[3];
rz(-1.8825718) q[3];
sx q[3];
rz(2.752302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
