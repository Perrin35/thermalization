OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(4.0806169) q[0];
sx q[0];
rz(9.4299849) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(-1.9519238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.386728) q[0];
sx q[0];
rz(-2.0804188) q[0];
sx q[0];
rz(1.1763563) q[0];
rz(-2.4742545) q[2];
sx q[2];
rz(-0.22775209) q[2];
sx q[2];
rz(-1.3078794) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78566879) q[1];
sx q[1];
rz(-0.34468109) q[1];
sx q[1];
rz(1.1204526) q[1];
x q[2];
rz(2.7351904) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(1.3960081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71620119) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(-0.11322583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229867) q[0];
sx q[0];
rz(-1.9959873) q[0];
sx q[0];
rz(3.0990764) q[0];
x q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5511302) q[2];
sx q[2];
rz(2.427223) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37564056) q[1];
sx q[1];
rz(-1.9687708) q[1];
sx q[1];
rz(-2.6938733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23785915) q[3];
sx q[3];
rz(-2.1616462) q[3];
sx q[3];
rz(1.8464586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.018628) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(2.9193027) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(2.2608742) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(3.0139794) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.688296) q[0];
sx q[0];
rz(-1.7444381) q[0];
sx q[0];
rz(2.2542473) q[0];
rz(-pi) q[1];
rz(-0.67655501) q[2];
sx q[2];
rz(-2.1495719) q[2];
sx q[2];
rz(0.81625953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4096654) q[1];
sx q[1];
rz(-1.8095784) q[1];
sx q[1];
rz(-3.0569397) q[1];
x q[2];
rz(2.1979638) q[3];
sx q[3];
rz(-1.3645932) q[3];
sx q[3];
rz(-2.9507153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3601274) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(0.310251) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(-0.37240949) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.303064) q[0];
sx q[0];
rz(-1.4493363) q[0];
sx q[0];
rz(-0.86304201) q[0];
x q[1];
rz(2.0580975) q[2];
sx q[2];
rz(-2.5022025) q[2];
sx q[2];
rz(2.1208178) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.403703) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(1.9562734) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72752556) q[3];
sx q[3];
rz(-1.993506) q[3];
sx q[3];
rz(2.5142575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-2.2327936) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(2.6866384) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(1.1313261) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(2.4647443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880496) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(-0.44102863) q[0];
rz(0.73111515) q[2];
sx q[2];
rz(-1.7321246) q[2];
sx q[2];
rz(2.4246755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10830282) q[1];
sx q[1];
rz(-1.8255594) q[1];
sx q[1];
rz(1.5441896) q[1];
x q[2];
rz(-1.5464209) q[3];
sx q[3];
rz(-1.1773603) q[3];
sx q[3];
rz(0.76373053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.1452902) q[2];
rz(-0.14906135) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(2.3964264) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(0.93313342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720848) q[0];
sx q[0];
rz(-2.5818995) q[0];
sx q[0];
rz(2.1493388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2007347) q[2];
sx q[2];
rz(-2.1413295) q[2];
sx q[2];
rz(2.7119315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17174013) q[1];
sx q[1];
rz(-1.2029359) q[1];
sx q[1];
rz(-0.0024585558) q[1];
rz(-0.17721456) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(2.9110416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(2.8549426) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(-1.6749143) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(-2.1405623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2350378) q[0];
sx q[0];
rz(-0.28309238) q[0];
sx q[0];
rz(1.4366158) q[0];
x q[1];
rz(1.2008576) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(1.5202886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96879362) q[1];
sx q[1];
rz(-1.7263004) q[1];
sx q[1];
rz(1.1771727) q[1];
x q[2];
rz(-1.3592968) q[3];
sx q[3];
rz(-0.930951) q[3];
sx q[3];
rz(2.0169472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(2.5296339) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(0.38988316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0355426) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(-0.57530595) q[0];
x q[1];
rz(0.6673442) q[2];
sx q[2];
rz(-2.5255425) q[2];
sx q[2];
rz(-1.4916071) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65615678) q[1];
sx q[1];
rz(-1.1404783) q[1];
sx q[1];
rz(-2.7780611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.893703) q[3];
sx q[3];
rz(-2.4363323) q[3];
sx q[3];
rz(-0.44912072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-0.53047696) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-1.4488688) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(0.77734787) q[0];
rz(0.84689394) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(0.27639595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8951176) q[0];
sx q[0];
rz(-0.44024375) q[0];
sx q[0];
rz(2.759139) q[0];
x q[1];
rz(0.70119621) q[2];
sx q[2];
rz(-0.29507911) q[2];
sx q[2];
rz(-1.8104749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2497219) q[1];
sx q[1];
rz(-2.8686214) q[1];
sx q[1];
rz(-2.194838) q[1];
x q[2];
rz(1.8180088) q[3];
sx q[3];
rz(-0.97655481) q[3];
sx q[3];
rz(1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(2.7888443) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(-2.1113077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.83537) q[0];
sx q[0];
rz(-1.3754002) q[0];
sx q[0];
rz(0.72619254) q[0];
x q[1];
rz(-0.68600168) q[2];
sx q[2];
rz(-1.8324319) q[2];
sx q[2];
rz(-2.5930282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13223091) q[1];
sx q[1];
rz(-1.8877601) q[1];
sx q[1];
rz(-0.46225458) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9029721) q[3];
sx q[3];
rz(-0.98364753) q[3];
sx q[3];
rz(1.7410994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-3.0497131) q[3];
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
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108903) q[0];
sx q[0];
rz(-1.1477092) q[0];
sx q[0];
rz(-0.64684091) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(-2.7394221) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(-2.4521811) q[3];
sx q[3];
rz(-2.0333615) q[3];
sx q[3];
rz(0.61483308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
