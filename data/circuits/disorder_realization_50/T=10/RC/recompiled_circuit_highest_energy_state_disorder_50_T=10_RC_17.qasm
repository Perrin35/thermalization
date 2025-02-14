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
rz(0.4515689) q[0];
sx q[0];
rz(-1.7234001) q[0];
sx q[0];
rz(0.42887846) q[0];
rz(-0.14313993) q[1];
sx q[1];
rz(-2.0506471) q[1];
sx q[1];
rz(-3.0615038) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9432491) q[0];
sx q[0];
rz(-0.08028537) q[0];
sx q[0];
rz(-0.16585089) q[0];
rz(-3.0386904) q[2];
sx q[2];
rz(-1.9299302) q[2];
sx q[2];
rz(-2.5996836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7987207) q[1];
sx q[1];
rz(-1.8340543) q[1];
sx q[1];
rz(2.2895471) q[1];
rz(-pi) q[2];
rz(2.56404) q[3];
sx q[3];
rz(-1.7087382) q[3];
sx q[3];
rz(2.4393314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.074778883) q[2];
sx q[2];
rz(-2.0902233) q[2];
sx q[2];
rz(1.7195513) q[2];
rz(1.413013) q[3];
sx q[3];
rz(-1.6001817) q[3];
sx q[3];
rz(-1.448267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(-2.3894892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1307169) q[0];
sx q[0];
rz(-0.97144043) q[0];
sx q[0];
rz(-0.54979558) q[0];
rz(-3.0235748) q[2];
sx q[2];
rz(-2.9213597) q[2];
sx q[2];
rz(1.3853467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32622329) q[1];
sx q[1];
rz(-1.0462251) q[1];
sx q[1];
rz(-0.12614819) q[1];
rz(-2.8798298) q[3];
sx q[3];
rz(-1.7883375) q[3];
sx q[3];
rz(-2.1037773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.077945262) q[2];
sx q[2];
rz(-0.91219488) q[2];
sx q[2];
rz(-0.65703854) q[2];
rz(-0.71197236) q[3];
sx q[3];
rz(-2.4896121) q[3];
sx q[3];
rz(2.3967192) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7065358) q[0];
sx q[0];
rz(-2.5223795) q[0];
sx q[0];
rz(1.0414498) q[0];
rz(-2.7667747) q[1];
sx q[1];
rz(-1.638214) q[1];
sx q[1];
rz(2.5846438) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0335046) q[0];
sx q[0];
rz(-1.9226941) q[0];
sx q[0];
rz(-2.3826588) q[0];
x q[1];
rz(-3.0663112) q[2];
sx q[2];
rz(-2.0213375) q[2];
sx q[2];
rz(2.0496429) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3279325) q[1];
sx q[1];
rz(-2.2371462) q[1];
sx q[1];
rz(1.440669) q[1];
rz(-pi) q[2];
rz(1.3339504) q[3];
sx q[3];
rz(-0.93333731) q[3];
sx q[3];
rz(-1.4887631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53691429) q[2];
sx q[2];
rz(-0.92427173) q[2];
sx q[2];
rz(-0.43517819) q[2];
rz(-0.52470454) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(3.0068126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404469) q[0];
sx q[0];
rz(-1.0643767) q[0];
sx q[0];
rz(1.5144298) q[0];
rz(1.0487522) q[1];
sx q[1];
rz(-2.2188413) q[1];
sx q[1];
rz(-1.4357766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.993282) q[0];
sx q[0];
rz(-1.7459501) q[0];
sx q[0];
rz(-1.6143285) q[0];
rz(-pi) q[1];
rz(-0.54331358) q[2];
sx q[2];
rz(-1.8738973) q[2];
sx q[2];
rz(1.4599279) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32325778) q[1];
sx q[1];
rz(-2.0916953) q[1];
sx q[1];
rz(-0.073831115) q[1];
rz(-pi) q[2];
rz(0.9463598) q[3];
sx q[3];
rz(-2.588996) q[3];
sx q[3];
rz(2.432807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7413896) q[2];
sx q[2];
rz(-1.3578537) q[2];
sx q[2];
rz(0.025156585) q[2];
rz(-1.4870421) q[3];
sx q[3];
rz(-0.39792037) q[3];
sx q[3];
rz(-0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095489278) q[0];
sx q[0];
rz(-2.4920721) q[0];
sx q[0];
rz(0.15897861) q[0];
rz(1.1051296) q[1];
sx q[1];
rz(-0.72560328) q[1];
sx q[1];
rz(-0.22430688) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0889322) q[0];
sx q[0];
rz(-1.1692246) q[0];
sx q[0];
rz(-2.5395655) q[0];
rz(1.5184899) q[2];
sx q[2];
rz(-2.5807267) q[2];
sx q[2];
rz(0.069151783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4558002) q[1];
sx q[1];
rz(-0.53736254) q[1];
sx q[1];
rz(-1.3440029) q[1];
rz(2.5644825) q[3];
sx q[3];
rz(-1.464586) q[3];
sx q[3];
rz(-2.4662227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20532456) q[2];
sx q[2];
rz(-1.779413) q[2];
sx q[2];
rz(2.3929907) q[2];
rz(0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1614138) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(0.73032105) q[0];
rz(-1.6397363) q[1];
sx q[1];
rz(-1.7514936) q[1];
sx q[1];
rz(-0.23434848) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47539179) q[0];
sx q[0];
rz(-2.0115543) q[0];
sx q[0];
rz(-0.73385629) q[0];
rz(-pi) q[1];
rz(2.2147398) q[2];
sx q[2];
rz(-2.0704464) q[2];
sx q[2];
rz(0.32678963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40962666) q[1];
sx q[1];
rz(-1.6183369) q[1];
sx q[1];
rz(-2.9788245) q[1];
x q[2];
rz(1.430351) q[3];
sx q[3];
rz(-1.5101119) q[3];
sx q[3];
rz(1.3827326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.948287) q[2];
sx q[2];
rz(-0.81393465) q[2];
sx q[2];
rz(1.2813655) q[2];
rz(2.5365601) q[3];
sx q[3];
rz(-1.0993967) q[3];
sx q[3];
rz(1.9090778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.989711) q[0];
sx q[0];
rz(-1.2622156) q[0];
sx q[0];
rz(0.61304098) q[0];
rz(-0.68060654) q[1];
sx q[1];
rz(-2.0304408) q[1];
sx q[1];
rz(1.8162762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2277057) q[0];
sx q[0];
rz(-1.837758) q[0];
sx q[0];
rz(-0.56446979) q[0];
rz(-pi) q[1];
rz(2.7737977) q[2];
sx q[2];
rz(-2.9136474) q[2];
sx q[2];
rz(-1.6942555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5449574) q[1];
sx q[1];
rz(-1.4442634) q[1];
sx q[1];
rz(0.56565779) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9716827) q[3];
sx q[3];
rz(-1.0260149) q[3];
sx q[3];
rz(2.0044897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94720542) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(0.80805937) q[2];
rz(-2.4957073) q[3];
sx q[3];
rz(-0.26791993) q[3];
sx q[3];
rz(2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5757669) q[0];
sx q[0];
rz(-2.4149826) q[0];
sx q[0];
rz(0.44276825) q[0];
rz(2.5495461) q[1];
sx q[1];
rz(-0.7213842) q[1];
sx q[1];
rz(-2.9659042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73978776) q[0];
sx q[0];
rz(-1.8462088) q[0];
sx q[0];
rz(2.7729183) q[0];
rz(-2.6900411) q[2];
sx q[2];
rz(-1.5349261) q[2];
sx q[2];
rz(-0.6605688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67461102) q[1];
sx q[1];
rz(-1.5276485) q[1];
sx q[1];
rz(-0.64137913) q[1];
rz(-pi) q[2];
rz(2.2288228) q[3];
sx q[3];
rz(-2.4173749) q[3];
sx q[3];
rz(0.66373827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(-2.234484) q[2];
rz(-0.90726566) q[3];
sx q[3];
rz(-1.3943358) q[3];
sx q[3];
rz(-3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9573145) q[0];
sx q[0];
rz(-0.90183455) q[0];
sx q[0];
rz(2.9515475) q[0];
rz(1.5826781) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(-1.83439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060576749) q[0];
sx q[0];
rz(-1.2071484) q[0];
sx q[0];
rz(-2.9878997) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8030012) q[2];
sx q[2];
rz(-2.2420852) q[2];
sx q[2];
rz(0.012635144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7293882) q[1];
sx q[1];
rz(-1.6784759) q[1];
sx q[1];
rz(3.1086224) q[1];
rz(-pi) q[2];
rz(0.43453026) q[3];
sx q[3];
rz(-1.4659681) q[3];
sx q[3];
rz(-1.9275582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.068453161) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(-2.9404409) q[2];
rz(2.1189832) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(-1.6020017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1744743) q[0];
sx q[0];
rz(-2.0956464) q[0];
sx q[0];
rz(2.2651267) q[0];
rz(-0.62249741) q[1];
sx q[1];
rz(-2.9298156) q[1];
sx q[1];
rz(-0.34271398) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22607813) q[0];
sx q[0];
rz(-2.5097339) q[0];
sx q[0];
rz(1.4591239) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12212716) q[2];
sx q[2];
rz(-1.9052278) q[2];
sx q[2];
rz(1.8129406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6500998) q[1];
sx q[1];
rz(-2.0180185) q[1];
sx q[1];
rz(2.4504721) q[1];
x q[2];
rz(2.5779419) q[3];
sx q[3];
rz(-2.6915801) q[3];
sx q[3];
rz(-1.2695509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9524625) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-0.44378898) q[2];
rz(-2.022838) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(-0.55383033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8583869) q[0];
sx q[0];
rz(-1.5844185) q[0];
sx q[0];
rz(1.6208741) q[0];
rz(-3.0781147) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-0.061894682) q[2];
sx q[2];
rz(-1.2409004) q[2];
sx q[2];
rz(0.093766669) q[2];
rz(-1.1984428) q[3];
sx q[3];
rz(-0.59881864) q[3];
sx q[3];
rz(-1.3074085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
