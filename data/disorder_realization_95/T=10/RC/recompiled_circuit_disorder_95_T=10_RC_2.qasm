OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449887) q[0];
sx q[0];
rz(-0.94192266) q[0];
sx q[0];
rz(-3.0552342) q[0];
x q[1];
rz(2.2034982) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(3.0245568) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.861046) q[1];
sx q[1];
rz(-1.7409054) q[1];
sx q[1];
rz(2.4982846) q[1];
x q[2];
rz(-0.1944794) q[3];
sx q[3];
rz(-0.76768657) q[3];
sx q[3];
rz(-0.0038113468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(2.2251341) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-0.8786456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4247596) q[0];
sx q[0];
rz(-2.3905907) q[0];
sx q[0];
rz(-2.1483833) q[0];
rz(-pi) q[1];
rz(0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-0.49371142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2616927) q[1];
sx q[1];
rz(-2.7163134) q[1];
sx q[1];
rz(0.17875032) q[1];
x q[2];
rz(-0.068275498) q[3];
sx q[3];
rz(-1.2618511) q[3];
sx q[3];
rz(2.6451049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8615222) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.52755255) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(-1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(-1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3504346) q[0];
sx q[0];
rz(-0.8236304) q[0];
sx q[0];
rz(-0.50823786) q[0];
rz(-pi) q[1];
rz(3.0776575) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(2.8901697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(-0.73168879) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(-0.63347647) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(-1.7184005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0013106) q[1];
sx q[1];
rz(-1.9754859) q[1];
sx q[1];
rz(-2.6005122) q[1];
rz(-pi) q[2];
rz(1.5192401) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(1.6980905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(0.25203618) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900433) q[0];
sx q[0];
rz(-0.94928375) q[0];
sx q[0];
rz(3.050699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0257452) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(-2.2124706) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7733113) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.7721304) q[1];
x q[2];
rz(-2.4162021) q[3];
sx q[3];
rz(-1.8728421) q[3];
sx q[3];
rz(-0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(0.71470913) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(0.2063624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9128742) q[0];
sx q[0];
rz(-1.657907) q[0];
sx q[0];
rz(-0.74761439) q[0];
x q[1];
rz(-2.2588737) q[2];
sx q[2];
rz(-1.025891) q[2];
sx q[2];
rz(-2.8628652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5620835) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(-3.0444006) q[1];
x q[2];
rz(-1.0422802) q[3];
sx q[3];
rz(-0.37502608) q[3];
sx q[3];
rz(-0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21748397) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(0.80362513) q[0];
x q[1];
rz(2.6518875) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(0.44039886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8163029) q[1];
sx q[1];
rz(-1.5029229) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(2.8190523) q[3];
sx q[3];
rz(-0.99223677) q[3];
sx q[3];
rz(1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3367735) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.221009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2894665) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(-0.37129398) q[0];
rz(2.3983725) q[2];
sx q[2];
rz(-2.9033702) q[2];
sx q[2];
rz(-1.5526349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80291058) q[1];
sx q[1];
rz(-1.7065587) q[1];
sx q[1];
rz(-1.6822097) q[1];
rz(-pi) q[2];
rz(2.8052748) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28875479) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98709479) q[0];
sx q[0];
rz(-2.8821766) q[0];
sx q[0];
rz(-2.3214066) q[0];
x q[1];
rz(-2.8081886) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(2.8508027) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3667664) q[1];
sx q[1];
rz(-2.3407196) q[1];
sx q[1];
rz(1.1714539) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3689234) q[3];
sx q[3];
rz(-1.8608421) q[3];
sx q[3];
rz(1.9162852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126497) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(1.2055231) q[0];
x q[1];
rz(0.26009772) q[2];
sx q[2];
rz(-1.1220699) q[2];
sx q[2];
rz(-1.5425494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1117489) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(2.7010121) q[1];
rz(2.3481028) q[3];
sx q[3];
rz(-1.3440545) q[3];
sx q[3];
rz(1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-0.5207516) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(-0.76277914) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(2.5557774) q[3];
sx q[3];
rz(-1.9184434) q[3];
sx q[3];
rz(-0.8003269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
