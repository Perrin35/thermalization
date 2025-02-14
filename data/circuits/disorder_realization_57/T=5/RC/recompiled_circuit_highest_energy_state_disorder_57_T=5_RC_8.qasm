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
rz(-1.3402101) q[0];
sx q[0];
rz(3.4184472) q[0];
sx q[0];
rz(10.463538) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2622803) q[0];
sx q[0];
rz(-1.6769772) q[0];
sx q[0];
rz(-0.72632974) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38484599) q[2];
sx q[2];
rz(-0.59078465) q[2];
sx q[2];
rz(-0.42921517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8326679) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(-0.69242386) q[1];
rz(2.9800426) q[3];
sx q[3];
rz(-1.5388425) q[3];
sx q[3];
rz(0.7940426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1699528) q[2];
sx q[2];
rz(-1.4687186) q[2];
sx q[2];
rz(-1.8876342) q[2];
rz(1.5818671) q[3];
sx q[3];
rz(-2.4032148) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935683) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(-0.76876202) q[0];
rz(-1.7747152) q[1];
sx q[1];
rz(-1.3221075) q[1];
sx q[1];
rz(-1.0381402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0739792) q[0];
sx q[0];
rz(-1.5609976) q[0];
sx q[0];
rz(0.0086179535) q[0];
rz(-1.9351472) q[2];
sx q[2];
rz(-0.88308217) q[2];
sx q[2];
rz(1.5503413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2778216) q[1];
sx q[1];
rz(-1.9367332) q[1];
sx q[1];
rz(-0.65524958) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6065381) q[3];
sx q[3];
rz(-2.6527191) q[3];
sx q[3];
rz(1.3066178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(0.25343728) q[3];
sx q[3];
rz(-0.35496747) q[3];
sx q[3];
rz(-0.36111116) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(-2.2542727) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-2.2128426) q[1];
sx q[1];
rz(-1.1393772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47592012) q[0];
sx q[0];
rz(-0.3753271) q[0];
sx q[0];
rz(0.39552839) q[0];
rz(-pi) q[1];
rz(1.0408632) q[2];
sx q[2];
rz(-1.8609253) q[2];
sx q[2];
rz(-2.0255956) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92147747) q[1];
sx q[1];
rz(-2.8954828) q[1];
sx q[1];
rz(-3.0819478) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5185131) q[3];
sx q[3];
rz(-2.0889335) q[3];
sx q[3];
rz(0.62925807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1246216) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-0.43928453) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.144217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-0.93337494) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(-3.0237517) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(1.3062564) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7457434) q[0];
sx q[0];
rz(-1.1538528) q[0];
sx q[0];
rz(-0.31503079) q[0];
rz(0.47331402) q[2];
sx q[2];
rz(-2.4243948) q[2];
sx q[2];
rz(-2.0462478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0769121) q[1];
sx q[1];
rz(-0.91935989) q[1];
sx q[1];
rz(-0.76131911) q[1];
rz(-1.4545069) q[3];
sx q[3];
rz(-0.57580417) q[3];
sx q[3];
rz(0.35279122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.304504) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(2.9869249) q[2];
rz(1.539544) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(-0.60607564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347572) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(-0.71594816) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(-0.11944019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855921) q[0];
sx q[0];
rz(-1.553926) q[0];
sx q[0];
rz(-1.5149679) q[0];
rz(-pi) q[1];
rz(-0.44942707) q[2];
sx q[2];
rz(-1.4327421) q[2];
sx q[2];
rz(-1.4269478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24158289) q[1];
sx q[1];
rz(-1.1092343) q[1];
sx q[1];
rz(-0.66630967) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64924134) q[3];
sx q[3];
rz(-0.63707817) q[3];
sx q[3];
rz(1.0604881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9014088) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(-3.0805947) q[2];
rz(0.048246233) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(0.72859305) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401684) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(0.4785969) q[1];
sx q[1];
rz(-2.5391948) q[1];
sx q[1];
rz(-1.7344249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76001747) q[0];
sx q[0];
rz(-1.572346) q[0];
sx q[0];
rz(-3.0880307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2931104) q[2];
sx q[2];
rz(-1.3568078) q[2];
sx q[2];
rz(3.1297562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4987405) q[1];
sx q[1];
rz(-1.3428709) q[1];
sx q[1];
rz(3.0251316) q[1];
rz(-2.9286372) q[3];
sx q[3];
rz(-1.3939438) q[3];
sx q[3];
rz(0.18421728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5799134) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(1.178406) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(-1.8572846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9310164) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(-1.3354906) q[0];
rz(1.3212851) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(0.30805045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9481407) q[0];
sx q[0];
rz(-1.3592255) q[0];
sx q[0];
rz(2.7443462) q[0];
rz(-pi) q[1];
rz(2.5045583) q[2];
sx q[2];
rz(-2.5278628) q[2];
sx q[2];
rz(-1.5368652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1738893) q[1];
sx q[1];
rz(-2.1427665) q[1];
sx q[1];
rz(0.58764761) q[1];
x q[2];
rz(-3.0940476) q[3];
sx q[3];
rz(-1.9564087) q[3];
sx q[3];
rz(-1.6507208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41219741) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3103631) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(2.6112153) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(2.4200965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24151267) q[0];
sx q[0];
rz(-2.3519313) q[0];
sx q[0];
rz(-2.4898743) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8558735) q[2];
sx q[2];
rz(-2.6556394) q[2];
sx q[2];
rz(3.1348117) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7466654) q[1];
sx q[1];
rz(-0.29005602) q[1];
sx q[1];
rz(-2.2941089) q[1];
rz(-pi) q[2];
rz(-2.3590478) q[3];
sx q[3];
rz(-1.492332) q[3];
sx q[3];
rz(-2.9786547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.4777769) q[2];
rz(-1.1635228) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(-1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32996938) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(2.5323618) q[0];
rz(-1.5513264) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824175) q[0];
sx q[0];
rz(-1.32138) q[0];
sx q[0];
rz(-1.9391869) q[0];
x q[1];
rz(0.48522075) q[2];
sx q[2];
rz(-1.5169797) q[2];
sx q[2];
rz(-2.9973928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(-1.9202597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2942469) q[3];
sx q[3];
rz(-1.0718126) q[3];
sx q[3];
rz(-0.33667281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41813254) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708165) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(0.01195512) q[0];
rz(1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(0.59648046) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76145455) q[0];
sx q[0];
rz(-1.589847) q[0];
sx q[0];
rz(0.46940243) q[0];
rz(-1.8305186) q[2];
sx q[2];
rz(-1.2190281) q[2];
sx q[2];
rz(-3.0992257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0833019) q[1];
sx q[1];
rz(-1.4920248) q[1];
sx q[1];
rz(1.0942671) q[1];
rz(0.47816737) q[3];
sx q[3];
rz(-0.33278012) q[3];
sx q[3];
rz(-1.5733583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.163588) q[2];
sx q[2];
rz(-0.96077335) q[2];
sx q[2];
rz(-2.9445924) q[2];
rz(1.152285) q[3];
sx q[3];
rz(-2.9983493) q[3];
sx q[3];
rz(-2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(0.20881431) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(0.88125689) q[2];
sx q[2];
rz(-1.643558) q[2];
sx q[2];
rz(0.087842077) q[2];
rz(2.1590334) q[3];
sx q[3];
rz(-2.0053902) q[3];
sx q[3];
rz(-2.7872661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
