OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(1.877797) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11557134) q[0];
sx q[0];
rz(-1.3260256) q[0];
sx q[0];
rz(0.90306247) q[0];
rz(0.99743263) q[2];
sx q[2];
rz(-1.9321529) q[2];
sx q[2];
rz(-2.685355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3448779) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(0.085436324) q[1];
rz(-pi) q[2];
rz(0.24407152) q[3];
sx q[3];
rz(-1.1435978) q[3];
sx q[3];
rz(2.0352767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(1.3686352) q[2];
rz(2.9030419) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(-0.11428782) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7411165) q[0];
sx q[0];
rz(-2.0023161) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(3.053983) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(-1.5709343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.621839) q[0];
sx q[0];
rz(-1.6770419) q[0];
sx q[0];
rz(2.6563717) q[0];
rz(-pi) q[1];
rz(-0.3377487) q[2];
sx q[2];
rz(-2.9523473) q[2];
sx q[2];
rz(-2.7836329) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9894692) q[1];
sx q[1];
rz(-1.996576) q[1];
sx q[1];
rz(-2.2683737) q[1];
rz(-2.0474954) q[3];
sx q[3];
rz(-2.3040651) q[3];
sx q[3];
rz(2.1708716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7972083) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(-0.79948419) q[3];
sx q[3];
rz(-1.5461642) q[3];
sx q[3];
rz(2.1605087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.28476533) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(-3.0322266) q[0];
rz(0.60375396) q[1];
sx q[1];
rz(-1.3394638) q[1];
sx q[1];
rz(-1.0341136) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4277735) q[0];
sx q[0];
rz(-1.3610916) q[0];
sx q[0];
rz(0.13596491) q[0];
rz(-pi) q[1];
rz(3.0523006) q[2];
sx q[2];
rz(-1.6411601) q[2];
sx q[2];
rz(-2.8825837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8173966) q[1];
sx q[1];
rz(-1.4651555) q[1];
sx q[1];
rz(0.28075851) q[1];
rz(0.77165551) q[3];
sx q[3];
rz(-1.4691741) q[3];
sx q[3];
rz(1.1751428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6094531) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(2.3405781) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48701778) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(2.6825478) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(-2.3150516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61454489) q[0];
sx q[0];
rz(-1.9226719) q[0];
sx q[0];
rz(1.1220758) q[0];
x q[1];
rz(-3.0439348) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(2.1922534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7956808) q[1];
sx q[1];
rz(-1.7670146) q[1];
sx q[1];
rz(2.1225342) q[1];
rz(-pi) q[2];
rz(1.6891805) q[3];
sx q[3];
rz(-1.9031798) q[3];
sx q[3];
rz(2.3207302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7894342) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-1.1226783) q[3];
sx q[3];
rz(-1.2238945) q[3];
sx q[3];
rz(0.96021715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66626755) q[0];
sx q[0];
rz(-2.1551977) q[0];
sx q[0];
rz(0.112003) q[0];
rz(0.22459596) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-0.78132838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3749237) q[0];
sx q[0];
rz(-1.5974853) q[0];
sx q[0];
rz(1.1666537) q[0];
x q[1];
rz(-2.8113643) q[2];
sx q[2];
rz(-1.4817186) q[2];
sx q[2];
rz(1.1116127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0722149) q[1];
sx q[1];
rz(-0.88000127) q[1];
sx q[1];
rz(-2.1363972) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9797793) q[3];
sx q[3];
rz(-2.1573632) q[3];
sx q[3];
rz(2.2590421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3096699) q[2];
sx q[2];
rz(-2.2396542) q[2];
sx q[2];
rz(2.208948) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340853) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.750741) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(-1.5501685) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0858147) q[0];
sx q[0];
rz(-1.5356488) q[0];
sx q[0];
rz(-0.9279535) q[0];
x q[1];
rz(-1.6076902) q[2];
sx q[2];
rz(-2.0301182) q[2];
sx q[2];
rz(-1.5609891) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6350123) q[1];
sx q[1];
rz(-1.4476579) q[1];
sx q[1];
rz(-1.413373) q[1];
x q[2];
rz(-1.1892631) q[3];
sx q[3];
rz(-1.23151) q[3];
sx q[3];
rz(1.8041704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3127689) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(-2.17365) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62436002) q[0];
sx q[0];
rz(-0.42651287) q[0];
sx q[0];
rz(2.970001) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(-1.4412057) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489451) q[0];
sx q[0];
rz(-1.2496557) q[0];
sx q[0];
rz(-1.4399066) q[0];
rz(2.9459459) q[2];
sx q[2];
rz(-0.52426978) q[2];
sx q[2];
rz(-1.9877246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4525057) q[1];
sx q[1];
rz(-1.2156665) q[1];
sx q[1];
rz(1.0853115) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1342742) q[3];
sx q[3];
rz(-1.451056) q[3];
sx q[3];
rz(0.086825018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8695996) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(2.5970411) q[2];
rz(-0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(0.47206363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.8498103) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(0.52919069) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(-1.1189438) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3854348) q[0];
sx q[0];
rz(-1.5498112) q[0];
sx q[0];
rz(0.051647112) q[0];
x q[1];
rz(2.697102) q[2];
sx q[2];
rz(-0.40676446) q[2];
sx q[2];
rz(3.0009746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3333867) q[1];
sx q[1];
rz(-1.1115326) q[1];
sx q[1];
rz(2.4833268) q[1];
x q[2];
rz(0.20967926) q[3];
sx q[3];
rz(-2.0884313) q[3];
sx q[3];
rz(2.8594494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90073663) q[2];
sx q[2];
rz(-0.46515981) q[2];
sx q[2];
rz(-0.71211234) q[2];
rz(1.5322878) q[3];
sx q[3];
rz(-0.94523793) q[3];
sx q[3];
rz(1.3473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79494548) q[0];
sx q[0];
rz(-1.1279339) q[0];
sx q[0];
rz(-2.0704863) q[0];
rz(-1.935293) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(1.6546904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95222118) q[0];
sx q[0];
rz(-0.69714386) q[0];
sx q[0];
rz(2.8482262) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63196147) q[2];
sx q[2];
rz(-0.8425396) q[2];
sx q[2];
rz(-1.8441083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35569977) q[1];
sx q[1];
rz(-1.3972056) q[1];
sx q[1];
rz(-2.8462571) q[1];
rz(-pi) q[2];
rz(-2.1110299) q[3];
sx q[3];
rz(-0.5117473) q[3];
sx q[3];
rz(2.668021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(3.0637975) q[2];
rz(-2.9142694) q[3];
sx q[3];
rz(-2.0096171) q[3];
sx q[3];
rz(2.0859065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530497) q[0];
sx q[0];
rz(-2.396614) q[0];
sx q[0];
rz(2.7689834) q[0];
rz(-0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(-2.2850697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055182) q[0];
sx q[0];
rz(-1.627632) q[0];
sx q[0];
rz(-1.6147805) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7046518) q[2];
sx q[2];
rz(-1.6758462) q[2];
sx q[2];
rz(-1.5488226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0792993) q[1];
sx q[1];
rz(-0.41626272) q[1];
sx q[1];
rz(1.2757311) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3395446) q[3];
sx q[3];
rz(-1.060876) q[3];
sx q[3];
rz(2.0039441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7591758) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(-2.5860533) q[2];
rz(0.38604745) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(0.66421318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0059218) q[0];
sx q[0];
rz(-1.4501403) q[0];
sx q[0];
rz(-1.8474664) q[0];
rz(-0.80264965) q[1];
sx q[1];
rz(-1.6812656) q[1];
sx q[1];
rz(0.38801286) q[1];
rz(-2.7886709) q[2];
sx q[2];
rz(-1.4581994) q[2];
sx q[2];
rz(-1.1117473) q[2];
rz(-2.3014746) q[3];
sx q[3];
rz(-1.6230604) q[3];
sx q[3];
rz(-0.75480672) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
