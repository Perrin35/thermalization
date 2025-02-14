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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(2.5084578) q[1];
sx q[1];
rz(8.8517744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0664803) q[0];
sx q[0];
rz(-1.5085992) q[0];
sx q[0];
rz(0.94554995) q[0];
rz(-pi) q[1];
rz(-1.139976) q[2];
sx q[2];
rz(-0.9316906) q[2];
sx q[2];
rz(-0.72102816) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1243368) q[1];
sx q[1];
rz(-1.5005497) q[1];
sx q[1];
rz(-2.6439366) q[1];
rz(-pi) q[2];
rz(-0.12270452) q[3];
sx q[3];
rz(-2.4413681) q[3];
sx q[3];
rz(2.6168857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28213349) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(2.0852883) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2317155) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6352954) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(0.70962972) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0906648) q[2];
sx q[2];
rz(-1.6279334) q[2];
sx q[2];
rz(0.47540755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9641407) q[1];
sx q[1];
rz(-2.3712082) q[1];
sx q[1];
rz(0.48078534) q[1];
rz(-2.7892465) q[3];
sx q[3];
rz(-2.5735613) q[3];
sx q[3];
rz(-1.2888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(0.44542584) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(0.87944952) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875558) q[0];
sx q[0];
rz(-1.7374514) q[0];
sx q[0];
rz(2.119407) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1398323) q[2];
sx q[2];
rz(-1.941701) q[2];
sx q[2];
rz(-2.6628464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5048557) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(1.4192753) q[1];
x q[2];
rz(1.2233569) q[3];
sx q[3];
rz(-1.8104324) q[3];
sx q[3];
rz(-2.6289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(1.5039911) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(-0.98201069) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0729436) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(1.9085931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017576305) q[0];
sx q[0];
rz(-1.4498382) q[0];
sx q[0];
rz(1.308754) q[0];
x q[1];
rz(2.2188051) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(3.0344935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0291042) q[1];
sx q[1];
rz(-1.8524516) q[1];
sx q[1];
rz(-2.2206578) q[1];
x q[2];
rz(-3.0844131) q[3];
sx q[3];
rz(-2.8897277) q[3];
sx q[3];
rz(-2.4185023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.5114463) q[2];
rz(1.1394507) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(-0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646506) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(1.1530217) q[0];
rz(2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(2.8797454) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7196822) q[0];
sx q[0];
rz(-1.6043264) q[0];
sx q[0];
rz(-1.6543341) q[0];
x q[1];
rz(-2.0551497) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(-0.5400368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56526977) q[1];
sx q[1];
rz(-2.5082631) q[1];
sx q[1];
rz(2.759658) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7747986) q[3];
sx q[3];
rz(-2.3424934) q[3];
sx q[3];
rz(0.71281709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61949817) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(-1.7363133) q[2];
rz(-0.74553472) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(-2.9702759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650314) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(1.879346) q[0];
x q[1];
rz(-1.2662884) q[2];
sx q[2];
rz(-2.300548) q[2];
sx q[2];
rz(-0.27308057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3503987) q[1];
sx q[1];
rz(-1.5916675) q[1];
sx q[1];
rz(0.69857614) q[1];
rz(-0.71368696) q[3];
sx q[3];
rz(-2.1296576) q[3];
sx q[3];
rz(1.9969459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0773086) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-2.1997931) q[2];
rz(-1.1549548) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44523859) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(0.44662961) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-0.68797025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94243455) q[0];
sx q[0];
rz(-0.73478886) q[0];
sx q[0];
rz(-0.25811974) q[0];
x q[1];
rz(0.61882682) q[2];
sx q[2];
rz(-1.8558673) q[2];
sx q[2];
rz(-3.0322591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0258642) q[1];
sx q[1];
rz(-2.0021571) q[1];
sx q[1];
rz(-0.23484767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3456081) q[3];
sx q[3];
rz(-2.2777657) q[3];
sx q[3];
rz(0.10315264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61116162) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.2215349) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68194836) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(-1.4247165) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(-0.43509126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1662707) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(-0.50384371) q[0];
rz(1.7516364) q[2];
sx q[2];
rz(-2.2308084) q[2];
sx q[2];
rz(-1.7969799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96174091) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(0.95203103) q[1];
rz(-pi) q[2];
rz(-3.1396542) q[3];
sx q[3];
rz(-2.7172305) q[3];
sx q[3];
rz(-2.3821609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6626176) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(-0.62620658) q[2];
rz(-0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5597252) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.8064226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30935449) q[0];
sx q[0];
rz(-1.0296427) q[0];
sx q[0];
rz(-0.65017976) q[0];
rz(1.0112052) q[2];
sx q[2];
rz(-1.8768684) q[2];
sx q[2];
rz(1.0572421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9250592) q[1];
sx q[1];
rz(-1.0581985) q[1];
sx q[1];
rz(2.53043) q[1];
rz(0.068633462) q[3];
sx q[3];
rz(-1.7083999) q[3];
sx q[3];
rz(0.28561628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6173031) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(0.35150251) q[2];
rz(1.3737804) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(2.3824298) q[0];
rz(1.0836541) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(-2.0738475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6747492) q[0];
sx q[0];
rz(-2.3567794) q[0];
sx q[0];
rz(2.270257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8482089) q[2];
sx q[2];
rz(-1.2870064) q[2];
sx q[2];
rz(-2.0744155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0615649) q[1];
sx q[1];
rz(-2.1127709) q[1];
sx q[1];
rz(-0.47596495) q[1];
x q[2];
rz(-2.8078733) q[3];
sx q[3];
rz(-0.90147831) q[3];
sx q[3];
rz(3.1143318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(-0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44337153) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(-1.632985) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(0.40661033) q[2];
sx q[2];
rz(-0.92401531) q[2];
sx q[2];
rz(0.32377908) q[2];
rz(-2.0785594) q[3];
sx q[3];
rz(-0.31217839) q[3];
sx q[3];
rz(0.97806539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
