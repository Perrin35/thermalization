OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(-2.5581869) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28461449) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(2.7906228) q[0];
rz(2.5298654) q[2];
sx q[2];
rz(-2.3731542) q[2];
sx q[2];
rz(2.7035463) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(2.288504) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5276018) q[3];
sx q[3];
rz(-0.82026635) q[3];
sx q[3];
rz(-2.299472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.974568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562247) q[0];
sx q[0];
rz(-1.313374) q[0];
sx q[0];
rz(-3.0811937) q[0];
rz(-pi) q[1];
rz(0.24658792) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.2621244) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0547304) q[1];
sx q[1];
rz(-0.80838258) q[1];
sx q[1];
rz(-0.085309172) q[1];
rz(-pi) q[2];
rz(-2.0608276) q[3];
sx q[3];
rz(-0.96066517) q[3];
sx q[3];
rz(-0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(2.1267557) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(-1.048676) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(-0.93333474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.240757) q[1];
sx q[1];
rz(-1.490331) q[1];
sx q[1];
rz(-0.97297538) q[1];
rz(-pi) q[2];
rz(2.2567743) q[3];
sx q[3];
rz(-1.3664477) q[3];
sx q[3];
rz(-2.200978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(0.63582173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0974554) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.9299279) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3859149) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(-2.3787969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(-0.26054392) q[1];
x q[2];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976945) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-2.863046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700096) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(2.5894126) q[0];
x q[1];
rz(-1.6804382) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.7313752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(-2.6890523) q[1];
rz(2.1513125) q[3];
sx q[3];
rz(-0.84458447) q[3];
sx q[3];
rz(2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(0.0090573514) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6271583) q[0];
sx q[0];
rz(-1.2540199) q[0];
sx q[0];
rz(0.47382521) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38185264) q[2];
sx q[2];
rz(-1.4421041) q[2];
sx q[2];
rz(-0.51634386) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84032816) q[1];
sx q[1];
rz(-2.5003308) q[1];
sx q[1];
rz(-2.458359) q[1];
rz(-2.1170656) q[3];
sx q[3];
rz(-1.3427991) q[3];
sx q[3];
rz(-0.25817623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-2.0097849) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(0.65761956) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(0.89362842) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052505715) q[0];
sx q[0];
rz(-0.94852175) q[0];
sx q[0];
rz(-1.1629348) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0412752) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(1.668001) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.212008) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(0.68540539) q[1];
rz(-2.3818447) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(-0.37978803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787518) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-2.0513127) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(-1.1539248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2147804) q[0];
sx q[0];
rz(-1.6048389) q[0];
sx q[0];
rz(-2.4916617) q[0];
rz(-0.55118982) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(-1.7984496) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1591783) q[1];
sx q[1];
rz(-2.3783861) q[1];
sx q[1];
rz(0.089341954) q[1];
x q[2];
rz(1.0941986) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.170084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3568048) q[0];
sx q[0];
rz(-0.49236449) q[0];
sx q[0];
rz(0.82226336) q[0];
rz(-pi) q[1];
rz(-2.3652472) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(-2.5459144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6744086) q[1];
sx q[1];
rz(-1.3896175) q[1];
sx q[1];
rz(0.49023899) q[1];
rz(-1.7438668) q[3];
sx q[3];
rz(-0.67865463) q[3];
sx q[3];
rz(-2.9242587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(0.66463566) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-0.47992596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23586789) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-0.41525526) q[0];
rz(-2.2357113) q[2];
sx q[2];
rz(-1.7582298) q[2];
sx q[2];
rz(0.71042176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(0.6154284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2009833) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(1.9711232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(-2.813415) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(0.27771523) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
