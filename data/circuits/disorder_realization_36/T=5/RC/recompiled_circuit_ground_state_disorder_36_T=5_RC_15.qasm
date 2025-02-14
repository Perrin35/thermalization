OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(-1.9195317) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9400222) q[0];
sx q[0];
rz(-2.4594677) q[0];
sx q[0];
rz(-1.8322102) q[0];
x q[1];
rz(1.0009345) q[2];
sx q[2];
rz(-2.9409932) q[2];
sx q[2];
rz(2.4127485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.192321) q[1];
sx q[1];
rz(-1.9931798) q[1];
sx q[1];
rz(1.3793089) q[1];
rz(2.7056115) q[3];
sx q[3];
rz(-1.1262622) q[3];
sx q[3];
rz(-1.7874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1808971) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(2.1107471) q[2];
rz(-2.893462) q[3];
sx q[3];
rz(-1.3458359) q[3];
sx q[3];
rz(0.26721755) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2613075) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(2.0174568) q[0];
rz(-2.0696438) q[1];
sx q[1];
rz(-1.5771461) q[1];
sx q[1];
rz(-2.7148278) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65703553) q[0];
sx q[0];
rz(-0.74321514) q[0];
sx q[0];
rz(1.5748181) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8324805) q[2];
sx q[2];
rz(-1.2361457) q[2];
sx q[2];
rz(-1.4317012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76600915) q[1];
sx q[1];
rz(-1.6787306) q[1];
sx q[1];
rz(2.4426779) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20189607) q[3];
sx q[3];
rz(-1.1888973) q[3];
sx q[3];
rz(2.524162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65853226) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(0.10022441) q[2];
rz(-2.0641067) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(1.235435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60233068) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(1.7072898) q[0];
rz(-2.3151248) q[1];
sx q[1];
rz(-1.6441328) q[1];
sx q[1];
rz(2.7109587) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.193482) q[0];
sx q[0];
rz(-1.3010345) q[0];
sx q[0];
rz(-1.7408235) q[0];
rz(-pi) q[1];
rz(3.011322) q[2];
sx q[2];
rz(-1.8415057) q[2];
sx q[2];
rz(0.66616466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34738008) q[1];
sx q[1];
rz(-0.93232226) q[1];
sx q[1];
rz(-2.641851) q[1];
rz(-pi) q[2];
rz(1.0266193) q[3];
sx q[3];
rz(-1.9154356) q[3];
sx q[3];
rz(-2.1066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2677801) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(2.7576647) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(2.4887776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225994) q[0];
sx q[0];
rz(-0.79177952) q[0];
sx q[0];
rz(1.0571887) q[0];
rz(2.3253697) q[1];
sx q[1];
rz(-0.99333251) q[1];
sx q[1];
rz(2.9073471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.809991) q[0];
sx q[0];
rz(-1.1561285) q[0];
sx q[0];
rz(0.11964397) q[0];
x q[1];
rz(-2.0969056) q[2];
sx q[2];
rz(-1.1981694) q[2];
sx q[2];
rz(-1.5043196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1409451) q[1];
sx q[1];
rz(-1.06303) q[1];
sx q[1];
rz(0.80054466) q[1];
x q[2];
rz(0.81984249) q[3];
sx q[3];
rz(-2.4944135) q[3];
sx q[3];
rz(2.1539713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8546042) q[2];
sx q[2];
rz(-1.8082666) q[2];
sx q[2];
rz(2.2387779) q[2];
rz(-1.3826987) q[3];
sx q[3];
rz(-0.32302502) q[3];
sx q[3];
rz(-0.84669101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803902) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(-0.70988208) q[0];
rz(-0.4711802) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(-3.0772298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854537) q[0];
sx q[0];
rz(-2.6366933) q[0];
sx q[0];
rz(-1.8527777) q[0];
rz(0.2336524) q[2];
sx q[2];
rz(-2.0548247) q[2];
sx q[2];
rz(0.40586995) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2272294) q[1];
sx q[1];
rz(-1.0681947) q[1];
sx q[1];
rz(-1.0881626) q[1];
rz(1.5803171) q[3];
sx q[3];
rz(-1.4839982) q[3];
sx q[3];
rz(-2.1824126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(-0.95332471) q[2];
rz(-0.55771762) q[3];
sx q[3];
rz(-1.4671114) q[3];
sx q[3];
rz(-0.52887708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271269) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(-0.99660981) q[0];
rz(-0.059545513) q[1];
sx q[1];
rz(-2.153502) q[1];
sx q[1];
rz(-1.7893808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4137581) q[0];
sx q[0];
rz(-2.067885) q[0];
sx q[0];
rz(2.1700806) q[0];
x q[1];
rz(0.33783317) q[2];
sx q[2];
rz(-2.3408836) q[2];
sx q[2];
rz(1.5986625) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3660197) q[1];
sx q[1];
rz(-1.8692353) q[1];
sx q[1];
rz(1.8587023) q[1];
x q[2];
rz(-2.2289946) q[3];
sx q[3];
rz(-1.8479947) q[3];
sx q[3];
rz(-2.4656221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(2.582029) q[2];
rz(-2.3061421) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(1.8387509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68592042) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(0.36668229) q[0];
rz(-2.2600251) q[1];
sx q[1];
rz(-0.63789788) q[1];
sx q[1];
rz(-1.4311904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0019497) q[0];
sx q[0];
rz(-0.77660034) q[0];
sx q[0];
rz(2.5703672) q[0];
rz(-pi) q[1];
rz(-2.7273799) q[2];
sx q[2];
rz(-0.72986947) q[2];
sx q[2];
rz(3.102906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8579432) q[1];
sx q[1];
rz(-1.6873296) q[1];
sx q[1];
rz(-1.5416508) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8058947) q[3];
sx q[3];
rz(-2.2155016) q[3];
sx q[3];
rz(-0.39261445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35393474) q[2];
sx q[2];
rz(-0.5160318) q[2];
sx q[2];
rz(-2.3616974) q[2];
rz(-0.53906131) q[3];
sx q[3];
rz(-1.6293679) q[3];
sx q[3];
rz(2.5299634) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014538177) q[0];
sx q[0];
rz(-2.4712565) q[0];
sx q[0];
rz(-2.3624453) q[0];
rz(0.51891333) q[1];
sx q[1];
rz(-1.2823558) q[1];
sx q[1];
rz(-0.40294495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084672734) q[0];
sx q[0];
rz(-2.2209327) q[0];
sx q[0];
rz(2.5385602) q[0];
x q[1];
rz(2.05415) q[2];
sx q[2];
rz(-1.8515808) q[2];
sx q[2];
rz(-1.3705685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.886595) q[1];
sx q[1];
rz(-2.1817921) q[1];
sx q[1];
rz(-2.0605036) q[1];
rz(-1.9698079) q[3];
sx q[3];
rz(-0.61099377) q[3];
sx q[3];
rz(-2.4965167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0275823) q[2];
sx q[2];
rz(-1.1677914) q[2];
sx q[2];
rz(1.1172509) q[2];
rz(0.69636017) q[3];
sx q[3];
rz(-2.0317234) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47525147) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(2.7269205) q[0];
rz(1.0617725) q[1];
sx q[1];
rz(-2.2551408) q[1];
sx q[1];
rz(-0.83962238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3436326) q[0];
sx q[0];
rz(-1.1251483) q[0];
sx q[0];
rz(0.4526505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62445415) q[2];
sx q[2];
rz(-0.82459282) q[2];
sx q[2];
rz(-2.7471682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.190358) q[1];
sx q[1];
rz(-2.2273185) q[1];
sx q[1];
rz(1.3912398) q[1];
x q[2];
rz(1.3247588) q[3];
sx q[3];
rz(-1.2241505) q[3];
sx q[3];
rz(-2.3510166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7825369) q[2];
sx q[2];
rz(-1.7724719) q[2];
sx q[2];
rz(-0.071652023) q[2];
rz(0.045529384) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(2.6825452) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59938207) q[0];
sx q[0];
rz(-2.5686503) q[0];
sx q[0];
rz(-1.0546767) q[0];
rz(0.032546267) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(0.5221101) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73682937) q[0];
sx q[0];
rz(-2.1407581) q[0];
sx q[0];
rz(-1.9362157) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5100461) q[2];
sx q[2];
rz(-0.8435404) q[2];
sx q[2];
rz(-0.90383672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73192683) q[1];
sx q[1];
rz(-2.5931011) q[1];
sx q[1];
rz(2.3540703) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22618146) q[3];
sx q[3];
rz(-1.0345248) q[3];
sx q[3];
rz(0.29335833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1471499) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(-1.408255) q[2];
rz(-1.555892) q[3];
sx q[3];
rz(-2.9097911) q[3];
sx q[3];
rz(-2.2304992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-2.945902) q[2];
sx q[2];
rz(-1.2441845) q[2];
sx q[2];
rz(1.2938538) q[2];
rz(0.23789858) q[3];
sx q[3];
rz(-2.4645505) q[3];
sx q[3];
rz(-0.6445618) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
