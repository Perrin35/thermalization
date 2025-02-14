OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(3.6977036) q[0];
sx q[0];
rz(7.2365427) q[0];
rz(-3.1218627) q[1];
sx q[1];
rz(-1.035773) q[1];
sx q[1];
rz(2.1929725) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.106038) q[0];
sx q[0];
rz(-2.4427938) q[0];
sx q[0];
rz(0.6684371) q[0];
rz(-pi) q[1];
rz(0.011587338) q[2];
sx q[2];
rz(-2.9040948) q[2];
sx q[2];
rz(-1.604014) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.030101209) q[1];
sx q[1];
rz(-1.9989357) q[1];
sx q[1];
rz(2.4011022) q[1];
x q[2];
rz(-2.7953732) q[3];
sx q[3];
rz(-1.0448714) q[3];
sx q[3];
rz(-0.093737515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76401508) q[2];
sx q[2];
rz(-0.074187584) q[2];
sx q[2];
rz(0.78262502) q[2];
rz(-3.022656) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(-2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9453732) q[0];
sx q[0];
rz(-1.1213028) q[0];
sx q[0];
rz(-2.8837606) q[0];
rz(0.091015426) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(-1.4965422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5071867) q[0];
sx q[0];
rz(-1.0189462) q[0];
sx q[0];
rz(0.080291434) q[0];
x q[1];
rz(0.74380959) q[2];
sx q[2];
rz(-2.653476) q[2];
sx q[2];
rz(2.0632191) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2014034) q[1];
sx q[1];
rz(-1.1588105) q[1];
sx q[1];
rz(-0.03279107) q[1];
x q[2];
rz(-2.6378651) q[3];
sx q[3];
rz(-2.4085975) q[3];
sx q[3];
rz(-0.021857787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1220793) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(-2.7369734) q[2];
rz(0.98958611) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(0.62304455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0290381) q[0];
sx q[0];
rz(-0.14415388) q[0];
sx q[0];
rz(-0.83830225) q[0];
rz(-0.84838947) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(-2.1489876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557195) q[0];
sx q[0];
rz(-0.8536754) q[0];
sx q[0];
rz(1.756987) q[0];
rz(-pi) q[1];
rz(0.16319947) q[2];
sx q[2];
rz(-0.33849785) q[2];
sx q[2];
rz(-1.6632572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6684718) q[1];
sx q[1];
rz(-1.3137914) q[1];
sx q[1];
rz(0.42196749) q[1];
rz(-pi) q[2];
rz(2.7843568) q[3];
sx q[3];
rz(-0.30870507) q[3];
sx q[3];
rz(-1.5964519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2720211) q[2];
sx q[2];
rz(-1.973899) q[2];
sx q[2];
rz(-1.1248379) q[2];
rz(0.62075067) q[3];
sx q[3];
rz(-2.1706332) q[3];
sx q[3];
rz(-3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.89897412) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(-2.9677891) q[0];
rz(0.68570343) q[1];
sx q[1];
rz(-1.4861636) q[1];
sx q[1];
rz(-2.3033843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7156775) q[0];
sx q[0];
rz(-2.0260149) q[0];
sx q[0];
rz(-2.0865738) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5694705) q[2];
sx q[2];
rz(-2.0184709) q[2];
sx q[2];
rz(1.8378225) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77828171) q[1];
sx q[1];
rz(-1.4290591) q[1];
sx q[1];
rz(2.4635386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.533314) q[3];
sx q[3];
rz(-1.1082543) q[3];
sx q[3];
rz(1.9328062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59914261) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(0.35169265) q[2];
rz(1.1951949) q[3];
sx q[3];
rz(-1.9545133) q[3];
sx q[3];
rz(-3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675932) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(2.7886673) q[0];
rz(0.30883166) q[1];
sx q[1];
rz(-1.0363204) q[1];
sx q[1];
rz(-0.028506361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3660748) q[0];
sx q[0];
rz(-0.99364508) q[0];
sx q[0];
rz(2.9329027) q[0];
x q[1];
rz(0.38178954) q[2];
sx q[2];
rz(-1.93474) q[2];
sx q[2];
rz(-0.84944968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5529768) q[1];
sx q[1];
rz(-1.7874831) q[1];
sx q[1];
rz(-0.79766794) q[1];
rz(2.0057553) q[3];
sx q[3];
rz(-1.3490145) q[3];
sx q[3];
rz(1.7580166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21294022) q[2];
sx q[2];
rz(-0.063455909) q[2];
sx q[2];
rz(2.1707936) q[2];
rz(-2.094723) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(0.48986062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8397119) q[0];
sx q[0];
rz(-1.7827001) q[0];
sx q[0];
rz(0.090959892) q[0];
rz(-1.2184527) q[1];
sx q[1];
rz(-2.8107042) q[1];
sx q[1];
rz(2.5023696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1341742) q[0];
sx q[0];
rz(-1.245541) q[0];
sx q[0];
rz(1.7132617) q[0];
rz(-0.12219723) q[2];
sx q[2];
rz(-1.7079884) q[2];
sx q[2];
rz(0.73681632) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6248577) q[1];
sx q[1];
rz(-1.0565041) q[1];
sx q[1];
rz(-0.74511294) q[1];
x q[2];
rz(1.2442144) q[3];
sx q[3];
rz(-1.5726798) q[3];
sx q[3];
rz(-2.7060946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0794534) q[2];
sx q[2];
rz(-0.97525758) q[2];
sx q[2];
rz(-0.47387588) q[2];
rz(-1.3300995) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(2.7661095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57629958) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(-2.9190049) q[0];
rz(-1.0635771) q[1];
sx q[1];
rz(-1.56366) q[1];
sx q[1];
rz(1.9482013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28119606) q[0];
sx q[0];
rz(-1.4544857) q[0];
sx q[0];
rz(-0.75006811) q[0];
rz(-pi) q[1];
rz(1.7980099) q[2];
sx q[2];
rz(-1.6889204) q[2];
sx q[2];
rz(2.3051777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2031415) q[1];
sx q[1];
rz(-1.4355772) q[1];
sx q[1];
rz(-0.44288992) q[1];
rz(-pi) q[2];
rz(1.0289331) q[3];
sx q[3];
rz(-2.8325084) q[3];
sx q[3];
rz(2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0003164) q[2];
sx q[2];
rz(-0.87339425) q[2];
sx q[2];
rz(2.5013962) q[2];
rz(0.021942465) q[3];
sx q[3];
rz(-1.7317737) q[3];
sx q[3];
rz(2.791361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6496277) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(-0.87977663) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(0.89281503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3443106) q[0];
sx q[0];
rz(-1.6160674) q[0];
sx q[0];
rz(0.08749732) q[0];
rz(-pi) q[1];
rz(2.3273349) q[2];
sx q[2];
rz(-0.99254464) q[2];
sx q[2];
rz(1.1922497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2485321) q[1];
sx q[1];
rz(-1.2888146) q[1];
sx q[1];
rz(-0.79459135) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51629169) q[3];
sx q[3];
rz(-1.1974575) q[3];
sx q[3];
rz(-0.0084127154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8799379) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(-0.086816303) q[2];
rz(2.1964729) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(-1.1300348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9574808) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(-3.0775253) q[0];
rz(1.13569) q[1];
sx q[1];
rz(-1.4402729) q[1];
sx q[1];
rz(-0.51220977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24146809) q[0];
sx q[0];
rz(-1.9808123) q[0];
sx q[0];
rz(1.2154237) q[0];
rz(-pi) q[1];
rz(1.1697328) q[2];
sx q[2];
rz(-2.3980015) q[2];
sx q[2];
rz(-2.0498073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2641702) q[1];
sx q[1];
rz(-1.3119196) q[1];
sx q[1];
rz(-1.6593462) q[1];
x q[2];
rz(3.104761) q[3];
sx q[3];
rz(-2.7324711) q[3];
sx q[3];
rz(-2.6927519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8168489) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(-1.0653488) q[2];
rz(1.6522853) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(0.092441946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5866933) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(-2.8107585) q[0];
rz(1.9031485) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(-0.27473658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58178066) q[0];
sx q[0];
rz(-1.1413478) q[0];
sx q[0];
rz(2.7665124) q[0];
rz(2.8682235) q[2];
sx q[2];
rz(-1.6022575) q[2];
sx q[2];
rz(0.69845573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52636056) q[1];
sx q[1];
rz(-0.52605275) q[1];
sx q[1];
rz(1.2895209) q[1];
x q[2];
rz(0.98311456) q[3];
sx q[3];
rz(-1.6310167) q[3];
sx q[3];
rz(-1.3863939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9773418) q[2];
sx q[2];
rz(-0.75426102) q[2];
sx q[2];
rz(-1.3573307) q[2];
rz(0.50061289) q[3];
sx q[3];
rz(-1.5631792) q[3];
sx q[3];
rz(-1.64465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.5951344) q[0];
sx q[0];
rz(-1.559579) q[0];
sx q[0];
rz(2.5478242) q[0];
rz(0.068269923) q[1];
sx q[1];
rz(-1.9207813) q[1];
sx q[1];
rz(-3.1178738) q[1];
rz(-1.4056397) q[2];
sx q[2];
rz(-1.9122304) q[2];
sx q[2];
rz(-0.0093218439) q[2];
rz(-1.7583634) q[3];
sx q[3];
rz(-1.3959342) q[3];
sx q[3];
rz(0.29831553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
