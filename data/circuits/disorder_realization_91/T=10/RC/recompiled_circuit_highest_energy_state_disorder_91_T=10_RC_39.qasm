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
rz(0.089601547) q[0];
sx q[0];
rz(-0.90553415) q[0];
sx q[0];
rz(2.9517458) q[0];
rz(2.7449961) q[1];
sx q[1];
rz(-0.34062579) q[1];
sx q[1];
rz(-0.61000282) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0002132) q[0];
sx q[0];
rz(-1.3249517) q[0];
sx q[0];
rz(-0.2443929) q[0];
rz(-pi) q[1];
rz(1.1930614) q[2];
sx q[2];
rz(-0.47171041) q[2];
sx q[2];
rz(-0.38208252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0335281) q[1];
sx q[1];
rz(-0.92108291) q[1];
sx q[1];
rz(2.482138) q[1];
x q[2];
rz(0.27746706) q[3];
sx q[3];
rz(-1.101081) q[3];
sx q[3];
rz(-1.0535002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4492599) q[2];
sx q[2];
rz(-1.2168987) q[2];
sx q[2];
rz(-2.8999691) q[2];
rz(-0.062945098) q[3];
sx q[3];
rz(-1.2017622) q[3];
sx q[3];
rz(2.3285749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45074201) q[0];
sx q[0];
rz(-1.2797322) q[0];
sx q[0];
rz(0.29020852) q[0];
rz(0.33503512) q[1];
sx q[1];
rz(-2.2033043) q[1];
sx q[1];
rz(0.32523528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9262518) q[0];
sx q[0];
rz(-1.6115973) q[0];
sx q[0];
rz(0.46567877) q[0];
x q[1];
rz(-1.2464855) q[2];
sx q[2];
rz(-0.70581573) q[2];
sx q[2];
rz(-1.2845662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57978153) q[1];
sx q[1];
rz(-2.9159286) q[1];
sx q[1];
rz(0.8439941) q[1];
rz(-pi) q[2];
rz(2.144037) q[3];
sx q[3];
rz(-1.0008423) q[3];
sx q[3];
rz(0.23373014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30514303) q[2];
sx q[2];
rz(-0.93905753) q[2];
sx q[2];
rz(1.1391501) q[2];
rz(-0.75974733) q[3];
sx q[3];
rz(-3.0436438) q[3];
sx q[3];
rz(-2.7636512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515795) q[0];
sx q[0];
rz(-2.4115998) q[0];
sx q[0];
rz(-2.3082025) q[0];
rz(-0.003412811) q[1];
sx q[1];
rz(-0.5178057) q[1];
sx q[1];
rz(1.980967) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072642751) q[0];
sx q[0];
rz(-2.8153689) q[0];
sx q[0];
rz(1.6618927) q[0];
x q[1];
rz(2.848804) q[2];
sx q[2];
rz(-0.66899111) q[2];
sx q[2];
rz(-2.3716253) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1091253) q[1];
sx q[1];
rz(-1.3360436) q[1];
sx q[1];
rz(0.0003628022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9367552) q[3];
sx q[3];
rz(-1.7410918) q[3];
sx q[3];
rz(-2.1619201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2680336) q[2];
sx q[2];
rz(-1.4159091) q[2];
sx q[2];
rz(0.1669008) q[2];
rz(1.1901101) q[3];
sx q[3];
rz(-0.24470617) q[3];
sx q[3];
rz(2.0580097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0077008) q[0];
sx q[0];
rz(-1.5950483) q[0];
sx q[0];
rz(-2.290945) q[0];
rz(-2.3386686) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(2.0096774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7106774) q[0];
sx q[0];
rz(-1.3899904) q[0];
sx q[0];
rz(2.0742832) q[0];
rz(-pi) q[1];
rz(2.1430127) q[2];
sx q[2];
rz(-2.7957186) q[2];
sx q[2];
rz(1.5550176) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4635361) q[1];
sx q[1];
rz(-0.66036036) q[1];
sx q[1];
rz(2.6956431) q[1];
x q[2];
rz(-0.72992562) q[3];
sx q[3];
rz(-1.2582964) q[3];
sx q[3];
rz(-0.70532986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5256646) q[2];
sx q[2];
rz(-1.5027081) q[2];
sx q[2];
rz(-0.28044236) q[2];
rz(2.7403455) q[3];
sx q[3];
rz(-0.29956996) q[3];
sx q[3];
rz(-1.7846599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.930437) q[0];
sx q[0];
rz(-1.7481952) q[0];
sx q[0];
rz(0.18779553) q[0];
rz(-2.4941749) q[1];
sx q[1];
rz(-2.1719666) q[1];
sx q[1];
rz(-0.59026778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2881516) q[0];
sx q[0];
rz(-0.36124215) q[0];
sx q[0];
rz(1.4807184) q[0];
rz(-1.9365385) q[2];
sx q[2];
rz(-2.3252333) q[2];
sx q[2];
rz(-2.2162645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94475585) q[1];
sx q[1];
rz(-1.6513255) q[1];
sx q[1];
rz(1.5989941) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3594317) q[3];
sx q[3];
rz(-2.2244033) q[3];
sx q[3];
rz(0.048649064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.65842998) q[2];
sx q[2];
rz(-1.1700609) q[2];
sx q[2];
rz(-0.17407334) q[2];
rz(0.46711323) q[3];
sx q[3];
rz(-2.4090448) q[3];
sx q[3];
rz(-2.9113801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9664522) q[0];
sx q[0];
rz(-1.3850965) q[0];
sx q[0];
rz(-1.0873644) q[0];
rz(-0.1611791) q[1];
sx q[1];
rz(-1.168707) q[1];
sx q[1];
rz(-2.2081614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0421222) q[0];
sx q[0];
rz(-0.5802896) q[0];
sx q[0];
rz(0.97162928) q[0];
rz(-pi) q[1];
rz(-1.1968139) q[2];
sx q[2];
rz(-1.6829964) q[2];
sx q[2];
rz(0.89857946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7644706) q[1];
sx q[1];
rz(-2.5896448) q[1];
sx q[1];
rz(0.66739018) q[1];
rz(-pi) q[2];
rz(-0.94901086) q[3];
sx q[3];
rz(-2.057926) q[3];
sx q[3];
rz(2.0355527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33438385) q[2];
sx q[2];
rz(-1.7024567) q[2];
sx q[2];
rz(1.8857694) q[2];
rz(0.0028336023) q[3];
sx q[3];
rz(-0.56003672) q[3];
sx q[3];
rz(-2.475256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6781219) q[0];
sx q[0];
rz(-3.1246298) q[0];
sx q[0];
rz(-0.28067881) q[0];
rz(2.8374788) q[1];
sx q[1];
rz(-2.3801443) q[1];
sx q[1];
rz(-1.4842518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2646177) q[0];
sx q[0];
rz(-1.5533981) q[0];
sx q[0];
rz(1.5511284) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4062469) q[2];
sx q[2];
rz(-0.79716792) q[2];
sx q[2];
rz(-2.2556502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6929984) q[1];
sx q[1];
rz(-1.0585551) q[1];
sx q[1];
rz(-2.1650141) q[1];
rz(0.13388195) q[3];
sx q[3];
rz(-1.1999793) q[3];
sx q[3];
rz(-1.6755392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7228399) q[2];
sx q[2];
rz(-2.3374228) q[2];
sx q[2];
rz(1.0214825) q[2];
rz(0.57182765) q[3];
sx q[3];
rz(-2.573206) q[3];
sx q[3];
rz(0.9182601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9224213) q[0];
sx q[0];
rz(-0.89184856) q[0];
sx q[0];
rz(-3.0944371) q[0];
rz(1.7258518) q[1];
sx q[1];
rz(-1.9784617) q[1];
sx q[1];
rz(-0.39173752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120868) q[0];
sx q[0];
rz(-1.1177309) q[0];
sx q[0];
rz(1.5231251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9525346) q[2];
sx q[2];
rz(-1.9044276) q[2];
sx q[2];
rz(2.5572973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9122767) q[1];
sx q[1];
rz(-3.1247093) q[1];
sx q[1];
rz(2.2364278) q[1];
rz(-1.3545808) q[3];
sx q[3];
rz(-1.9554227) q[3];
sx q[3];
rz(-0.011354488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91531104) q[2];
sx q[2];
rz(-0.87258029) q[2];
sx q[2];
rz(-2.937781) q[2];
rz(-0.63621825) q[3];
sx q[3];
rz(-1.3602942) q[3];
sx q[3];
rz(-3.0740331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5485452) q[0];
sx q[0];
rz(-3.1199582) q[0];
sx q[0];
rz(1.1170603) q[0];
rz(-1.2720269) q[1];
sx q[1];
rz(-2.2950324) q[1];
sx q[1];
rz(0.55714947) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8654738) q[0];
sx q[0];
rz(-2.0003002) q[0];
sx q[0];
rz(-0.19709023) q[0];
rz(-pi) q[1];
rz(-1.7366055) q[2];
sx q[2];
rz(-1.4784383) q[2];
sx q[2];
rz(-0.044008642) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7751366) q[1];
sx q[1];
rz(-1.8231099) q[1];
sx q[1];
rz(1.2229678) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45175868) q[3];
sx q[3];
rz(-2.5521899) q[3];
sx q[3];
rz(1.9077993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40994) q[2];
sx q[2];
rz(-0.2356379) q[2];
sx q[2];
rz(-1.635599) q[2];
rz(-0.19009863) q[3];
sx q[3];
rz(-2.5458769) q[3];
sx q[3];
rz(0.070040919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6079123) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(-2.876907) q[0];
rz(-0.39868042) q[1];
sx q[1];
rz(-0.52972263) q[1];
sx q[1];
rz(2.9369489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35025233) q[0];
sx q[0];
rz(-2.5584284) q[0];
sx q[0];
rz(-1.9406609) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25252931) q[2];
sx q[2];
rz(-2.1871532) q[2];
sx q[2];
rz(-1.7231307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.252534) q[1];
sx q[1];
rz(-1.5100579) q[1];
sx q[1];
rz(-0.82734682) q[1];
rz(0.47765215) q[3];
sx q[3];
rz(-1.700145) q[3];
sx q[3];
rz(2.1952598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22145049) q[2];
sx q[2];
rz(-1.3416938) q[2];
sx q[2];
rz(1.3915001) q[2];
rz(0.57735389) q[3];
sx q[3];
rz(-2.5632863) q[3];
sx q[3];
rz(0.81380832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53265424) q[0];
sx q[0];
rz(-1.2464936) q[0];
sx q[0];
rz(-2.0410224) q[0];
rz(0.34250034) q[1];
sx q[1];
rz(-0.96440146) q[1];
sx q[1];
rz(-1.0920116) q[1];
rz(-2.9663646) q[2];
sx q[2];
rz(-1.2675076) q[2];
sx q[2];
rz(-2.8855973) q[2];
rz(0.92855056) q[3];
sx q[3];
rz(-2.3026569) q[3];
sx q[3];
rz(2.3893366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
