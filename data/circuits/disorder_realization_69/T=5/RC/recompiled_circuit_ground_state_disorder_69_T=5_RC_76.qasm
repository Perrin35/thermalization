OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2562113) q[0];
sx q[0];
rz(-0.71330944) q[0];
sx q[0];
rz(-1.2137086) q[0];
rz(1.775939) q[1];
sx q[1];
rz(-0.27888271) q[1];
sx q[1];
rz(-2.9409148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047223481) q[0];
sx q[0];
rz(-0.46877623) q[0];
sx q[0];
rz(-1.3865406) q[0];
rz(-1.7332478) q[2];
sx q[2];
rz(-1.7853569) q[2];
sx q[2];
rz(-1.868737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37113443) q[1];
sx q[1];
rz(-1.374303) q[1];
sx q[1];
rz(-1.3875828) q[1];
x q[2];
rz(2.8332769) q[3];
sx q[3];
rz(-1.4065969) q[3];
sx q[3];
rz(-2.7621244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0783477) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(-0.71301785) q[2];
rz(-1.8578505) q[3];
sx q[3];
rz(-0.72834891) q[3];
sx q[3];
rz(-2.242105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9707608) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(-1.4738039) q[0];
rz(0.57227349) q[1];
sx q[1];
rz(-1.0621366) q[1];
sx q[1];
rz(-1.7639814) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.361918) q[0];
sx q[0];
rz(-0.55680823) q[0];
sx q[0];
rz(2.1614055) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1123231) q[2];
sx q[2];
rz(-1.0857669) q[2];
sx q[2];
rz(0.84225076) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.936539) q[1];
sx q[1];
rz(-1.2779826) q[1];
sx q[1];
rz(1.9163301) q[1];
x q[2];
rz(2.1491941) q[3];
sx q[3];
rz(-0.36250329) q[3];
sx q[3];
rz(0.9622919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1231692) q[2];
sx q[2];
rz(-1.6320684) q[2];
sx q[2];
rz(1.8040166) q[2];
rz(-2.7959974) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(-2.2071297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551986) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(-0.58919543) q[0];
rz(-2.9335754) q[1];
sx q[1];
rz(-0.55744019) q[1];
sx q[1];
rz(-0.20588188) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34768057) q[0];
sx q[0];
rz(-0.97647515) q[0];
sx q[0];
rz(-2.0012614) q[0];
rz(-1.0127064) q[2];
sx q[2];
rz(-0.79540247) q[2];
sx q[2];
rz(0.87855065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9045453) q[1];
sx q[1];
rz(-2.3694607) q[1];
sx q[1];
rz(0.44293483) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7725189) q[3];
sx q[3];
rz(-2.7274787) q[3];
sx q[3];
rz(0.61455807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0516633) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(0.54452407) q[2];
rz(-0.45474592) q[3];
sx q[3];
rz(-2.5179722) q[3];
sx q[3];
rz(3.1375942) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398448) q[0];
sx q[0];
rz(-0.54773098) q[0];
sx q[0];
rz(2.4416583) q[0];
rz(-1.8327389) q[1];
sx q[1];
rz(-1.0700285) q[1];
sx q[1];
rz(1.7226137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15407059) q[0];
sx q[0];
rz(-0.51263499) q[0];
sx q[0];
rz(1.1521167) q[0];
rz(-0.019446418) q[2];
sx q[2];
rz(-2.0692973) q[2];
sx q[2];
rz(-2.4011103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1261038) q[1];
sx q[1];
rz(-2.2601193) q[1];
sx q[1];
rz(1.8828805) q[1];
rz(1.254564) q[3];
sx q[3];
rz(-1.4319311) q[3];
sx q[3];
rz(2.3854802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2122638) q[2];
sx q[2];
rz(-0.45390359) q[2];
sx q[2];
rz(-1.7960499) q[2];
rz(0.71873194) q[3];
sx q[3];
rz(-0.74145442) q[3];
sx q[3];
rz(-2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6253925) q[0];
sx q[0];
rz(-1.9560408) q[0];
sx q[0];
rz(0.54310435) q[0];
rz(0.25554666) q[1];
sx q[1];
rz(-0.4823904) q[1];
sx q[1];
rz(-1.7593174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10777625) q[0];
sx q[0];
rz(-1.9703426) q[0];
sx q[0];
rz(2.5679213) q[0];
x q[1];
rz(0.76680317) q[2];
sx q[2];
rz(-1.4652952) q[2];
sx q[2];
rz(-2.2720384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68866036) q[1];
sx q[1];
rz(-2.1218461) q[1];
sx q[1];
rz(2.8351985) q[1];
rz(-pi) q[2];
rz(2.7543254) q[3];
sx q[3];
rz(-2.7125689) q[3];
sx q[3];
rz(2.9128592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9690659) q[2];
sx q[2];
rz(-1.4719937) q[2];
sx q[2];
rz(-2.772061) q[2];
rz(-1.8116123) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(1.3548405) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0238277) q[0];
sx q[0];
rz(-0.33482877) q[0];
sx q[0];
rz(-3.0552926) q[0];
rz(1.6465126) q[1];
sx q[1];
rz(-1.1777271) q[1];
sx q[1];
rz(-2.7118447) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54329007) q[0];
sx q[0];
rz(-1.6409281) q[0];
sx q[0];
rz(-1.381641) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4835973) q[2];
sx q[2];
rz(-1.832973) q[2];
sx q[2];
rz(3.0397751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6657998) q[1];
sx q[1];
rz(-2.0415039) q[1];
sx q[1];
rz(2.0657971) q[1];
x q[2];
rz(-1.2901575) q[3];
sx q[3];
rz(-2.2220816) q[3];
sx q[3];
rz(0.7574581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0325844) q[2];
sx q[2];
rz(-2.0592368) q[2];
sx q[2];
rz(-1.9411795) q[2];
rz(0.16610185) q[3];
sx q[3];
rz(-2.3717272) q[3];
sx q[3];
rz(-0.86336819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85668856) q[0];
sx q[0];
rz(-0.76304522) q[0];
sx q[0];
rz(-0.82409182) q[0];
rz(1.9169982) q[1];
sx q[1];
rz(-1.9173744) q[1];
sx q[1];
rz(-2.1276988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4201502) q[0];
sx q[0];
rz(-0.80810302) q[0];
sx q[0];
rz(-1.1366913) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43724202) q[2];
sx q[2];
rz(-1.8595942) q[2];
sx q[2];
rz(2.0721958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6142004) q[1];
sx q[1];
rz(-2.1062615) q[1];
sx q[1];
rz(2.0245069) q[1];
x q[2];
rz(-0.39172642) q[3];
sx q[3];
rz(-1.0684895) q[3];
sx q[3];
rz(1.4823748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27807221) q[2];
sx q[2];
rz(-2.4188228) q[2];
sx q[2];
rz(1.2449167) q[2];
rz(0.23166367) q[3];
sx q[3];
rz(-2.1004227) q[3];
sx q[3];
rz(-0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(2.046655) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(2.9723883) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(-1.4134891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.468576) q[0];
sx q[0];
rz(-1.4961188) q[0];
sx q[0];
rz(-0.59665307) q[0];
rz(-2.1229073) q[2];
sx q[2];
rz(-1.3231965) q[2];
sx q[2];
rz(-1.1820861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1991829) q[1];
sx q[1];
rz(-1.8687667) q[1];
sx q[1];
rz(1.5375141) q[1];
x q[2];
rz(0.19988536) q[3];
sx q[3];
rz(-2.0428786) q[3];
sx q[3];
rz(2.1560644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32470545) q[2];
sx q[2];
rz(-1.4810666) q[2];
sx q[2];
rz(1.7186349) q[2];
rz(0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(0.62937361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9519575) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(0.54963175) q[0];
rz(0.59898392) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(-1.6302861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6236512) q[0];
sx q[0];
rz(-1.7474801) q[0];
sx q[0];
rz(1.5672632) q[0];
rz(1.8725996) q[2];
sx q[2];
rz(-1.0827218) q[2];
sx q[2];
rz(-2.7583964) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0907184) q[1];
sx q[1];
rz(-1.9719187) q[1];
sx q[1];
rz(-2.119333) q[1];
x q[2];
rz(-0.49137791) q[3];
sx q[3];
rz(-1.2983606) q[3];
sx q[3];
rz(2.2898554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0281684) q[2];
sx q[2];
rz(-1.7999962) q[2];
sx q[2];
rz(0.33451954) q[2];
rz(-2.4962375) q[3];
sx q[3];
rz(-2.1367475) q[3];
sx q[3];
rz(0.2373124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0765814) q[0];
sx q[0];
rz(-2.8028221) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(-2.8399641) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-1.041144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89462507) q[0];
sx q[0];
rz(-1.1389705) q[0];
sx q[0];
rz(-2.4869842) q[0];
rz(-pi) q[1];
rz(-1.663346) q[2];
sx q[2];
rz(-2.0299533) q[2];
sx q[2];
rz(-0.29806229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54339441) q[1];
sx q[1];
rz(-0.80401995) q[1];
sx q[1];
rz(-1.6135741) q[1];
rz(-pi) q[2];
rz(0.21386336) q[3];
sx q[3];
rz(-1.920061) q[3];
sx q[3];
rz(-1.6598778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4447896) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(0.94318548) q[2];
rz(-1.7275564) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(0.79677719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.4010451) q[0];
sx q[0];
rz(-2.7475806) q[0];
sx q[0];
rz(-0.077234118) q[0];
rz(-2.7229591) q[1];
sx q[1];
rz(-1.5822462) q[1];
sx q[1];
rz(0.15204631) q[1];
rz(-0.05337333) q[2];
sx q[2];
rz(-2.5659701) q[2];
sx q[2];
rz(0.26502668) q[2];
rz(-0.43155117) q[3];
sx q[3];
rz(-1.3234371) q[3];
sx q[3];
rz(0.32061843) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
