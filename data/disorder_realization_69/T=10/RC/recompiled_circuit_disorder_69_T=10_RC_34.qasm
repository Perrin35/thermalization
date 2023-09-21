OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(-0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(-0.68038124) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8898534) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(-2.2962909) q[0];
rz(-pi) q[1];
rz(-0.99256398) q[2];
sx q[2];
rz(-0.34878584) q[2];
sx q[2];
rz(-2.0860096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.8036098) q[1];
sx q[1];
rz(-1.0531055) q[1];
rz(-pi) q[2];
rz(2.666113) q[3];
sx q[3];
rz(-0.37215044) q[3];
sx q[3];
rz(-1.1572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(-0.087466784) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-2.6175245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.24633) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(2.2569879) q[0];
x q[1];
rz(2.2696813) q[2];
sx q[2];
rz(-1.2362091) q[2];
sx q[2];
rz(-1.3427693) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9973035) q[1];
sx q[1];
rz(-1.803949) q[1];
sx q[1];
rz(0.38645978) q[1];
rz(-0.25604053) q[3];
sx q[3];
rz(-2.5898993) q[3];
sx q[3];
rz(1.3224758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15741631) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7150869) q[0];
sx q[0];
rz(-2.3622892) q[0];
sx q[0];
rz(-2.7515609) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37522845) q[2];
sx q[2];
rz(-1.4725176) q[2];
sx q[2];
rz(0.74429846) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3526926) q[1];
sx q[1];
rz(-2.8955724) q[1];
sx q[1];
rz(0.23060631) q[1];
x q[2];
rz(-1.3174921) q[3];
sx q[3];
rz(-0.6164624) q[3];
sx q[3];
rz(-2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65961924) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(2.0514354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6369612) q[0];
sx q[0];
rz(-1.22389) q[0];
sx q[0];
rz(-0.42731254) q[0];
rz(-pi) q[1];
rz(0.81972576) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(2.0476598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66149368) q[1];
sx q[1];
rz(-1.9839994) q[1];
sx q[1];
rz(2.5109629) q[1];
rz(-1.3642717) q[3];
sx q[3];
rz(-1.1518475) q[3];
sx q[3];
rz(1.2057613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1426992) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(0.91317552) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(0.14973101) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67384185) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(-0.90676102) q[0];
rz(-0.32355002) q[2];
sx q[2];
rz(-2.2252482) q[2];
sx q[2];
rz(0.3277342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7829195) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(2.6637117) q[1];
rz(-1.7647469) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(-0.84996163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85149375) q[0];
sx q[0];
rz(-0.58304542) q[0];
sx q[0];
rz(1.3821938) q[0];
x q[1];
rz(-0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(-1.6753472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7676047) q[1];
sx q[1];
rz(-1.3629706) q[1];
sx q[1];
rz(-2.5307104) q[1];
x q[2];
rz(-2.9665885) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(-1.1328896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3488397) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7715493) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(-2.7605961) q[0];
x q[1];
rz(0.24537556) q[2];
sx q[2];
rz(-0.78005314) q[2];
sx q[2];
rz(-1.1749554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84343108) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(1.4069188) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7756895) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(-2.7555335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(0.24958615) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(2.6002398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2490059) q[0];
sx q[0];
rz(-1.5595946) q[0];
sx q[0];
rz(-2.5773994) q[0];
rz(-pi) q[1];
rz(-2.5112553) q[2];
sx q[2];
rz(-0.54140831) q[2];
sx q[2];
rz(0.88592096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1047157) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(-0.55286644) q[1];
rz(-0.0046644966) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(-0.62957803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(3.1316277) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.4452971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709258) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(-0.70113457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5601013) q[2];
sx q[2];
rz(-0.98052374) q[2];
sx q[2];
rz(1.7828538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0547486) q[1];
sx q[1];
rz(-0.63042414) q[1];
sx q[1];
rz(3.1151506) q[1];
rz(2.290756) q[3];
sx q[3];
rz(-2.3265504) q[3];
sx q[3];
rz(1.1834061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-2.5342069) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(2.0954258) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(2.7493431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(-2.4368068) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4105083) q[2];
sx q[2];
rz(-0.77027551) q[2];
sx q[2];
rz(-2.0172271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85234648) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(-0.5457408) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0030389) q[3];
sx q[3];
rz(-2.9346653) q[3];
sx q[3];
rz(-3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4338715) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-2.773182) q[2];
sx q[2];
rz(-0.80130063) q[2];
sx q[2];
rz(-0.038566312) q[2];
rz(-0.546904) q[3];
sx q[3];
rz(-2.1756267) q[3];
sx q[3];
rz(-0.970943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];