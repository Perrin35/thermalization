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
rz(2.6740958) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677833) q[0];
sx q[0];
rz(-2.4066204) q[0];
sx q[0];
rz(1.3761139) q[0];
x q[1];
rz(-1.8664076) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(1.0653898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28847028) q[1];
sx q[1];
rz(-0.56325699) q[1];
sx q[1];
rz(1.1239777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47547961) q[3];
sx q[3];
rz(-2.7694422) q[3];
sx q[3];
rz(-1.9843769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(-3.0541259) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(2.1402284) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090416106) q[0];
sx q[0];
rz(-2.2220082) q[0];
sx q[0];
rz(-0.37474664) q[0];
x q[1];
rz(2.0662202) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(2.5410595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9973035) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(2.7551329) q[1];
x q[2];
rz(1.4161795) q[3];
sx q[3];
rz(-2.1025476) q[3];
sx q[3];
rz(-1.0242517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(0.18951167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127842) q[0];
sx q[0];
rz(-1.3002987) q[0];
sx q[0];
rz(-0.74032797) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37522845) q[2];
sx q[2];
rz(-1.4725176) q[2];
sx q[2];
rz(0.74429846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1358007) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(-0.23975753) q[1];
rz(-pi) q[2];
rz(0.17574163) q[3];
sx q[3];
rz(-2.1648241) q[3];
sx q[3];
rz(-0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1069964) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(-1.0901573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91256234) q[0];
sx q[0];
rz(-1.9711442) q[0];
sx q[0];
rz(1.1926665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3218669) q[2];
sx q[2];
rz(-1.1851386) q[2];
sx q[2];
rz(-2.0476598) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5174597) q[1];
sx q[1];
rz(-1.0003261) q[1];
sx q[1];
rz(1.073451) q[1];
rz(-0.4269883) q[3];
sx q[3];
rz(-1.7592351) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1426992) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.9715462) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65483246) q[0];
sx q[0];
rz(-0.94615422) q[0];
sx q[0];
rz(0.39958582) q[0];
rz(-pi) q[1];
rz(2.8180426) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(2.8138585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.425881) q[1];
sx q[1];
rz(-0.60957805) q[1];
sx q[1];
rz(0.73519911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53797651) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(-2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.871792) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(-3.0043547) q[2];
rz(1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(-2.4687016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580299) q[0];
sx q[0];
rz(-1.467388) q[0];
sx q[0];
rz(-0.99594492) q[0];
x q[1];
rz(-0.85645533) q[2];
sx q[2];
rz(-1.5521142) q[2];
sx q[2];
rz(3.0532388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3402965) q[1];
sx q[1];
rz(-2.166689) q[1];
sx q[1];
rz(-1.8227541) q[1];
x q[2];
rz(-0.17500413) q[3];
sx q[3];
rz(-2.724218) q[3];
sx q[3];
rz(-1.1328896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(0.43506452) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(1.1116213) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.4020845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60915011) q[0];
sx q[0];
rz(-0.46450588) q[0];
sx q[0];
rz(-2.4971278) q[0];
x q[1];
rz(-2.8962171) q[2];
sx q[2];
rz(-2.3615395) q[2];
sx q[2];
rz(-1.9666372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2981616) q[1];
sx q[1];
rz(-1.9318252) q[1];
sx q[1];
rz(1.7346738) q[1];
x q[2];
rz(1.3659031) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(-2.7555335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4222251) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(0.68022234) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.359882) q[3];
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
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66051018) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(3.1206467) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9114242) q[2];
sx q[2];
rz(-1.1414141) q[2];
sx q[2];
rz(-2.9609749) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9882422) q[1];
sx q[1];
rz(-1.9820947) q[1];
sx q[1];
rz(2.3565355) q[1];
x q[2];
rz(2.7612711) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(2.3642335) q[2];
rz(-0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.6962956) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9706668) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(-2.4404581) q[0];
rz(-pi) q[1];
rz(2.5601013) q[2];
sx q[2];
rz(-2.1610689) q[2];
sx q[2];
rz(-1.7828538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.054113764) q[1];
sx q[1];
rz(-0.94062727) q[1];
sx q[1];
rz(-1.5515045) q[1];
rz(-pi) q[2];
x q[2];
rz(2.290756) q[3];
sx q[3];
rz(-0.81504226) q[3];
sx q[3];
rz(1.9581865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(-2.5342069) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059175) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(-2.5277396) q[0];
rz(-2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0961571) q[0];
sx q[0];
rz(-0.7924315) q[0];
sx q[0];
rz(2.5655454) q[0];
rz(-pi) q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(-0.12550437) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0543921) q[1];
sx q[1];
rz(-1.1322081) q[1];
sx q[1];
rz(-0.88263504) q[1];
x q[2];
rz(-3.0030389) q[3];
sx q[3];
rz(-0.20692736) q[3];
sx q[3];
rz(3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(-1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(-1.2148576) q[2];
sx q[2];
rz(-2.3050953) q[2];
sx q[2];
rz(-0.54511025) q[2];
rz(0.546904) q[3];
sx q[3];
rz(-0.96596598) q[3];
sx q[3];
rz(2.1706497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];