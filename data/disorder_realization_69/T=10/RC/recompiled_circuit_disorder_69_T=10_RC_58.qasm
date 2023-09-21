OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(-0.68038124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17380938) q[0];
sx q[0];
rz(-0.73497226) q[0];
sx q[0];
rz(-1.7654788) q[0];
rz(-1.8664076) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(-2.0762028) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28847028) q[1];
sx q[1];
rz(-2.5783357) q[1];
sx q[1];
rz(1.1239777) q[1];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(-0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9900069) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(-2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4366348) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.24633) q[0];
sx q[0];
rz(-1.8662013) q[0];
sx q[0];
rz(-0.88460478) q[0];
rz(1.0753724) q[2];
sx q[2];
rz(-0.76250695) q[2];
sx q[2];
rz(2.5410595) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9429051) q[1];
sx q[1];
rz(-0.44829475) q[1];
sx q[1];
rz(-0.56221902) q[1];
rz(0.53701138) q[3];
sx q[3];
rz(-1.4376663) q[3];
sx q[3];
rz(2.6739124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-0.036858233) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79725093) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(2.952081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0441168) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(-1.2114899) q[0];
x q[1];
rz(1.465221) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(-2.2764652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7889001) q[1];
sx q[1];
rz(-0.2460203) q[1];
sx q[1];
rz(2.9109863) q[1];
x q[2];
rz(-1.3174921) q[3];
sx q[3];
rz(-0.6164624) q[3];
sx q[3];
rz(-2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.754388) q[2];
rz(2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-2.0514354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2290303) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(-1.9489261) q[0];
x q[1];
rz(-1.0341187) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(0.85988753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66149368) q[1];
sx q[1];
rz(-1.9839994) q[1];
sx q[1];
rz(2.5109629) q[1];
rz(2.7146043) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-0.45005902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(-0.85582716) q[2];
rz(1.1936197) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(-1.9715462) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8613375) q[0];
sx q[0];
rz(-0.72685188) q[0];
sx q[0];
rz(-2.0656385) q[0];
x q[1];
rz(-1.1779551) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(2.9658085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7331446) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(-1.1327729) q[1];
x q[2];
rz(-2.6036161) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18734922) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(1.3504008) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2900989) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(-1.3821938) q[0];
x q[1];
rz(2.2851373) q[2];
sx q[2];
rz(-1.5894784) q[2];
sx q[2];
rz(0.088353889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91025483) q[1];
sx q[1];
rz(-2.5006223) q[1];
sx q[1];
rz(0.35229589) q[1];
rz(-0.17500413) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(-2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(-1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(1.4020845) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37004334) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(-0.38099654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24537556) q[2];
sx q[2];
rz(-2.3615395) q[2];
sx q[2];
rz(-1.9666372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3558823) q[1];
sx q[1];
rz(-1.4175698) q[1];
sx q[1];
rz(-0.36550891) q[1];
rz(3.0038463) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(-3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4222251) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(0.68022234) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(-1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.592214) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(2.6002398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-0.56429243) q[0];
sx q[0];
rz(-0.020945992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9114242) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(-2.9609749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1047157) q[1];
sx q[1];
rz(-2.2762205) q[1];
sx q[1];
rz(-0.55286644) q[1];
x q[2];
rz(-1.5689315) q[3];
sx q[3];
rz(-1.1904753) q[3];
sx q[3];
rz(2.5170381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580496) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(-3.1316277) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709258) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(2.4404581) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89500918) q[2];
sx q[2];
rz(-2.0446606) q[2];
sx q[2];
rz(0.5627788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5053133) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(0.63025766) q[1];
x q[2];
rz(0.61049283) q[3];
sx q[3];
rz(-2.1497512) q[3];
sx q[3];
rz(1.0510774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69892591) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
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
rz(-2.0454355) q[0];
sx q[0];
rz(-0.7924315) q[0];
sx q[0];
rz(-0.57604726) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(3.0160883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0543921) q[1];
sx q[1];
rz(-1.1322081) q[1];
sx q[1];
rz(-2.2589576) q[1];
rz(-1.5418105) q[3];
sx q[3];
rz(-1.7757123) q[3];
sx q[3];
rz(2.8857185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-0.60047853) q[2];
rz(2.0742119) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-1.2148576) q[2];
sx q[2];
rz(-2.3050953) q[2];
sx q[2];
rz(-0.54511025) q[2];
rz(-2.2511803) q[3];
sx q[3];
rz(-1.1287516) q[3];
sx q[3];
rz(0.26656084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
