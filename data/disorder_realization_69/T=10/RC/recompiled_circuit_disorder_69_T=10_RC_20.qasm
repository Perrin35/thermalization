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
rz(5.2678582) q[0];
sx q[0];
rz(9.8922748) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.4407053) q[0];
sx q[0];
rz(0.84530172) q[0];
rz(-pi) q[1];
rz(0.99256398) q[2];
sx q[2];
rz(-0.34878584) q[2];
sx q[2];
rz(2.0860096) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28847028) q[1];
sx q[1];
rz(-0.56325699) q[1];
sx q[1];
rz(2.017615) q[1];
x q[2];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(-2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(-2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4366348) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090416106) q[0];
sx q[0];
rz(-0.91958445) q[0];
sx q[0];
rz(0.37474664) q[0];
rz(-1.0753724) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(2.5410595) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
x q[2];
rz(-2.8855521) q[3];
sx q[3];
rz(-0.55169332) q[3];
sx q[3];
rz(-1.8191169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15741631) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-3.1047344) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-0.18951167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150869) q[0];
sx q[0];
rz(-2.3622892) q[0];
sx q[0];
rz(2.7515609) q[0];
rz(-pi) q[1];
x q[1];
rz(1.465221) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(-2.2764652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7889001) q[1];
sx q[1];
rz(-0.2460203) q[1];
sx q[1];
rz(2.9109863) q[1];
rz(-pi) q[2];
rz(0.96954815) q[3];
sx q[3];
rz(-1.7161955) q[3];
sx q[3];
rz(1.5090514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1069964) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(1.754388) q[2];
rz(-0.43131367) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(-2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-1.0901573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50463146) q[0];
sx q[0];
rz(-1.9177027) q[0];
sx q[0];
rz(-2.7142801) q[0];
rz(-pi) q[1];
rz(-0.50699373) q[2];
sx q[2];
rz(-0.88627964) q[2];
sx q[2];
rz(0.13912858) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7294457) q[1];
sx q[1];
rz(-2.4033961) q[1];
sx q[1];
rz(-0.63936887) q[1];
x q[2];
rz(1.777321) q[3];
sx q[3];
rz(-1.1518475) q[3];
sx q[3];
rz(-1.9358313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(0.91317552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(0.14973101) q[0];
rz(0.99114746) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(-1.9715462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2802551) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(-2.0656385) q[0];
rz(-pi) q[1];
rz(0.89044517) q[2];
sx q[2];
rz(-1.8257942) q[2];
sx q[2];
rz(-1.6971708) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7829195) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(2.6637117) q[1];
x q[2];
rz(1.7647469) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(-2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(-2.4687016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653172) q[0];
sx q[0];
rz(-0.99940171) q[0];
sx q[0];
rz(-3.0185643) q[0];
x q[1];
rz(-0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(1.4662454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8012961) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(-1.3188386) q[1];
rz(-pi) q[2];
rz(0.17500413) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.4020845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8334478) q[0];
sx q[0];
rz(-1.937056) q[0];
sx q[0];
rz(1.863198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8066605) q[2];
sx q[2];
rz(-2.3216322) q[2];
sx q[2];
rz(1.5136528) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40560383) q[1];
sx q[1];
rz(-0.39499184) q[1];
sx q[1];
rz(0.40785457) q[1];
x q[2];
rz(2.1590809) q[3];
sx q[3];
rz(-1.6856098) q[3];
sx q[3];
rz(-1.7867775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(0.68022234) q[2];
rz(2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.5493786) q[0];
rz(-2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(0.54135281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2490059) q[0];
sx q[0];
rz(-1.5595946) q[0];
sx q[0];
rz(2.5773994) q[0];
rz(0.63033732) q[2];
sx q[2];
rz(-0.54140831) q[2];
sx q[2];
rz(0.88592096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1047157) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(-2.5887262) q[1];
x q[2];
rz(-2.7612711) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(-2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(-3.1316277) q[0];
rz(-2.126157) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.6962956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89307071) q[0];
sx q[0];
rz(-2.0914441) q[0];
sx q[0];
rz(2.394886) q[0];
x q[1];
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
rz(-3.0874789) q[1];
sx q[1];
rz(-2.2009654) q[1];
sx q[1];
rz(-1.5515045) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.290756) q[3];
sx q[3];
rz(-2.3265504) q[3];
sx q[3];
rz(1.9581865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(2.7493431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046767226) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(0.70478583) q[0];
x q[1];
rz(2.5160772) q[2];
sx q[2];
rz(-2.054347) q[2];
sx q[2];
rz(3.0160883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85234648) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(0.5457408) q[1];
x q[2];
rz(3.0030389) q[3];
sx q[3];
rz(-2.9346653) q[3];
sx q[3];
rz(3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-1.2148576) q[2];
sx q[2];
rz(-2.3050953) q[2];
sx q[2];
rz(-0.54511025) q[2];
rz(-2.5946887) q[3];
sx q[3];
rz(-0.96596598) q[3];
sx q[3];
rz(2.1706497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
