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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(-0.5156762) q[0];
rz(-2.9134143) q[1];
sx q[1];
rz(-2.394634) q[1];
sx q[1];
rz(2.7153314) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.212972) q[0];
sx q[0];
rz(-2.3317695) q[0];
sx q[0];
rz(2.9024603) q[0];
x q[1];
rz(-0.49246712) q[2];
sx q[2];
rz(-2.5829007) q[2];
sx q[2];
rz(-2.69489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8350153) q[1];
sx q[1];
rz(-1.5131498) q[1];
sx q[1];
rz(2.539544) q[1];
x q[2];
rz(1.0996885) q[3];
sx q[3];
rz(-2.503241) q[3];
sx q[3];
rz(-1.1717351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8969741) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(-2.8986325) q[2];
rz(2.7729559) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(1.9093556) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826137) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(2.7742703) q[0];
rz(-2.3452554) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(-0.64250362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2655661) q[0];
sx q[0];
rz(-0.6548223) q[0];
sx q[0];
rz(-0.73237082) q[0];
x q[1];
rz(2.5902469) q[2];
sx q[2];
rz(-1.2228106) q[2];
sx q[2];
rz(-0.47737338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31256235) q[1];
sx q[1];
rz(-2.1305269) q[1];
sx q[1];
rz(1.0360495) q[1];
rz(1.0890555) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(0.58603906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.7712234) q[2];
rz(-2.9141407) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74533904) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(1.1981717) q[0];
rz(-1.0802957) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(-3.0467196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085233) q[0];
sx q[0];
rz(-1.0329536) q[0];
sx q[0];
rz(-2.8406124) q[0];
x q[1];
rz(-2.2302365) q[2];
sx q[2];
rz(-1.240864) q[2];
sx q[2];
rz(-2.099769) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66902924) q[1];
sx q[1];
rz(-1.5298784) q[1];
sx q[1];
rz(1.0340235) q[1];
rz(-pi) q[2];
rz(1.1867939) q[3];
sx q[3];
rz(-1.3130762) q[3];
sx q[3];
rz(0.4674165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0723116) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.1627722) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56226319) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(-2.80559) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(0.10428183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606303) q[0];
sx q[0];
rz(-1.917836) q[0];
sx q[0];
rz(2.7913559) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5566795) q[2];
sx q[2];
rz(-1.4143362) q[2];
sx q[2];
rz(-0.9597646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6809685) q[1];
sx q[1];
rz(-1.9379341) q[1];
sx q[1];
rz(-2.2904094) q[1];
rz(-pi) q[2];
rz(1.453425) q[3];
sx q[3];
rz(-1.1201829) q[3];
sx q[3];
rz(0.71244682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61863724) q[2];
sx q[2];
rz(-2.0433661) q[2];
sx q[2];
rz(1.9970419) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(-2.0539637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2523786) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(2.8401781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2269991) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(-0.27497681) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3630939) q[2];
sx q[2];
rz(-1.1900717) q[2];
sx q[2];
rz(-2.8247339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.099286) q[1];
sx q[1];
rz(-1.9507244) q[1];
sx q[1];
rz(0.2289339) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2603792) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(2.415433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(0.0069228355) q[2];
rz(2.210468) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865006) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(-1.6078) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-0.30219561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93389171) q[0];
sx q[0];
rz(-0.91312983) q[0];
sx q[0];
rz(-3.012863) q[0];
rz(-pi) q[1];
rz(2.7415761) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(-3.0885901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8540878) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(-1.8604398) q[1];
rz(-pi) q[2];
rz(-1.2421397) q[3];
sx q[3];
rz(-0.16020003) q[3];
sx q[3];
rz(-0.70806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2493784) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(1.9269491) q[3];
sx q[3];
rz(-2.6483783) q[3];
sx q[3];
rz(-2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183384) q[0];
sx q[0];
rz(-0.76699081) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(2.8077937) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(2.0557859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54667771) q[0];
sx q[0];
rz(-1.6502066) q[0];
sx q[0];
rz(-2.9360442) q[0];
rz(-pi) q[1];
rz(-2.9606107) q[2];
sx q[2];
rz(-1.9662734) q[2];
sx q[2];
rz(0.68812319) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27583308) q[1];
sx q[1];
rz(-1.0414904) q[1];
sx q[1];
rz(-2.5342026) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2074039) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(0.31588376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51599017) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(1.9773352) q[2];
rz(2.8734251) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(-0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.5371573) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.1489768) q[0];
rz(2.8474498) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(-2.099096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4228758) q[0];
sx q[0];
rz(-1.9008753) q[0];
sx q[0];
rz(0.80083682) q[0];
rz(-pi) q[1];
rz(-2.3518219) q[2];
sx q[2];
rz(-1.7347448) q[2];
sx q[2];
rz(0.81610926) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6175872) q[1];
sx q[1];
rz(-1.397314) q[1];
sx q[1];
rz(2.7702727) q[1];
rz(-pi) q[2];
rz(-0.028826272) q[3];
sx q[3];
rz(-2.2346063) q[3];
sx q[3];
rz(0.038906038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59757549) q[2];
sx q[2];
rz(-2.47051) q[2];
sx q[2];
rz(-1.3647122) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.515601) q[0];
sx q[0];
rz(-2.7643272) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(-3.1248202) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(0.15388547) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34129225) q[0];
sx q[0];
rz(-1.253739) q[0];
sx q[0];
rz(-1.6503667) q[0];
x q[1];
rz(-2.4272404) q[2];
sx q[2];
rz(-1.0754943) q[2];
sx q[2];
rz(-2.60204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9388403) q[1];
sx q[1];
rz(-1.2766663) q[1];
sx q[1];
rz(0.37962706) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4948972) q[3];
sx q[3];
rz(-1.5468949) q[3];
sx q[3];
rz(2.0233512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1880356) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(2.6970862) q[2];
rz(0.19017531) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1403777) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(-1.0706527) q[0];
rz(0.045914687) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(2.3505223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.175638) q[0];
sx q[0];
rz(-1.6583867) q[0];
sx q[0];
rz(-1.9071294) q[0];
rz(-pi) q[1];
rz(2.5612381) q[2];
sx q[2];
rz(-0.46541801) q[2];
sx q[2];
rz(-0.67829715) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2998878) q[1];
sx q[1];
rz(-1.1763402) q[1];
sx q[1];
rz(-1.9576555) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8077352) q[3];
sx q[3];
rz(-0.51722368) q[3];
sx q[3];
rz(-2.5631529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(-1.4727288) q[2];
rz(-1.6139) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2033757) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(-2.6279502) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(1.8088874) q[2];
sx q[2];
rz(-1.0793964) q[2];
sx q[2];
rz(2.4894077) q[2];
rz(-0.41889965) q[3];
sx q[3];
rz(-2.2875026) q[3];
sx q[3];
rz(0.84411375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
