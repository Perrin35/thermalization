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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17380938) q[0];
sx q[0];
rz(-0.73497226) q[0];
sx q[0];
rz(1.3761139) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9453968) q[2];
sx q[2];
rz(-1.2805403) q[2];
sx q[2];
rz(-2.6930075) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28847028) q[1];
sx q[1];
rz(-0.56325699) q[1];
sx q[1];
rz(2.017615) q[1];
rz(2.8075571) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(-0.087466784) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66602089) q[0];
sx q[0];
rz(-0.7374987) q[0];
sx q[0];
rz(-2.0185508) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42627724) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(-0.041235812) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14428917) q[1];
sx q[1];
rz(-1.803949) q[1];
sx q[1];
rz(-2.7551329) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53701138) q[3];
sx q[3];
rz(-1.7039263) q[3];
sx q[3];
rz(2.6739124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(-0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(0.18951167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127842) q[0];
sx q[0];
rz(-1.841294) q[0];
sx q[0];
rz(-2.4012647) q[0];
rz(-pi) q[1];
rz(-1.465221) q[2];
sx q[2];
rz(-1.1974679) q[2];
sx q[2];
rz(0.8651274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.0057919766) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(-0.23975753) q[1];
rz(-pi) q[2];
rz(2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(-0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074993) q[0];
sx q[0];
rz(-0.543569) q[0];
sx q[0];
rz(2.4242633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3218669) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(2.0476598) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66149368) q[1];
sx q[1];
rz(-1.1575932) q[1];
sx q[1];
rz(-2.5109629) q[1];
rz(-pi) q[2];
rz(-0.43153901) q[3];
sx q[3];
rz(-0.46436897) q[3];
sx q[3];
rz(-2.4114256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-2.1676846) q[0];
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
rz(-2.1954384) q[0];
sx q[0];
rz(2.7420068) q[0];
rz(2.2511475) q[2];
sx q[2];
rz(-1.3157985) q[2];
sx q[2];
rz(1.4444218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3586732) q[1];
sx q[1];
rz(-1.9649319) q[1];
sx q[1];
rz(-0.47788099) q[1];
x q[2];
rz(1.3768457) q[3];
sx q[3];
rz(-1.0411106) q[3];
sx q[3];
rz(2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(2.2987135) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(-2.4687016) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653172) q[0];
sx q[0];
rz(-0.99940171) q[0];
sx q[0];
rz(-3.0185643) q[0];
x q[1];
rz(-2.2851373) q[2];
sx q[2];
rz(-1.5521142) q[2];
sx q[2];
rz(-3.0532388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3402965) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(1.3188386) q[1];
x q[2];
rz(-2.9665885) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(-1.7815636) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.7395082) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8334478) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(-1.863198) q[0];
rz(-pi) q[1];
rz(-2.8962171) q[2];
sx q[2];
rz(-0.78005314) q[2];
sx q[2];
rz(-1.1749554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7359888) q[1];
sx q[1];
rz(-2.7466008) q[1];
sx q[1];
rz(0.40785457) q[1];
rz(-3.0038463) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.54135281) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68529785) q[0];
sx q[0];
rz(-2.1349499) q[0];
sx q[0];
rz(-1.5575404) q[0];
rz(-1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(0.18061772) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9882422) q[1];
sx q[1];
rz(-1.9820947) q[1];
sx q[1];
rz(0.78505713) q[1];
rz(-pi) q[2];
rz(-0.0046644966) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(-0.62957803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(-2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709258) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(-2.4404581) q[0];
x q[1];
rz(0.88405774) q[2];
sx q[2];
rz(-0.80342711) q[2];
sx q[2];
rz(0.49056177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0874789) q[1];
sx q[1];
rz(-0.94062727) q[1];
sx q[1];
rz(-1.5515045) q[1];
rz(-pi) q[2];
rz(2.290756) q[3];
sx q[3];
rz(-2.3265504) q[3];
sx q[3];
rz(1.1834061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33567515) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(-2.5277396) q[0];
rz(-2.0954258) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(-2.7493431) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0961571) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(-0.57604726) q[0];
rz(-pi) q[1];
rz(2.4105083) q[2];
sx q[2];
rz(-2.3713171) q[2];
sx q[2];
rz(-1.1243656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0543921) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(2.2589576) q[1];
rz(-pi) q[2];
rz(1.5997821) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7077211) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(2.375013) q[2];
sx q[2];
rz(-1.8324413) q[2];
sx q[2];
rz(1.2698297) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];