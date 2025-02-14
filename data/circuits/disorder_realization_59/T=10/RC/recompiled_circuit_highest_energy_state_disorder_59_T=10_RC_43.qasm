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
rz(1.264313) q[0];
sx q[0];
rz(6.6067196) q[0];
sx q[0];
rz(8.9759448) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(1.356025) q[1];
sx q[1];
rz(7.4577509) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313707) q[0];
sx q[0];
rz(-2.1315984) q[0];
sx q[0];
rz(-1.9813374) q[0];
rz(2.9401112) q[2];
sx q[2];
rz(-0.32989943) q[2];
sx q[2];
rz(-0.13452521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.9372896) q[1];
sx q[1];
rz(-2.1912592) q[1];
sx q[1];
rz(-3.0985249) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9644792) q[3];
sx q[3];
rz(-0.52059725) q[3];
sx q[3];
rz(0.89730747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.590098) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-2.2287915) q[2];
rz(-1.270795) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(-0.84996581) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368644) q[0];
sx q[0];
rz(-0.57924634) q[0];
sx q[0];
rz(2.7194523) q[0];
rz(-2.7719851) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(-0.94013989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5053609) q[0];
sx q[0];
rz(-2.8121037) q[0];
sx q[0];
rz(1.5630653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11905382) q[2];
sx q[2];
rz(-1.040689) q[2];
sx q[2];
rz(1.334545) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9782422) q[1];
sx q[1];
rz(-2.2903195) q[1];
sx q[1];
rz(0.6759599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89408447) q[3];
sx q[3];
rz(-1.7853123) q[3];
sx q[3];
rz(2.0443975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1993572) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(2.2352236) q[2];
rz(2.1622315) q[3];
sx q[3];
rz(-0.40458471) q[3];
sx q[3];
rz(-1.2649068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449428) q[0];
sx q[0];
rz(-0.066983797) q[0];
sx q[0];
rz(1.8970733) q[0];
rz(-2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-0.32274524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30814442) q[0];
sx q[0];
rz(-2.2287773) q[0];
sx q[0];
rz(-1.4288155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8333146) q[2];
sx q[2];
rz(-1.0325422) q[2];
sx q[2];
rz(2.8652712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6420329) q[1];
sx q[1];
rz(-2.316058) q[1];
sx q[1];
rz(-2.8241629) q[1];
rz(-pi) q[2];
rz(2.7904872) q[3];
sx q[3];
rz(-1.5241844) q[3];
sx q[3];
rz(-0.24735115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9129703) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(1.7886394) q[2];
rz(-0.27255034) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6005741) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(1.898265) q[0];
rz(-0.86878949) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(-1.0702081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58067034) q[0];
sx q[0];
rz(-1.0698478) q[0];
sx q[0];
rz(0.83606677) q[0];
x q[1];
rz(-1.1203499) q[2];
sx q[2];
rz(-2.3131158) q[2];
sx q[2];
rz(0.38249215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3096042) q[1];
sx q[1];
rz(-1.2491396) q[1];
sx q[1];
rz(3.020806) q[1];
rz(0.34262212) q[3];
sx q[3];
rz(-1.7512158) q[3];
sx q[3];
rz(-2.8460549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5493245) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(2.2198524) q[2];
rz(-0.24374572) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(0.29269472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.9859966) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(-0.50088125) q[0];
rz(-2.0241375) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(-0.31455988) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4582493) q[0];
sx q[0];
rz(-1.9106081) q[0];
sx q[0];
rz(-0.37796867) q[0];
rz(-pi) q[1];
rz(-0.0044745046) q[2];
sx q[2];
rz(-1.6142077) q[2];
sx q[2];
rz(-0.89148607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3479473) q[1];
sx q[1];
rz(-2.1160908) q[1];
sx q[1];
rz(0.46635638) q[1];
rz(-3.0231832) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(1.6550696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.306281) q[2];
sx q[2];
rz(-2.0302782) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(0.48526192) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(-3.0418292) q[0];
rz(-1.2707204) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(1.7521923) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64369624) q[0];
sx q[0];
rz(-0.31539105) q[0];
sx q[0];
rz(2.5853378) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6011141) q[2];
sx q[2];
rz(-2.9394144) q[2];
sx q[2];
rz(2.0433189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2459316) q[1];
sx q[1];
rz(-0.97890039) q[1];
sx q[1];
rz(-2.3227957) q[1];
rz(-1.4445969) q[3];
sx q[3];
rz(-0.79168237) q[3];
sx q[3];
rz(2.0711957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(0.89034447) q[2];
rz(-1.9503615) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(2.2395535) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(-3.1088767) q[0];
rz(2.6383767) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(3.018766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3204725) q[0];
sx q[0];
rz(-0.80424612) q[0];
sx q[0];
rz(0.45566092) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2314447) q[2];
sx q[2];
rz(-2.0464346) q[2];
sx q[2];
rz(0.95107691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6325153) q[1];
sx q[1];
rz(-1.6700524) q[1];
sx q[1];
rz(-2.0429913) q[1];
rz(-0.34103532) q[3];
sx q[3];
rz(-0.37217316) q[3];
sx q[3];
rz(0.91590524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7922625) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(1.9897423) q[2];
rz(2.2564015) q[3];
sx q[3];
rz(-1.7722292) q[3];
sx q[3];
rz(-1.4256029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952878) q[0];
sx q[0];
rz(-0.9335683) q[0];
sx q[0];
rz(0.92920148) q[0];
rz(0.62905351) q[1];
sx q[1];
rz(-1.7864497) q[1];
sx q[1];
rz(0.74310511) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72637314) q[0];
sx q[0];
rz(-0.45436828) q[0];
sx q[0];
rz(1.7512064) q[0];
rz(-2.9628108) q[2];
sx q[2];
rz(-0.4950854) q[2];
sx q[2];
rz(2.924233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9070702) q[1];
sx q[1];
rz(-1.6369695) q[1];
sx q[1];
rz(-2.727134) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1984897) q[3];
sx q[3];
rz(-2.4022958) q[3];
sx q[3];
rz(-3.0969248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5762081) q[2];
sx q[2];
rz(-0.78018633) q[2];
sx q[2];
rz(0.72594491) q[2];
rz(-2.7706326) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(-2.7788739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87042701) q[0];
sx q[0];
rz(-0.83036625) q[0];
sx q[0];
rz(0.56394947) q[0];
rz(-1.14934) q[1];
sx q[1];
rz(-0.96200689) q[1];
sx q[1];
rz(-0.74251485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1127204) q[0];
sx q[0];
rz(-2.0322587) q[0];
sx q[0];
rz(2.0966778) q[0];
rz(-pi) q[1];
rz(-0.60301247) q[2];
sx q[2];
rz(-0.49076395) q[2];
sx q[2];
rz(-0.025207504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2143605) q[1];
sx q[1];
rz(-0.4402658) q[1];
sx q[1];
rz(0.10378598) q[1];
x q[2];
rz(-2.2015195) q[3];
sx q[3];
rz(-2.3204062) q[3];
sx q[3];
rz(-1.8668777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62770647) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(0.34454301) q[2];
rz(-0.44376093) q[3];
sx q[3];
rz(-1.0756451) q[3];
sx q[3];
rz(1.8189323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(2.479582) q[0];
rz(2.8061891) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(-0.9368771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360817) q[0];
sx q[0];
rz(-0.72692211) q[0];
sx q[0];
rz(1.0127823) q[0];
rz(-pi) q[1];
rz(1.4858021) q[2];
sx q[2];
rz(-0.27117816) q[2];
sx q[2];
rz(-0.32395872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7171927) q[1];
sx q[1];
rz(-0.48312995) q[1];
sx q[1];
rz(2.4880954) q[1];
rz(2.9659445) q[3];
sx q[3];
rz(-1.3743168) q[3];
sx q[3];
rz(2.548747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-2.9073538) q[2];
sx q[2];
rz(0.49918175) q[2];
rz(-1.4623803) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(1.6818989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41596169) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(2.536643) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(-1.4251322) q[2];
sx q[2];
rz(-2.1389516) q[2];
sx q[2];
rz(0.4898796) q[2];
rz(2.6400186) q[3];
sx q[3];
rz(-1.4091103) q[3];
sx q[3];
rz(1.4179358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
