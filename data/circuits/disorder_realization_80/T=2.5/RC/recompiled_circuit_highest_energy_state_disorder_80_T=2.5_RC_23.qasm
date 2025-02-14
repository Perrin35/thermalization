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
rz(0.12240527) q[0];
sx q[0];
rz(2.2221017) q[0];
sx q[0];
rz(10.623951) q[0];
rz(0.1872669) q[1];
sx q[1];
rz(-2.5993102) q[1];
sx q[1];
rz(-1.5195001) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9334979) q[0];
sx q[0];
rz(-1.6483432) q[0];
sx q[0];
rz(0.61570517) q[0];
rz(0.89251065) q[2];
sx q[2];
rz(-2.5830088) q[2];
sx q[2];
rz(0.19772274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28827661) q[1];
sx q[1];
rz(-1.1697993) q[1];
sx q[1];
rz(-0.53944352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19272645) q[3];
sx q[3];
rz(-0.24860074) q[3];
sx q[3];
rz(-2.1603487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8444933) q[2];
sx q[2];
rz(-0.62701925) q[2];
sx q[2];
rz(1.373488) q[2];
rz(2.3376236) q[3];
sx q[3];
rz(-1.6425902) q[3];
sx q[3];
rz(1.9544301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74541575) q[0];
sx q[0];
rz(-0.71393037) q[0];
sx q[0];
rz(-0.3748689) q[0];
rz(-0.51775852) q[1];
sx q[1];
rz(-1.2034143) q[1];
sx q[1];
rz(-1.7727324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33297172) q[0];
sx q[0];
rz(-0.00095168984) q[0];
sx q[0];
rz(3.0406221) q[0];
x q[1];
rz(-0.19646074) q[2];
sx q[2];
rz(-1.5443373) q[2];
sx q[2];
rz(0.65633869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9437286) q[1];
sx q[1];
rz(-1.9790589) q[1];
sx q[1];
rz(0.54710435) q[1];
rz(-pi) q[2];
rz(-1.7852704) q[3];
sx q[3];
rz(-1.4376336) q[3];
sx q[3];
rz(-1.5758297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0181197) q[2];
sx q[2];
rz(-0.70088434) q[2];
sx q[2];
rz(2.4269721) q[2];
rz(0.2054275) q[3];
sx q[3];
rz(-1.3222062) q[3];
sx q[3];
rz(1.8429168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648723) q[0];
sx q[0];
rz(-2.1045852) q[0];
sx q[0];
rz(-0.49764693) q[0];
rz(-0.63703713) q[1];
sx q[1];
rz(-0.62890816) q[1];
sx q[1];
rz(3.0348437) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15539385) q[0];
sx q[0];
rz(-1.8698911) q[0];
sx q[0];
rz(0.22217447) q[0];
x q[1];
rz(-1.5791513) q[2];
sx q[2];
rz(-2.4718577) q[2];
sx q[2];
rz(2.1363897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1615352) q[1];
sx q[1];
rz(-1.5173459) q[1];
sx q[1];
rz(-2.2845099) q[1];
rz(-0.40092881) q[3];
sx q[3];
rz(-1.2967921) q[3];
sx q[3];
rz(-1.9131434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5367624) q[2];
sx q[2];
rz(-2.3894775) q[2];
sx q[2];
rz(2.9467648) q[2];
rz(3.1365862) q[3];
sx q[3];
rz(-2.5275793) q[3];
sx q[3];
rz(-0.79202882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0358589) q[0];
sx q[0];
rz(-0.53426131) q[0];
sx q[0];
rz(2.1015097) q[0];
rz(2.9478574) q[1];
sx q[1];
rz(-1.5015142) q[1];
sx q[1];
rz(-0.74660444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116007) q[0];
sx q[0];
rz(-2.009729) q[0];
sx q[0];
rz(-0.12719391) q[0];
x q[1];
rz(1.7393635) q[2];
sx q[2];
rz(-2.4344846) q[2];
sx q[2];
rz(-0.76729028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0185701) q[1];
sx q[1];
rz(-2.445652) q[1];
sx q[1];
rz(3.0607515) q[1];
rz(-pi) q[2];
rz(-1.6799404) q[3];
sx q[3];
rz(-1.8092844) q[3];
sx q[3];
rz(1.3772688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1497583) q[2];
sx q[2];
rz(-1.4703625) q[2];
sx q[2];
rz(0.20450083) q[2];
rz(-1.1931984) q[3];
sx q[3];
rz(-2.2843993) q[3];
sx q[3];
rz(-2.1804555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1292773) q[0];
sx q[0];
rz(-0.17450541) q[0];
sx q[0];
rz(-1.0804863) q[0];
rz(2.7642545) q[1];
sx q[1];
rz(-1.6308547) q[1];
sx q[1];
rz(-0.89471716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0102826) q[0];
sx q[0];
rz(-1.9246411) q[0];
sx q[0];
rz(1.8132339) q[0];
rz(-pi) q[1];
rz(0.9099877) q[2];
sx q[2];
rz(-0.89471451) q[2];
sx q[2];
rz(-2.8582339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1929568) q[1];
sx q[1];
rz(-1.0851477) q[1];
sx q[1];
rz(-0.57106496) q[1];
x q[2];
rz(-1.0977488) q[3];
sx q[3];
rz(-0.69907197) q[3];
sx q[3];
rz(0.11156532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47199029) q[2];
sx q[2];
rz(-0.44595465) q[2];
sx q[2];
rz(-0.10776821) q[2];
rz(-2.9660411) q[3];
sx q[3];
rz(-1.2551509) q[3];
sx q[3];
rz(-0.56041437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046722978) q[0];
sx q[0];
rz(-0.51114285) q[0];
sx q[0];
rz(-2.129659) q[0];
rz(-1.5049505) q[1];
sx q[1];
rz(-2.4220146) q[1];
sx q[1];
rz(-2.124427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0295231) q[0];
sx q[0];
rz(-1.3991303) q[0];
sx q[0];
rz(-1.8708234) q[0];
x q[1];
rz(-2.5565113) q[2];
sx q[2];
rz(-0.85137109) q[2];
sx q[2];
rz(0.78851267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1348844) q[1];
sx q[1];
rz(-1.5794288) q[1];
sx q[1];
rz(2.951202) q[1];
x q[2];
rz(0.15729842) q[3];
sx q[3];
rz(-1.1264115) q[3];
sx q[3];
rz(-0.97287662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4096058) q[2];
sx q[2];
rz(-0.66183949) q[2];
sx q[2];
rz(0.99676639) q[2];
rz(3.0766727) q[3];
sx q[3];
rz(-1.4109979) q[3];
sx q[3];
rz(1.9916649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0885334) q[0];
sx q[0];
rz(-0.94045883) q[0];
sx q[0];
rz(0.9084107) q[0];
rz(0.21993318) q[1];
sx q[1];
rz(-1.5989774) q[1];
sx q[1];
rz(-1.8904846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.661099) q[0];
sx q[0];
rz(-1.1014043) q[0];
sx q[0];
rz(-0.29222699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.471453) q[2];
sx q[2];
rz(-1.8548428) q[2];
sx q[2];
rz(2.5104475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0323769) q[1];
sx q[1];
rz(-2.2046979) q[1];
sx q[1];
rz(-0.79297592) q[1];
rz(-pi) q[2];
rz(1.7102107) q[3];
sx q[3];
rz(-2.406139) q[3];
sx q[3];
rz(0.212634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0327586) q[2];
sx q[2];
rz(-2.3075576) q[2];
sx q[2];
rz(1.533482) q[2];
rz(1.9882625) q[3];
sx q[3];
rz(-1.0554375) q[3];
sx q[3];
rz(-1.4920894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.15025) q[0];
sx q[0];
rz(-2.7582176) q[0];
sx q[0];
rz(-2.2008994) q[0];
rz(3.0062145) q[1];
sx q[1];
rz(-1.7372513) q[1];
sx q[1];
rz(-2.1167596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.301046) q[0];
sx q[0];
rz(-2.4468166) q[0];
sx q[0];
rz(2.0639117) q[0];
rz(1.5103136) q[2];
sx q[2];
rz(-0.96394682) q[2];
sx q[2];
rz(1.902193) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7861745) q[1];
sx q[1];
rz(-0.23705951) q[1];
sx q[1];
rz(2.5393344) q[1];
rz(0.41127326) q[3];
sx q[3];
rz(-1.5233293) q[3];
sx q[3];
rz(-2.2293343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3465053) q[2];
sx q[2];
rz(-2.2515209) q[2];
sx q[2];
rz(0.65004641) q[2];
rz(0.58642379) q[3];
sx q[3];
rz(-1.6610049) q[3];
sx q[3];
rz(-0.75064269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.2670249) q[0];
sx q[0];
rz(-2.6331007) q[0];
sx q[0];
rz(1.9961927) q[0];
rz(1.3151431) q[1];
sx q[1];
rz(-1.9086842) q[1];
sx q[1];
rz(-2.4536536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8696339) q[0];
sx q[0];
rz(-0.42068538) q[0];
sx q[0];
rz(-2.0482333) q[0];
rz(1.0264977) q[2];
sx q[2];
rz(-0.67152464) q[2];
sx q[2];
rz(1.6467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44837828) q[1];
sx q[1];
rz(-2.3168987) q[1];
sx q[1];
rz(-1.34971) q[1];
x q[2];
rz(2.2086772) q[3];
sx q[3];
rz(-2.5543) q[3];
sx q[3];
rz(-1.6577394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73100662) q[2];
sx q[2];
rz(-2.035391) q[2];
sx q[2];
rz(0.86622396) q[2];
rz(0.90803641) q[3];
sx q[3];
rz(-0.76483813) q[3];
sx q[3];
rz(-2.0456555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753801) q[0];
sx q[0];
rz(-1.9575653) q[0];
sx q[0];
rz(2.7259735) q[0];
rz(-2.3497154) q[1];
sx q[1];
rz(-1.917058) q[1];
sx q[1];
rz(-0.22722879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89657518) q[0];
sx q[0];
rz(-1.554721) q[0];
sx q[0];
rz(1.5825558) q[0];
x q[1];
rz(1.4899026) q[2];
sx q[2];
rz(-2.3757907) q[2];
sx q[2];
rz(-1.3882989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38293441) q[1];
sx q[1];
rz(-0.99731748) q[1];
sx q[1];
rz(1.2632779) q[1];
x q[2];
rz(-2.4070508) q[3];
sx q[3];
rz(-2.4895222) q[3];
sx q[3];
rz(-1.2800686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8772584) q[2];
sx q[2];
rz(-2.4206968) q[2];
sx q[2];
rz(2.7016675) q[2];
rz(-2.544493) q[3];
sx q[3];
rz(-0.66949451) q[3];
sx q[3];
rz(1.3574903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4933585) q[0];
sx q[0];
rz(-1.5536722) q[0];
sx q[0];
rz(-1.8552725) q[0];
rz(1.7276806) q[1];
sx q[1];
rz(-1.021011) q[1];
sx q[1];
rz(0.11722142) q[1];
rz(2.734388) q[2];
sx q[2];
rz(-1.5199678) q[2];
sx q[2];
rz(0.54553568) q[2];
rz(2.8942378) q[3];
sx q[3];
rz(-2.3079899) q[3];
sx q[3];
rz(0.079903825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
