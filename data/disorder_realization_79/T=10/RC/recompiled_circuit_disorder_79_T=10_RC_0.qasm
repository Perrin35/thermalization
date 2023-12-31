OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642759) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(3.0741865) q[0];
rz(0.97857742) q[2];
sx q[2];
rz(-2.7263612) q[2];
sx q[2];
rz(-0.65939553) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4343623) q[1];
sx q[1];
rz(-2.011236) q[1];
sx q[1];
rz(0.67727725) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1300415) q[3];
sx q[3];
rz(-0.64239255) q[3];
sx q[3];
rz(-1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(2.0430298) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.6764486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6786878) q[0];
sx q[0];
rz(-1.5683187) q[0];
sx q[0];
rz(0.72420995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0589774) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(-0.28085923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1190471) q[1];
sx q[1];
rz(-0.76821583) q[1];
sx q[1];
rz(0.70336282) q[1];
rz(-1.2259237) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8974509) q[0];
sx q[0];
rz(-1.7139009) q[0];
sx q[0];
rz(1.1600526) q[0];
x q[1];
rz(0.82528798) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(0.7691783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.264818) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(0.012621565) q[1];
rz(-pi) q[2];
rz(2.839746) q[3];
sx q[3];
rz(-2.2491124) q[3];
sx q[3];
rz(-0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.7165855) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(-1.0850614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59180075) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(-0.53996284) q[0];
rz(-0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(-2.0476066) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9492053) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(-2.5237571) q[1];
x q[2];
rz(-2.5892341) q[3];
sx q[3];
rz(-0.68867749) q[3];
sx q[3];
rz(-2.3104582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(2.4604649) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-2.1599105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0504793) q[0];
sx q[0];
rz(-2.1149153) q[0];
sx q[0];
rz(-1.0247466) q[0];
rz(-pi) q[1];
rz(-2.1443411) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(-2.2821102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.423646) q[1];
sx q[1];
rz(-1.4232993) q[1];
sx q[1];
rz(-0.46848483) q[1];
x q[2];
rz(1.6793628) q[3];
sx q[3];
rz(-0.6001937) q[3];
sx q[3];
rz(-1.997228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.3283407) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(-1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-0.16528027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344791) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(1.5414184) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6406306) q[2];
sx q[2];
rz(-2.2572821) q[2];
sx q[2];
rz(-1.7718466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.936441) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(-2.5564297) q[1];
rz(-pi) q[2];
rz(1.1214439) q[3];
sx q[3];
rz(-1.3580139) q[3];
sx q[3];
rz(1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(-0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1897141) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(-0.017107054) q[0];
x q[1];
rz(-1.9095124) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(-0.92600694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2939824) q[1];
sx q[1];
rz(-0.81969205) q[1];
sx q[1];
rz(1.3241029) q[1];
rz(1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(-1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(-2.2198026) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(0.69358967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.2459754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31729749) q[1];
sx q[1];
rz(-1.824914) q[1];
sx q[1];
rz(1.7829249) q[1];
x q[2];
rz(0.17955762) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(-1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.3681108) q[0];
sx q[0];
rz(-0.52635877) q[0];
rz(0.87551261) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(-1.306844) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(-0.70396522) q[1];
x q[2];
rz(0.67919517) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1264964) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-2.2311907) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(-1.1766599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991853) q[0];
sx q[0];
rz(-0.17051324) q[0];
sx q[0];
rz(-1.8605581) q[0];
rz(-pi) q[1];
rz(2.2048336) q[2];
sx q[2];
rz(-2.0271795) q[2];
sx q[2];
rz(-0.87181834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6224222) q[1];
sx q[1];
rz(-1.4960438) q[1];
sx q[1];
rz(2.1291817) q[1];
rz(-pi) q[2];
rz(-0.15575274) q[3];
sx q[3];
rz(-2.6253346) q[3];
sx q[3];
rz(1.2092276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(0.30504967) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(0.022396537) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
