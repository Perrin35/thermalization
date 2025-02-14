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
rz(1.3266069) q[0];
sx q[0];
rz(3.366037) q[0];
sx q[0];
rz(11.233575) q[0];
rz(-0.83238554) q[1];
sx q[1];
rz(4.8424911) q[1];
sx q[1];
rz(10.535156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8081759) q[0];
sx q[0];
rz(-1.7000755) q[0];
sx q[0];
rz(-0.16462932) q[0];
rz(-pi) q[1];
rz(-0.30958561) q[2];
sx q[2];
rz(-0.88304115) q[2];
sx q[2];
rz(1.0716455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79506766) q[1];
sx q[1];
rz(-1.5988886) q[1];
sx q[1];
rz(-0.053140716) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2625804) q[3];
sx q[3];
rz(-2.1860414) q[3];
sx q[3];
rz(-0.32549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2276459) q[2];
sx q[2];
rz(-0.0089184428) q[2];
sx q[2];
rz(-1.439636) q[2];
rz(-1.4134183) q[3];
sx q[3];
rz(-0.011912502) q[3];
sx q[3];
rz(2.8720776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.444376) q[0];
sx q[0];
rz(-1.603729) q[0];
sx q[0];
rz(0.61475301) q[0];
rz(0.53601021) q[1];
sx q[1];
rz(-0.025608048) q[1];
sx q[1];
rz(-0.33686179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1/(6*pi)) q[0];
sx q[0];
rz(-1.7053002) q[0];
sx q[0];
rz(0.10945871) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.031177715) q[2];
sx q[2];
rz(-0.47652361) q[2];
sx q[2];
rz(-0.88821238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.555111) q[1];
sx q[1];
rz(-0.046849061) q[1];
sx q[1];
rz(1.2270801) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9514489) q[3];
sx q[3];
rz(-1.8845673) q[3];
sx q[3];
rz(-2.6699175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2577995) q[2];
sx q[2];
rz(-0.012848583) q[2];
sx q[2];
rz(0.69965714) q[2];
rz(2.9720225) q[3];
sx q[3];
rz(-0.01378672) q[3];
sx q[3];
rz(0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4366142) q[0];
sx q[0];
rz(-0.545937) q[0];
sx q[0];
rz(-2.8563232) q[0];
rz(-2.8778695) q[1];
sx q[1];
rz(-3.141005) q[1];
sx q[1];
rz(-2.347351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4141615) q[0];
sx q[0];
rz(-2.2530956) q[0];
sx q[0];
rz(-1.4201565) q[0];
x q[1];
rz(1.8828234) q[2];
sx q[2];
rz(-2.8030382) q[2];
sx q[2];
rz(-2.3046913) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7020828) q[1];
sx q[1];
rz(-3.1065662) q[1];
sx q[1];
rz(-0.18478025) q[1];
rz(-0.66624347) q[3];
sx q[3];
rz(-1.0292205) q[3];
sx q[3];
rz(0.2592087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11702015) q[2];
sx q[2];
rz(-3.0954376) q[2];
sx q[2];
rz(-0.864492) q[2];
rz(-0.74508673) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(0.043070506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87557781) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(-2.2752046) q[0];
rz(0.59151793) q[1];
sx q[1];
rz(-0.94647995) q[1];
sx q[1];
rz(2.0902858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503904) q[0];
sx q[0];
rz(-1.5728381) q[0];
sx q[0];
rz(-1.5709086) q[0];
x q[1];
rz(0.00045108081) q[2];
sx q[2];
rz(-1.5685097) q[2];
sx q[2];
rz(-3.1368983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2308553) q[1];
sx q[1];
rz(-1.1713109) q[1];
sx q[1];
rz(2.3302156) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0962001) q[3];
sx q[3];
rz(-2.0936493) q[3];
sx q[3];
rz(-0.16663545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6136578) q[2];
sx q[2];
rz(-0.064585678) q[2];
sx q[2];
rz(-2.2876372) q[2];
rz(0.69771403) q[3];
sx q[3];
rz(-1.8055975) q[3];
sx q[3];
rz(0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39603221) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(-2.7295617) q[0];
rz(1.283006) q[1];
sx q[1];
rz(-2.367815) q[1];
sx q[1];
rz(2.3903019) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728516) q[0];
sx q[0];
rz(-1.6882467) q[0];
sx q[0];
rz(2.9020578) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22789006) q[2];
sx q[2];
rz(-2.2376559) q[2];
sx q[2];
rz(1.6713072) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3739864) q[1];
sx q[1];
rz(-2.3298378) q[1];
sx q[1];
rz(-2.4381258) q[1];
rz(-pi) q[2];
rz(0.86814084) q[3];
sx q[3];
rz(-2.9538547) q[3];
sx q[3];
rz(2.7526698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6777307) q[2];
sx q[2];
rz(-0.025721392) q[2];
sx q[2];
rz(0.99979293) q[2];
rz(2.540588) q[3];
sx q[3];
rz(-0.060636245) q[3];
sx q[3];
rz(2.291688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8782225) q[0];
sx q[0];
rz(-3.0510986) q[0];
sx q[0];
rz(-2.7409842) q[0];
rz(-1.7349617) q[1];
sx q[1];
rz(-2.7901283) q[1];
sx q[1];
rz(1.7469143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050848518) q[0];
sx q[0];
rz(-1.5511449) q[0];
sx q[0];
rz(-0.0046363884) q[0];
x q[1];
rz(-1.2591061) q[2];
sx q[2];
rz(-2.3505728) q[2];
sx q[2];
rz(-0.64779161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7575824) q[1];
sx q[1];
rz(-1.7613125) q[1];
sx q[1];
rz(-1.5662929) q[1];
rz(-2.9506545) q[3];
sx q[3];
rz(-0.74511601) q[3];
sx q[3];
rz(2.7332954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74581915) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(-2.0664717) q[2];
rz(-2.6221258) q[3];
sx q[3];
rz(-0.1778917) q[3];
sx q[3];
rz(1.547267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2125242) q[0];
sx q[0];
rz(-1.9269749) q[0];
sx q[0];
rz(1.6125096) q[0];
rz(-2.5802338) q[1];
sx q[1];
rz(-3.1293271) q[1];
sx q[1];
rz(-0.57048172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5382696) q[0];
sx q[0];
rz(-0.13165671) q[0];
sx q[0];
rz(-1.9719453) q[0];
x q[1];
rz(0.95593234) q[2];
sx q[2];
rz(-2.758965) q[2];
sx q[2];
rz(1.0997694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9614547) q[1];
sx q[1];
rz(-1.5552863) q[1];
sx q[1];
rz(0.0011891459) q[1];
rz(-pi) q[2];
rz(-0.99659427) q[3];
sx q[3];
rz(-0.83296597) q[3];
sx q[3];
rz(1.0723237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.059171112) q[2];
sx q[2];
rz(-2.8736281) q[2];
sx q[2];
rz(1.9783665) q[2];
rz(-1.5393114) q[3];
sx q[3];
rz(-3.0395165) q[3];
sx q[3];
rz(-0.99360895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7546643) q[0];
sx q[0];
rz(-0.29351497) q[0];
sx q[0];
rz(-0.64025229) q[0];
rz(-2.2786268) q[1];
sx q[1];
rz(-2.9997365) q[1];
sx q[1];
rz(-0.77756768) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3024738) q[0];
sx q[0];
rz(-0.57838744) q[0];
sx q[0];
rz(-2.3443787) q[0];
rz(1.338358) q[2];
sx q[2];
rz(-0.82045499) q[2];
sx q[2];
rz(-2.6636843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9915087) q[1];
sx q[1];
rz(-1.5822295) q[1];
sx q[1];
rz(1.4978133) q[1];
rz(-pi) q[2];
rz(-0.74109091) q[3];
sx q[3];
rz(-2.5020385) q[3];
sx q[3];
rz(0.049413817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61206478) q[2];
sx q[2];
rz(-2.7175588) q[2];
sx q[2];
rz(0.5303793) q[2];
rz(3.1139167) q[3];
sx q[3];
rz(-0.045462463) q[3];
sx q[3];
rz(-1.8224705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46488047) q[0];
sx q[0];
rz(-2.9806529) q[0];
sx q[0];
rz(-0.93223923) q[0];
rz(-1.5981916) q[1];
sx q[1];
rz(-2.3379969) q[1];
sx q[1];
rz(0.34206051) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0823183) q[0];
sx q[0];
rz(-3.0558028) q[0];
sx q[0];
rz(2.9640574) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6342966) q[2];
sx q[2];
rz(-2.3076153) q[2];
sx q[2];
rz(-0.33499107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4885713) q[1];
sx q[1];
rz(-2.584051) q[1];
sx q[1];
rz(-2.7776021) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4985648) q[3];
sx q[3];
rz(-2.9119303) q[3];
sx q[3];
rz(0.1230474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0590234) q[2];
sx q[2];
rz(-0.00073585357) q[2];
sx q[2];
rz(-1.9334582) q[2];
rz(-0.9515323) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(2.4212196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35219881) q[0];
sx q[0];
rz(-0.77704) q[0];
sx q[0];
rz(3.0598031) q[0];
rz(-2.8261322) q[1];
sx q[1];
rz(-0.052736484) q[1];
sx q[1];
rz(-1.9116521) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30222505) q[0];
sx q[0];
rz(-1.7291843) q[0];
sx q[0];
rz(-1.4119638) q[0];
rz(2.4794078) q[2];
sx q[2];
rz(-2.9176468) q[2];
sx q[2];
rz(2.0020773) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1217864) q[1];
sx q[1];
rz(-2.9933194) q[1];
sx q[1];
rz(2.2231977) q[1];
rz(1.337215) q[3];
sx q[3];
rz(-2.5982435) q[3];
sx q[3];
rz(-0.41239967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6007467) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(1.1136327) q[2];
rz(3.1019548) q[3];
sx q[3];
rz(-0.0097291917) q[3];
sx q[3];
rz(-0.59687328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.6780846) q[0];
sx q[0];
rz(-1.2152553) q[0];
sx q[0];
rz(-1.8334462) q[0];
rz(-2.5254163) q[1];
sx q[1];
rz(-1.0721075) q[1];
sx q[1];
rz(0.20872605) q[1];
rz(0.19113541) q[2];
sx q[2];
rz(-1.8046817) q[2];
sx q[2];
rz(2.6006151) q[2];
rz(-3.1022648) q[3];
sx q[3];
rz(-1.3500742) q[3];
sx q[3];
rz(-3.1108656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
