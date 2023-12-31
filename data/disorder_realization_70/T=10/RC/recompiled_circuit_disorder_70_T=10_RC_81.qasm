OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(-2.5075066) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4579826) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(-0.81575127) q[0];
rz(-pi) q[1];
rz(3.0787266) q[2];
sx q[2];
rz(-0.91744423) q[2];
sx q[2];
rz(-1.9444998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7635599) q[1];
sx q[1];
rz(-2.0114007) q[1];
sx q[1];
rz(-0.9159169) q[1];
x q[2];
rz(-1.5473458) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.263164) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-2.8485956) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(0.43757004) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76925812) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.9702205) q[0];
rz(-pi) q[1];
rz(-1.2220076) q[2];
sx q[2];
rz(-1.6752599) q[2];
sx q[2];
rz(-1.2806569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.931103) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(-1.7060075) q[1];
rz(-pi) q[2];
rz(1.4731746) q[3];
sx q[3];
rz(-2.0711581) q[3];
sx q[3];
rz(-0.2122768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-2.5533) q[2];
rz(0.44899392) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(-2.589059) q[0];
rz(-2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(-2.8734145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0231087) q[1];
sx q[1];
rz(-2.1493836) q[1];
sx q[1];
rz(-0.56225496) q[1];
x q[2];
rz(-2.5030701) q[3];
sx q[3];
rz(-1.4811852) q[3];
sx q[3];
rz(1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(1.3550718) q[0];
rz(2.3782016) q[2];
sx q[2];
rz(-1.1995458) q[2];
sx q[2];
rz(0.064918092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7698344) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(3.1410602) q[1];
x q[2];
rz(0.10947157) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(0.68144875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-0.50450528) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38406661) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(-1.9996044) q[0];
x q[1];
rz(2.675823) q[2];
sx q[2];
rz(-1.2301187) q[2];
sx q[2];
rz(2.6054232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5704755) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(1.8829324) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3584064) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.3822314) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81297368) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(2.5286753) q[0];
rz(-1.8421474) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(3.0999822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61154762) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(0.24286119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7536229) q[3];
sx q[3];
rz(-1.332924) q[3];
sx q[3];
rz(-0.49926234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.3151273) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3430816) q[0];
sx q[0];
rz(-1.8146975) q[0];
sx q[0];
rz(-3.0573465) q[0];
rz(-pi) q[1];
rz(-0.46160134) q[2];
sx q[2];
rz(-1.9399376) q[2];
sx q[2];
rz(-0.75418562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0566237) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(-1.4003537) q[1];
rz(-pi) q[2];
rz(0.58052766) q[3];
sx q[3];
rz(-1.0167828) q[3];
sx q[3];
rz(1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.8483298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9613567) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(0.21813099) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0027577) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(2.5790737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.075899374) q[1];
sx q[1];
rz(-1.7095487) q[1];
sx q[1];
rz(-1.1636415) q[1];
rz(-pi) q[2];
x q[2];
rz(1.771365) q[3];
sx q[3];
rz(-1.9631533) q[3];
sx q[3];
rz(-1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(-1.8539799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33289136) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-2.0602134) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.655517) q[2];
sx q[2];
rz(-1.7239778) q[2];
sx q[2];
rz(-1.4391522) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-2.8833564) q[1];
x q[2];
rz(1.6279531) q[3];
sx q[3];
rz(-2.3944629) q[3];
sx q[3];
rz(-0.77880083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(0.36718711) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8866667) q[0];
sx q[0];
rz(-1.5228132) q[0];
sx q[0];
rz(1.6019437) q[0];
x q[1];
rz(-2.514421) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2314184) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(-0.66399666) q[1];
rz(-0.96079798) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(1.0661351) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-0.74942855) q[2];
sx q[2];
rz(-1.6288169) q[2];
sx q[2];
rz(-1.5734869) q[2];
rz(-1.7442262) q[3];
sx q[3];
rz(-0.72931029) q[3];
sx q[3];
rz(1.9839877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
