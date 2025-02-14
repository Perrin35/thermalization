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
rz(1.2443378) q[0];
sx q[0];
rz(2.3490348) q[0];
sx q[0];
rz(9.6777182) q[0];
rz(-0.77415544) q[1];
sx q[1];
rz(-0.3781265) q[1];
sx q[1];
rz(0.3869431) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9772676) q[0];
sx q[0];
rz(-1.8935793) q[0];
sx q[0];
rz(2.8708145) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5637881) q[2];
sx q[2];
rz(-1.6434323) q[2];
sx q[2];
rz(2.5355123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36648892) q[1];
sx q[1];
rz(-1.6548043) q[1];
sx q[1];
rz(1.4909527) q[1];
x q[2];
rz(2.4693842) q[3];
sx q[3];
rz(-0.70939964) q[3];
sx q[3];
rz(-1.0649452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0953377) q[2];
sx q[2];
rz(-0.75901186) q[2];
sx q[2];
rz(-1.4026027) q[2];
rz(1.7272353) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6801179) q[0];
sx q[0];
rz(-0.48826996) q[0];
sx q[0];
rz(0.71037355) q[0];
rz(2.12517) q[1];
sx q[1];
rz(-1.8417532) q[1];
sx q[1];
rz(-2.3764835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7920096) q[0];
sx q[0];
rz(-0.63711626) q[0];
sx q[0];
rz(2.5883915) q[0];
x q[1];
rz(2.4508453) q[2];
sx q[2];
rz(-1.0260858) q[2];
sx q[2];
rz(-0.70554698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5580235) q[1];
sx q[1];
rz(-2.2127732) q[1];
sx q[1];
rz(-2.7318673) q[1];
x q[2];
rz(-1.160687) q[3];
sx q[3];
rz(-2.7795459) q[3];
sx q[3];
rz(-0.99560478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1171099) q[2];
sx q[2];
rz(-0.64579248) q[2];
sx q[2];
rz(2.4369241) q[2];
rz(-0.17413983) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(-0.38159889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0891377) q[0];
sx q[0];
rz(-1.313504) q[0];
sx q[0];
rz(0.88743368) q[0];
rz(-1.2499836) q[1];
sx q[1];
rz(-1.2108112) q[1];
sx q[1];
rz(-1.3608305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3850573) q[0];
sx q[0];
rz(-2.4891315) q[0];
sx q[0];
rz(2.4908309) q[0];
rz(-pi) q[1];
rz(-3.0615882) q[2];
sx q[2];
rz(-1.7081877) q[2];
sx q[2];
rz(-0.72890711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2584784) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(-1.4264631) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.674747) q[3];
sx q[3];
rz(-0.39109215) q[3];
sx q[3];
rz(2.0147188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7736194) q[2];
sx q[2];
rz(-0.36586389) q[2];
sx q[2];
rz(1.492929) q[2];
rz(-2.6922373) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0538977) q[0];
sx q[0];
rz(-2.5515285) q[0];
sx q[0];
rz(1.7328523) q[0];
rz(-2.159481) q[1];
sx q[1];
rz(-1.5487919) q[1];
sx q[1];
rz(1.7353479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236693) q[0];
sx q[0];
rz(-1.6817998) q[0];
sx q[0];
rz(1.7864321) q[0];
x q[1];
rz(-0.039578118) q[2];
sx q[2];
rz(-1.5180032) q[2];
sx q[2];
rz(-2.9556208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2212053) q[1];
sx q[1];
rz(-1.8169329) q[1];
sx q[1];
rz(-2.1046776) q[1];
x q[2];
rz(-2.8729183) q[3];
sx q[3];
rz(-2.2019535) q[3];
sx q[3];
rz(-2.1674066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1534319) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(2.9856258) q[2];
rz(-0.65034741) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(-0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0680189) q[0];
sx q[0];
rz(-1.8719712) q[0];
sx q[0];
rz(-2.9408348) q[0];
rz(-1.9643895) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(1.1600561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748799) q[0];
sx q[0];
rz(-1.3466517) q[0];
sx q[0];
rz(2.969541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54520901) q[2];
sx q[2];
rz(-3.0070947) q[2];
sx q[2];
rz(-2.4328977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0811942) q[1];
sx q[1];
rz(-1.8738998) q[1];
sx q[1];
rz(-0.36784192) q[1];
x q[2];
rz(-1.0409768) q[3];
sx q[3];
rz(-1.3966832) q[3];
sx q[3];
rz(1.8903738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7572299) q[2];
sx q[2];
rz(-0.94567662) q[2];
sx q[2];
rz(2.6833656) q[2];
rz(2.5107757) q[3];
sx q[3];
rz(-1.9061371) q[3];
sx q[3];
rz(2.6423776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45796564) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(-3.1347347) q[0];
rz(0.231617) q[1];
sx q[1];
rz(-1.7624785) q[1];
sx q[1];
rz(2.0793656) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4214913) q[0];
sx q[0];
rz(-1.4473697) q[0];
sx q[0];
rz(-0.31436679) q[0];
x q[1];
rz(0.49906667) q[2];
sx q[2];
rz(-1.3842513) q[2];
sx q[2];
rz(-1.7507493) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0343218) q[1];
sx q[1];
rz(-2.2244029) q[1];
sx q[1];
rz(1.2493709) q[1];
rz(-1.7627349) q[3];
sx q[3];
rz(-2.4232691) q[3];
sx q[3];
rz(0.11107055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0007533) q[2];
sx q[2];
rz(-1.2522298) q[2];
sx q[2];
rz(1.8737277) q[2];
rz(-0.86152348) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93290257) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(-0.21555899) q[0];
rz(-2.6453099) q[1];
sx q[1];
rz(-1.1880621) q[1];
sx q[1];
rz(2.9821679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27860363) q[0];
sx q[0];
rz(-1.5681603) q[0];
sx q[0];
rz(-2.2056666) q[0];
rz(-1.764545) q[2];
sx q[2];
rz(-2.0341349) q[2];
sx q[2];
rz(3.0510356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9558952) q[1];
sx q[1];
rz(-1.5296827) q[1];
sx q[1];
rz(-1.3354882) q[1];
rz(-pi) q[2];
rz(-2.255329) q[3];
sx q[3];
rz(-2.4733287) q[3];
sx q[3];
rz(-1.9192426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40519199) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(-0.51652017) q[2];
rz(1.9308331) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(-1.3794544) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4844168) q[0];
sx q[0];
rz(-3.0653937) q[0];
sx q[0];
rz(-2.6013689) q[0];
rz(1.6172488) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(-1.2670955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57813533) q[0];
sx q[0];
rz(-1.3523914) q[0];
sx q[0];
rz(2.6294623) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77719633) q[2];
sx q[2];
rz(-2.1146449) q[2];
sx q[2];
rz(-1.0842619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9766751) q[1];
sx q[1];
rz(-1.4265713) q[1];
sx q[1];
rz(0.37101908) q[1];
rz(-pi) q[2];
rz(2.0319875) q[3];
sx q[3];
rz(-2.4465585) q[3];
sx q[3];
rz(2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8871062) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(-0.17733388) q[2];
rz(2.0717428) q[3];
sx q[3];
rz(-1.0515352) q[3];
sx q[3];
rz(2.9593318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3607445) q[0];
sx q[0];
rz(-1.1540664) q[0];
sx q[0];
rz(0.10072197) q[0];
rz(1.992647) q[1];
sx q[1];
rz(-1.1958242) q[1];
sx q[1];
rz(-1.9675868) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3927395) q[0];
sx q[0];
rz(-0.97087384) q[0];
sx q[0];
rz(-1.8003182) q[0];
rz(0.10020013) q[2];
sx q[2];
rz(-2.3728307) q[2];
sx q[2];
rz(2.5454846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8619949) q[1];
sx q[1];
rz(-1.5902385) q[1];
sx q[1];
rz(-1.5330381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2618746) q[3];
sx q[3];
rz(-1.6501325) q[3];
sx q[3];
rz(0.21896958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9775057) q[2];
sx q[2];
rz(-1.4848494) q[2];
sx q[2];
rz(-0.40317765) q[2];
rz(0.66344231) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(3.1373533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73096257) q[0];
sx q[0];
rz(-2.0616489) q[0];
sx q[0];
rz(-3.0626815) q[0];
rz(-0.90947378) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(0.39631072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8412668) q[0];
sx q[0];
rz(-1.5864276) q[0];
sx q[0];
rz(0.0090740694) q[0];
rz(3.134507) q[2];
sx q[2];
rz(-2.3478697) q[2];
sx q[2];
rz(-2.7929579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8974298) q[1];
sx q[1];
rz(-1.8051935) q[1];
sx q[1];
rz(2.6525556) q[1];
rz(-pi) q[2];
x q[2];
rz(2.590606) q[3];
sx q[3];
rz(-2.1175623) q[3];
sx q[3];
rz(-1.8282229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69184715) q[2];
sx q[2];
rz(-0.88374603) q[2];
sx q[2];
rz(-2.3425102) q[2];
rz(0.07829047) q[3];
sx q[3];
rz(-1.8872063) q[3];
sx q[3];
rz(0.6453132) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4770724) q[0];
sx q[0];
rz(-1.7417396) q[0];
sx q[0];
rz(-2.59792) q[0];
rz(1.9480582) q[1];
sx q[1];
rz(-1.2763034) q[1];
sx q[1];
rz(-2.6313849) q[1];
rz(-0.1706201) q[2];
sx q[2];
rz(-1.1191875) q[2];
sx q[2];
rz(2.1714581) q[2];
rz(2.4166475) q[3];
sx q[3];
rz(-1.9564387) q[3];
sx q[3];
rz(2.6993668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
