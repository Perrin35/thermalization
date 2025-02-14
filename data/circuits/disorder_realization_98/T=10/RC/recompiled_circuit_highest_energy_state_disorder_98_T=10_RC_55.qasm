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
rz(-2.383411) q[0];
sx q[0];
rz(-0.36769205) q[0];
sx q[0];
rz(2.0453069) q[0];
rz(-2.3310989) q[1];
sx q[1];
rz(-2.9089622) q[1];
sx q[1];
rz(-1.1059603) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152412) q[0];
sx q[0];
rz(-2.5509074) q[0];
sx q[0];
rz(2.3870941) q[0];
rz(-2.8739086) q[2];
sx q[2];
rz(-1.4776728) q[2];
sx q[2];
rz(-0.47506079) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80820647) q[1];
sx q[1];
rz(-1.5667586) q[1];
sx q[1];
rz(-1.6081078) q[1];
x q[2];
rz(3.1177103) q[3];
sx q[3];
rz(-2.8571354) q[3];
sx q[3];
rz(2.6989355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7652863) q[2];
sx q[2];
rz(-2.6708965) q[2];
sx q[2];
rz(0.37962309) q[2];
rz(-1.5073353) q[3];
sx q[3];
rz(-1.0999271) q[3];
sx q[3];
rz(-0.96989337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955502) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(0.90091339) q[0];
rz(0.13149978) q[1];
sx q[1];
rz(-1.200518) q[1];
sx q[1];
rz(-2.6995755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1832804) q[0];
sx q[0];
rz(-1.4889476) q[0];
sx q[0];
rz(0.32572066) q[0];
rz(0.0046185812) q[2];
sx q[2];
rz(-1.5717662) q[2];
sx q[2];
rz(2.0624954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9513248) q[1];
sx q[1];
rz(-2.3735078) q[1];
sx q[1];
rz(-2.4332341) q[1];
x q[2];
rz(2.7412299) q[3];
sx q[3];
rz(-1.4678848) q[3];
sx q[3];
rz(2.3617552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3758731) q[2];
sx q[2];
rz(-1.7123545) q[2];
sx q[2];
rz(-0.6089375) q[2];
rz(-2.7385312) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(-2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(2.8159863) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(-2.0035279) q[0];
rz(1.552938) q[1];
sx q[1];
rz(-2.373003) q[1];
sx q[1];
rz(0.80702153) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8391221) q[0];
sx q[0];
rz(-0.10286504) q[0];
sx q[0];
rz(-2.8737322) q[0];
rz(-pi) q[1];
rz(1.7337695) q[2];
sx q[2];
rz(-1.1312684) q[2];
sx q[2];
rz(-1.7874315) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.096232) q[1];
sx q[1];
rz(-1.7590471) q[1];
sx q[1];
rz(1.4493311) q[1];
rz(-0.93385796) q[3];
sx q[3];
rz(-0.45972201) q[3];
sx q[3];
rz(-0.62593725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.011586172) q[2];
sx q[2];
rz(-2.5277972) q[2];
sx q[2];
rz(-1.9278795) q[2];
rz(0.064149292) q[3];
sx q[3];
rz(-1.6545273) q[3];
sx q[3];
rz(0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.5478058) q[0];
sx q[0];
rz(-1.8812027) q[0];
sx q[0];
rz(-2.2547145) q[0];
rz(0.15874323) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(1.278272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0905545) q[0];
sx q[0];
rz(-2.0060385) q[0];
sx q[0];
rz(-0.79751517) q[0];
rz(-pi) q[1];
rz(1.099894) q[2];
sx q[2];
rz(-1.0702025) q[2];
sx q[2];
rz(2.4125227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51047303) q[1];
sx q[1];
rz(-0.87952956) q[1];
sx q[1];
rz(3.0635628) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3298755) q[3];
sx q[3];
rz(-0.56801012) q[3];
sx q[3];
rz(-2.3532075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0025803) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(-0.95376897) q[2];
rz(1.9468797) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(-2.6660624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139528) q[0];
sx q[0];
rz(-2.4251745) q[0];
sx q[0];
rz(-2.5415976) q[0];
rz(-2.4586239) q[1];
sx q[1];
rz(-0.78790793) q[1];
sx q[1];
rz(-2.8111828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415057) q[0];
sx q[0];
rz(-1.3887902) q[0];
sx q[0];
rz(-2.299593) q[0];
x q[1];
rz(2.9783713) q[2];
sx q[2];
rz(-1.3484058) q[2];
sx q[2];
rz(1.6878273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40881316) q[1];
sx q[1];
rz(-2.6169852) q[1];
sx q[1];
rz(1.1779768) q[1];
rz(-2.2721259) q[3];
sx q[3];
rz(-1.2463797) q[3];
sx q[3];
rz(-0.78044915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43634513) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(-2.516563) q[2];
rz(0.69508067) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(3.0888016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92796749) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(1.3731765) q[0];
rz(-1.7534509) q[1];
sx q[1];
rz(-0.87166798) q[1];
sx q[1];
rz(1.53481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4563303) q[0];
sx q[0];
rz(-1.7182516) q[0];
sx q[0];
rz(3.1059424) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8839726) q[2];
sx q[2];
rz(-1.3968236) q[2];
sx q[2];
rz(-0.18944511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36027137) q[1];
sx q[1];
rz(-0.92642537) q[1];
sx q[1];
rz(2.086198) q[1];
x q[2];
rz(1.9717384) q[3];
sx q[3];
rz(-1.6982659) q[3];
sx q[3];
rz(-1.0893703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(1.6004174) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(2.8405564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.130126) q[0];
sx q[0];
rz(-3.0045894) q[0];
sx q[0];
rz(0.89114183) q[0];
rz(-0.64542422) q[1];
sx q[1];
rz(-1.7452469) q[1];
sx q[1];
rz(-0.83271629) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7695205) q[0];
sx q[0];
rz(-1.4345048) q[0];
sx q[0];
rz(-0.46532495) q[0];
x q[1];
rz(-3.0751257) q[2];
sx q[2];
rz(-0.94777135) q[2];
sx q[2];
rz(-2.6084171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5834076) q[1];
sx q[1];
rz(-1.7654357) q[1];
sx q[1];
rz(2.6411112) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5359573) q[3];
sx q[3];
rz(-0.6303595) q[3];
sx q[3];
rz(1.039884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58925313) q[2];
sx q[2];
rz(-2.4746042) q[2];
sx q[2];
rz(-1.5480631) q[2];
rz(0.42406905) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(2.8467395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9645204) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(2.3181584) q[0];
rz(-1.9442762) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(1.6808602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6525424) q[0];
sx q[0];
rz(-0.33777896) q[0];
sx q[0];
rz(-1.9073639) q[0];
rz(-pi) q[1];
rz(0.97779556) q[2];
sx q[2];
rz(-1.4306746) q[2];
sx q[2];
rz(1.8210379) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6951235) q[1];
sx q[1];
rz(-2.2317076) q[1];
sx q[1];
rz(1.5170694) q[1];
rz(1.8917055) q[3];
sx q[3];
rz(-1.9139288) q[3];
sx q[3];
rz(0.58716256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0494277) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(3.0180422) q[2];
rz(0.83856797) q[3];
sx q[3];
rz(-1.1127915) q[3];
sx q[3];
rz(0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0274444) q[0];
sx q[0];
rz(-1.0076948) q[0];
sx q[0];
rz(1.7281519) q[0];
rz(0.89883262) q[1];
sx q[1];
rz(-0.90679589) q[1];
sx q[1];
rz(-2.5487505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44594793) q[0];
sx q[0];
rz(-0.88291016) q[0];
sx q[0];
rz(-0.77151973) q[0];
x q[1];
rz(2.9594764) q[2];
sx q[2];
rz(-1.4859293) q[2];
sx q[2];
rz(0.1402459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.034380091) q[1];
sx q[1];
rz(-1.980956) q[1];
sx q[1];
rz(-0.35663794) q[1];
rz(-pi) q[2];
rz(0.70609181) q[3];
sx q[3];
rz(-0.81940813) q[3];
sx q[3];
rz(-3.135526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8574519) q[2];
sx q[2];
rz(-0.74791869) q[2];
sx q[2];
rz(-2.9202374) q[2];
rz(1.0434693) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(-0.38844696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2986044) q[0];
sx q[0];
rz(-2.8589111) q[0];
sx q[0];
rz(-0.84841949) q[0];
rz(0.09659718) q[1];
sx q[1];
rz(-2.0674457) q[1];
sx q[1];
rz(-2.3883147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71614186) q[0];
sx q[0];
rz(-1.6707067) q[0];
sx q[0];
rz(-1.3425171) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6069372) q[2];
sx q[2];
rz(-0.74191626) q[2];
sx q[2];
rz(-0.20797543) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92434525) q[1];
sx q[1];
rz(-2.0567523) q[1];
sx q[1];
rz(0.43083453) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6274458) q[3];
sx q[3];
rz(-0.78350583) q[3];
sx q[3];
rz(1.0422106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0640556) q[2];
sx q[2];
rz(-2.3040743) q[2];
sx q[2];
rz(2.688664) q[2];
rz(0.60106599) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-2.1452904) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30408981) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(2.9651463) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(-2.7929581) q[2];
sx q[2];
rz(-1.5713816) q[2];
sx q[2];
rz(-1.1272507) q[2];
rz(-2.4963958) q[3];
sx q[3];
rz(-2.3020036) q[3];
sx q[3];
rz(0.19097542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
