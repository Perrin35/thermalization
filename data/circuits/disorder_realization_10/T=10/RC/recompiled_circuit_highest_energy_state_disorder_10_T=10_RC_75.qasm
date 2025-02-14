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
rz(0.34974521) q[0];
sx q[0];
rz(5.3164696) q[0];
sx q[0];
rz(9.7860019) q[0];
rz(-0.53520441) q[1];
sx q[1];
rz(1.9898131) q[1];
sx q[1];
rz(10.726816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068548) q[0];
sx q[0];
rz(-1.2465205) q[0];
sx q[0];
rz(1.2756596) q[0];
x q[1];
rz(-0.041000276) q[2];
sx q[2];
rz(-0.71913856) q[2];
sx q[2];
rz(0.96895167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.227461) q[1];
sx q[1];
rz(-1.1644851) q[1];
sx q[1];
rz(1.8309469) q[1];
rz(2.8252772) q[3];
sx q[3];
rz(-0.22238734) q[3];
sx q[3];
rz(-2.5393644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5754622) q[2];
sx q[2];
rz(-0.29511109) q[2];
sx q[2];
rz(2.5556514) q[2];
rz(1.640004) q[3];
sx q[3];
rz(-1.274704) q[3];
sx q[3];
rz(1.2007825) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457299) q[0];
sx q[0];
rz(-2.0450617) q[0];
sx q[0];
rz(-2.0134266) q[0];
rz(0.7496756) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(-1.658176) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2735398) q[0];
sx q[0];
rz(-0.85208396) q[0];
sx q[0];
rz(2.5107288) q[0];
rz(0.60083484) q[2];
sx q[2];
rz(-1.2640096) q[2];
sx q[2];
rz(0.50150774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0771328) q[1];
sx q[1];
rz(-1.049548) q[1];
sx q[1];
rz(1.8690994) q[1];
rz(0.97256487) q[3];
sx q[3];
rz(-0.56275193) q[3];
sx q[3];
rz(-2.5431354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38773203) q[2];
sx q[2];
rz(-0.78247491) q[2];
sx q[2];
rz(-2.1369047) q[2];
rz(0.0082155148) q[3];
sx q[3];
rz(-1.3722082) q[3];
sx q[3];
rz(2.7042702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.8179271) q[0];
sx q[0];
rz(-1.265047) q[0];
sx q[0];
rz(1.0311968) q[0];
rz(-2.9464856) q[1];
sx q[1];
rz(-2.52774) q[1];
sx q[1];
rz(-0.88417792) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3613286) q[0];
sx q[0];
rz(-1.5710281) q[0];
sx q[0];
rz(0.0032629691) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6818037) q[2];
sx q[2];
rz(-1.3693235) q[2];
sx q[2];
rz(-0.10290621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14890524) q[1];
sx q[1];
rz(-0.42154781) q[1];
sx q[1];
rz(-1.804638) q[1];
rz(0.98434358) q[3];
sx q[3];
rz(-1.801681) q[3];
sx q[3];
rz(1.6572052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6241793) q[2];
sx q[2];
rz(-1.2531345) q[2];
sx q[2];
rz(0.45428983) q[2];
rz(-2.1381461) q[3];
sx q[3];
rz(-2.5277621) q[3];
sx q[3];
rz(-1.412089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91785947) q[0];
sx q[0];
rz(-2.7847325) q[0];
sx q[0];
rz(-2.3339363) q[0];
rz(-1.3937021) q[1];
sx q[1];
rz(-1.7180157) q[1];
sx q[1];
rz(-1.4180988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7515441) q[0];
sx q[0];
rz(-1.1466197) q[0];
sx q[0];
rz(2.9817355) q[0];
rz(-pi) q[1];
rz(0.22252616) q[2];
sx q[2];
rz(-2.4234802) q[2];
sx q[2];
rz(-1.1873174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0102301) q[1];
sx q[1];
rz(-1.5585962) q[1];
sx q[1];
rz(0.38174816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0491761) q[3];
sx q[3];
rz(-1.9289013) q[3];
sx q[3];
rz(1.988609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35004804) q[2];
sx q[2];
rz(-1.9106661) q[2];
sx q[2];
rz(3.0080504) q[2];
rz(0.40889016) q[3];
sx q[3];
rz(-3.0372527) q[3];
sx q[3];
rz(1.0094118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.908476) q[0];
sx q[0];
rz(-2.3148843) q[0];
sx q[0];
rz(0.46434656) q[0];
rz(0.50114477) q[1];
sx q[1];
rz(-2.8906288) q[1];
sx q[1];
rz(-0.73890013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66128816) q[0];
sx q[0];
rz(-1.4562074) q[0];
sx q[0];
rz(0.25729723) q[0];
rz(-pi) q[1];
rz(0.22245714) q[2];
sx q[2];
rz(-2.5875666) q[2];
sx q[2];
rz(1.2740327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0245783) q[1];
sx q[1];
rz(-2.0110333) q[1];
sx q[1];
rz(-1.8135032) q[1];
rz(-2.49315) q[3];
sx q[3];
rz(-1.1243226) q[3];
sx q[3];
rz(-0.97383271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.240694) q[2];
sx q[2];
rz(-1.4229166) q[2];
sx q[2];
rz(2.9071729) q[2];
rz(-2.8068986) q[3];
sx q[3];
rz(-0.60939279) q[3];
sx q[3];
rz(1.7083683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1471106) q[0];
sx q[0];
rz(-2.2587977) q[0];
sx q[0];
rz(2.2301646) q[0];
rz(1.8270739) q[1];
sx q[1];
rz(-0.85845566) q[1];
sx q[1];
rz(1.3988769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0056038) q[0];
sx q[0];
rz(-2.4244747) q[0];
sx q[0];
rz(0.75053458) q[0];
x q[1];
rz(-1.5406467) q[2];
sx q[2];
rz(-1.4316403) q[2];
sx q[2];
rz(2.2517532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28820747) q[1];
sx q[1];
rz(-2.7010272) q[1];
sx q[1];
rz(-1.8024995) q[1];
rz(-2.1067922) q[3];
sx q[3];
rz(-1.7977996) q[3];
sx q[3];
rz(1.5569527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92480245) q[2];
sx q[2];
rz(-0.81764597) q[2];
sx q[2];
rz(-1.372288) q[2];
rz(-1.4745845) q[3];
sx q[3];
rz(-1.0627012) q[3];
sx q[3];
rz(0.32349989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.737616) q[0];
sx q[0];
rz(-2.1196899) q[0];
sx q[0];
rz(1.2603941) q[0];
rz(0.34126392) q[1];
sx q[1];
rz(-0.9287467) q[1];
sx q[1];
rz(2.0492679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3288984) q[0];
sx q[0];
rz(-0.51691953) q[0];
sx q[0];
rz(-3.10269) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7997362) q[2];
sx q[2];
rz(-1.3876788) q[2];
sx q[2];
rz(-1.7464856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1398692) q[1];
sx q[1];
rz(-2.4635063) q[1];
sx q[1];
rz(-1.1820611) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.722205) q[3];
sx q[3];
rz(-1.2817146) q[3];
sx q[3];
rz(1.4037286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1616538) q[2];
sx q[2];
rz(-1.7540163) q[2];
sx q[2];
rz(-2.2765344) q[2];
rz(3.0480399) q[3];
sx q[3];
rz(-1.4253987) q[3];
sx q[3];
rz(3.0430999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.426067) q[0];
sx q[0];
rz(-2.1676097) q[0];
sx q[0];
rz(-1.4564212) q[0];
rz(2.9086225) q[1];
sx q[1];
rz(-1.7590245) q[1];
sx q[1];
rz(1.9219386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7160952) q[0];
sx q[0];
rz(-0.78491917) q[0];
sx q[0];
rz(-0.4161448) q[0];
x q[1];
rz(0.43841534) q[2];
sx q[2];
rz(-1.9389956) q[2];
sx q[2];
rz(0.52939289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8792272) q[1];
sx q[1];
rz(-2.7272537) q[1];
sx q[1];
rz(0.44466876) q[1];
x q[2];
rz(0.69583054) q[3];
sx q[3];
rz(-1.1828109) q[3];
sx q[3];
rz(1.4881021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.098027078) q[2];
sx q[2];
rz(-1.7003912) q[2];
sx q[2];
rz(-3.0230057) q[2];
rz(0.81973997) q[3];
sx q[3];
rz(-2.6974862) q[3];
sx q[3];
rz(2.8290101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176625) q[0];
sx q[0];
rz(-0.95159641) q[0];
sx q[0];
rz(-0.96624017) q[0];
rz(0.36695925) q[1];
sx q[1];
rz(-2.1789813) q[1];
sx q[1];
rz(0.56698322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2174187) q[0];
sx q[0];
rz(-1.770562) q[0];
sx q[0];
rz(-2.6471247) q[0];
x q[1];
rz(2.1927715) q[2];
sx q[2];
rz(-0.89911956) q[2];
sx q[2];
rz(-0.27510422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3388515) q[1];
sx q[1];
rz(-0.88715645) q[1];
sx q[1];
rz(1.834447) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8993239) q[3];
sx q[3];
rz(-1.2905811) q[3];
sx q[3];
rz(-1.3485419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6734267) q[2];
sx q[2];
rz(-1.4665073) q[2];
sx q[2];
rz(1.0183938) q[2];
rz(2.9112725) q[3];
sx q[3];
rz(-0.46263328) q[3];
sx q[3];
rz(3.0751626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224943) q[0];
sx q[0];
rz(-1.007217) q[0];
sx q[0];
rz(-1.3315211) q[0];
rz(1.9271756) q[1];
sx q[1];
rz(-1.6969095) q[1];
sx q[1];
rz(-2.725098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0016704069) q[0];
sx q[0];
rz(-0.71437144) q[0];
sx q[0];
rz(-1.2561965) q[0];
rz(-pi) q[1];
rz(1.824786) q[2];
sx q[2];
rz(-1.2653627) q[2];
sx q[2];
rz(2.0485725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0076722) q[1];
sx q[1];
rz(-0.67613542) q[1];
sx q[1];
rz(2.424404) q[1];
rz(-2.0216398) q[3];
sx q[3];
rz(-1.6623673) q[3];
sx q[3];
rz(-1.4705603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3930964) q[2];
sx q[2];
rz(-1.2747719) q[2];
sx q[2];
rz(1.3416802) q[2];
rz(-2.0210576) q[3];
sx q[3];
rz(-1.7399961) q[3];
sx q[3];
rz(0.11693461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31711598) q[0];
sx q[0];
rz(-1.8335637) q[0];
sx q[0];
rz(1.2807922) q[0];
rz(-0.9723797) q[1];
sx q[1];
rz(-2.0449816) q[1];
sx q[1];
rz(2.4349946) q[1];
rz(2.7823051) q[2];
sx q[2];
rz(-2.4854598) q[2];
sx q[2];
rz(1.796915) q[2];
rz(1.9119143) q[3];
sx q[3];
rz(-0.096466168) q[3];
sx q[3];
rz(-1.1763541) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
