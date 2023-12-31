OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(6.531125) q[0];
sx q[0];
rz(8.6046435) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983842) q[0];
sx q[0];
rz(-2.0294667) q[0];
sx q[0];
rz(-2.0660603) q[0];
x q[1];
rz(-1.0432265) q[2];
sx q[2];
rz(-1.9448148) q[2];
sx q[2];
rz(-2.5803215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80884113) q[1];
sx q[1];
rz(-1.0665227) q[1];
sx q[1];
rz(-1.5682975) q[1];
x q[2];
rz(-3.1357364) q[3];
sx q[3];
rz(-1.3121038) q[3];
sx q[3];
rz(2.4394025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(-1.7830085) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434718) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(-1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.9083317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-0.92139771) q[0];
sx q[0];
rz(0.68769023) q[0];
x q[1];
rz(-1.6987259) q[2];
sx q[2];
rz(-0.90020858) q[2];
sx q[2];
rz(2.6562481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.997628) q[1];
sx q[1];
rz(-1.3647172) q[1];
sx q[1];
rz(0.016184316) q[1];
rz(-0.66744653) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(-0.28039704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(1.8921651) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1059589) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(-1.3282302) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0257646) q[2];
sx q[2];
rz(-1.8075382) q[2];
sx q[2];
rz(0.67982212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0752807) q[1];
sx q[1];
rz(-1.7048786) q[1];
sx q[1];
rz(0.65384298) q[1];
rz(-pi) q[2];
rz(-2.4124868) q[3];
sx q[3];
rz(-2.2359072) q[3];
sx q[3];
rz(0.8641181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.694214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3576413) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(0.95143239) q[0];
rz(3.1324603) q[2];
sx q[2];
rz(-1.3924686) q[2];
sx q[2];
rz(0.4724617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4923258) q[1];
sx q[1];
rz(-2.4267303) q[1];
sx q[1];
rz(1.3267172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57595423) q[3];
sx q[3];
rz(-0.42612694) q[3];
sx q[3];
rz(-1.061561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(2.9966667) q[2];
rz(2.1302917) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(0.68403912) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.2129983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1367462) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(0.61704163) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2239387) q[2];
sx q[2];
rz(-1.0461763) q[2];
sx q[2];
rz(1.7078924) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8634833) q[1];
sx q[1];
rz(-0.28243318) q[1];
sx q[1];
rz(-2.7062098) q[1];
rz(-pi) q[2];
rz(-2.9641987) q[3];
sx q[3];
rz(-1.8284441) q[3];
sx q[3];
rz(0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1155788) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-0.58475959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023914) q[0];
sx q[0];
rz(-2.4524134) q[0];
sx q[0];
rz(-0.50424772) q[0];
rz(1.5445968) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(0.40528497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6962291) q[1];
sx q[1];
rz(-2.3949957) q[1];
sx q[1];
rz(-1.1397051) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8486409) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(1.1374744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.6284846) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(-2.5351977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8618146) q[0];
sx q[0];
rz(-1.8068411) q[0];
sx q[0];
rz(1.0311014) q[0];
rz(-pi) q[1];
rz(-2.0386001) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(1.4008092) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7402621) q[1];
sx q[1];
rz(-1.4583424) q[1];
sx q[1];
rz(-2.3349891) q[1];
x q[2];
rz(-0.86319478) q[3];
sx q[3];
rz(-0.51897012) q[3];
sx q[3];
rz(-2.2591345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(2.9139013) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75922155) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(1.9842149) q[0];
rz(1.20485) q[2];
sx q[2];
rz(-1.1507251) q[2];
sx q[2];
rz(-0.36047381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3701054) q[1];
sx q[1];
rz(-0.58413726) q[1];
sx q[1];
rz(0.41569709) q[1];
rz(1.6119192) q[3];
sx q[3];
rz(-2.1293981) q[3];
sx q[3];
rz(-2.5839644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(0.84028876) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(0.80378419) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(1.1484336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0949806) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(-1.6427965) q[0];
rz(-pi) q[1];
rz(-1.0461651) q[2];
sx q[2];
rz(-1.8881919) q[2];
sx q[2];
rz(1.5898926) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0828404) q[1];
sx q[1];
rz(-1.0995004) q[1];
sx q[1];
rz(-2.6575762) q[1];
rz(-2.7896342) q[3];
sx q[3];
rz(-1.6037914) q[3];
sx q[3];
rz(1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(0.69828066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822599) q[0];
sx q[0];
rz(-1.8730622) q[0];
sx q[0];
rz(-1.6460101) q[0];
rz(-1.6449498) q[2];
sx q[2];
rz(-1.5955131) q[2];
sx q[2];
rz(-2.2029684) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0638949) q[1];
sx q[1];
rz(-2.6297036) q[1];
sx q[1];
rz(-1.1179395) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9395589) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-0.84021604) q[2];
rz(0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-1.2226979) q[2];
sx q[2];
rz(-2.1445027) q[2];
sx q[2];
rz(2.4377844) q[2];
rz(1.7520262) q[3];
sx q[3];
rz(-0.69835865) q[3];
sx q[3];
rz(-2.9327304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
