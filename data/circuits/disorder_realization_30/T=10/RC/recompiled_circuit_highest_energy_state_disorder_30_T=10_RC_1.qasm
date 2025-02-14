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
rz(0.93267814) q[0];
sx q[0];
rz(-2.7490766) q[0];
sx q[0];
rz(2.0445332) q[0];
rz(1.2708083) q[1];
sx q[1];
rz(-1.1975809) q[1];
sx q[1];
rz(1.4374011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5839489) q[0];
sx q[0];
rz(-0.93867477) q[0];
sx q[0];
rz(1.1313788) q[0];
rz(-pi) q[1];
rz(2.9095638) q[2];
sx q[2];
rz(-1.4995607) q[2];
sx q[2];
rz(0.89398065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37482191) q[1];
sx q[1];
rz(-1.8954191) q[1];
sx q[1];
rz(0.9379106) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1278856) q[3];
sx q[3];
rz(-0.64798149) q[3];
sx q[3];
rz(1.614384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8561594) q[2];
sx q[2];
rz(-2.9069052) q[2];
sx q[2];
rz(2.6402546) q[2];
rz(-1.8173789) q[3];
sx q[3];
rz(-0.54558498) q[3];
sx q[3];
rz(1.4668303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.685574) q[0];
sx q[0];
rz(-2.7360003) q[0];
sx q[0];
rz(1.4612041) q[0];
rz(-2.9623518) q[1];
sx q[1];
rz(-0.76911175) q[1];
sx q[1];
rz(0.63757149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3831275) q[0];
sx q[0];
rz(-0.33184856) q[0];
sx q[0];
rz(-0.76803523) q[0];
rz(-1.7098956) q[2];
sx q[2];
rz(-0.31965986) q[2];
sx q[2];
rz(-2.5676198) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3481787) q[1];
sx q[1];
rz(-1.2169588) q[1];
sx q[1];
rz(-2.2458312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9470975) q[3];
sx q[3];
rz(-2.199389) q[3];
sx q[3];
rz(-0.95910536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3011938) q[2];
sx q[2];
rz(-2.0638778) q[2];
sx q[2];
rz(3.0974498) q[2];
rz(0.61331493) q[3];
sx q[3];
rz(-1.3215348) q[3];
sx q[3];
rz(0.80673748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089207) q[0];
sx q[0];
rz(-0.044600211) q[0];
sx q[0];
rz(-1.3432304) q[0];
rz(-1.4187468) q[1];
sx q[1];
rz(-0.90444618) q[1];
sx q[1];
rz(0.80352965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4889824) q[0];
sx q[0];
rz(-1.6296415) q[0];
sx q[0];
rz(0.68025689) q[0];
rz(-0.5838238) q[2];
sx q[2];
rz(-1.1808504) q[2];
sx q[2];
rz(2.4175298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60119286) q[1];
sx q[1];
rz(-1.4638888) q[1];
sx q[1];
rz(1.8501297) q[1];
rz(-0.66763287) q[3];
sx q[3];
rz(-0.61452319) q[3];
sx q[3];
rz(0.4947084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58328763) q[2];
sx q[2];
rz(-0.95197695) q[2];
sx q[2];
rz(2.1236911) q[2];
rz(-1.2244276) q[3];
sx q[3];
rz(-1.8095576) q[3];
sx q[3];
rz(2.0126191) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99854904) q[0];
sx q[0];
rz(-0.7989378) q[0];
sx q[0];
rz(2.1813188) q[0];
rz(-0.1943365) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(-0.0050553102) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968648) q[0];
sx q[0];
rz(-0.93679777) q[0];
sx q[0];
rz(0.89374884) q[0];
rz(-pi) q[1];
rz(-1.5326535) q[2];
sx q[2];
rz(-1.3873709) q[2];
sx q[2];
rz(-2.9596552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3761085) q[1];
sx q[1];
rz(-1.6836201) q[1];
sx q[1];
rz(2.9214431) q[1];
rz(-pi) q[2];
x q[2];
rz(1.051733) q[3];
sx q[3];
rz(-1.6994275) q[3];
sx q[3];
rz(-0.41384709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.888835) q[2];
sx q[2];
rz(-1.5450666) q[2];
sx q[2];
rz(1.4463536) q[2];
rz(1.3085922) q[3];
sx q[3];
rz(-1.3819709) q[3];
sx q[3];
rz(2.5368209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62110353) q[0];
sx q[0];
rz(-0.72431722) q[0];
sx q[0];
rz(-2.5526168) q[0];
rz(-0.0079872459) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(2.0569107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0953012) q[0];
sx q[0];
rz(-1.8043552) q[0];
sx q[0];
rz(0.26136847) q[0];
x q[1];
rz(-3.1118664) q[2];
sx q[2];
rz(-0.27966248) q[2];
sx q[2];
rz(3.0936808) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1639203) q[1];
sx q[1];
rz(-2.2798988) q[1];
sx q[1];
rz(0.52859882) q[1];
rz(-pi) q[2];
rz(1.6570377) q[3];
sx q[3];
rz(-0.6644333) q[3];
sx q[3];
rz(2.2905615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7320554) q[2];
sx q[2];
rz(-1.0169225) q[2];
sx q[2];
rz(-1.766905) q[2];
rz(-0.095666766) q[3];
sx q[3];
rz(-1.5868264) q[3];
sx q[3];
rz(-2.5110551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(-1.3404982) q[0];
rz(-0.66572491) q[1];
sx q[1];
rz(-1.3434854) q[1];
sx q[1];
rz(-0.7241157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511696) q[0];
sx q[0];
rz(-1.2688338) q[0];
sx q[0];
rz(-1.9504551) q[0];
rz(-pi) q[1];
rz(-2.3409136) q[2];
sx q[2];
rz(-0.73298798) q[2];
sx q[2];
rz(-1.4633601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.648702) q[1];
sx q[1];
rz(-0.48149432) q[1];
sx q[1];
rz(0.29701155) q[1];
rz(1.6203262) q[3];
sx q[3];
rz(-0.096442744) q[3];
sx q[3];
rz(0.24094757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93806997) q[2];
sx q[2];
rz(-1.4504434) q[2];
sx q[2];
rz(-0.098048992) q[2];
rz(0.34052643) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(-0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53511867) q[0];
sx q[0];
rz(-0.70347324) q[0];
sx q[0];
rz(0.055543609) q[0];
rz(0.37560383) q[1];
sx q[1];
rz(-2.3665078) q[1];
sx q[1];
rz(1.4449878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9268413) q[0];
sx q[0];
rz(-2.3390649) q[0];
sx q[0];
rz(1.7178128) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3705424) q[2];
sx q[2];
rz(-2.8828388) q[2];
sx q[2];
rz(3.0508397) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0241131) q[1];
sx q[1];
rz(-1.584519) q[1];
sx q[1];
rz(-2.7639469) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4395326) q[3];
sx q[3];
rz(-0.5834879) q[3];
sx q[3];
rz(1.7398011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.043341788) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(-3.1116389) q[3];
sx q[3];
rz(-1.5104537) q[3];
sx q[3];
rz(-2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444262) q[0];
sx q[0];
rz(-2.501896) q[0];
sx q[0];
rz(-1.8901012) q[0];
rz(-2.331612) q[1];
sx q[1];
rz(-0.73695838) q[1];
sx q[1];
rz(1.6862148) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90235119) q[0];
sx q[0];
rz(-2.3555814) q[0];
sx q[0];
rz(2.8186574) q[0];
rz(1.6258442) q[2];
sx q[2];
rz(-2.0924797) q[2];
sx q[2];
rz(1.59984) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7500893) q[1];
sx q[1];
rz(-2.1632101) q[1];
sx q[1];
rz(-0.77880145) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2627566) q[3];
sx q[3];
rz(-2.2744826) q[3];
sx q[3];
rz(-1.7233985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84888419) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.1513618) q[2];
rz(0.014160841) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8779811) q[0];
sx q[0];
rz(-2.4284555) q[0];
sx q[0];
rz(-2.1942595) q[0];
rz(0.5961279) q[1];
sx q[1];
rz(-1.7201741) q[1];
sx q[1];
rz(-1.5790342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2939908) q[0];
sx q[0];
rz(-2.2819073) q[0];
sx q[0];
rz(1.343665) q[0];
x q[1];
rz(-2.8364592) q[2];
sx q[2];
rz(-2.6465979) q[2];
sx q[2];
rz(-2.309805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4217675) q[1];
sx q[1];
rz(-0.9523069) q[1];
sx q[1];
rz(2.5773125) q[1];
rz(2.4189255) q[3];
sx q[3];
rz(-2.8407466) q[3];
sx q[3];
rz(-2.2900555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24937135) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(0.48386827) q[2];
rz(-0.59746915) q[3];
sx q[3];
rz(-1.4474844) q[3];
sx q[3];
rz(-0.57815236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.4018965) q[0];
sx q[0];
rz(-1.5902436) q[0];
sx q[0];
rz(2.7751112) q[0];
rz(1.81555) q[1];
sx q[1];
rz(-0.96249023) q[1];
sx q[1];
rz(-2.0749626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410927) q[0];
sx q[0];
rz(-2.0632732) q[0];
sx q[0];
rz(-0.89333138) q[0];
rz(-0.60619715) q[2];
sx q[2];
rz(-1.6087039) q[2];
sx q[2];
rz(-2.4647922) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45056709) q[1];
sx q[1];
rz(-0.98916173) q[1];
sx q[1];
rz(-2.8514991) q[1];
x q[2];
rz(1.0935547) q[3];
sx q[3];
rz(-1.541409) q[3];
sx q[3];
rz(-0.060197006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8804973) q[2];
sx q[2];
rz(-2.3636621) q[2];
sx q[2];
rz(-2.0757389) q[2];
rz(2.8737658) q[3];
sx q[3];
rz(-1.6949751) q[3];
sx q[3];
rz(0.48521391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601226) q[0];
sx q[0];
rz(-2.1385834) q[0];
sx q[0];
rz(-2.9271097) q[0];
rz(-2.8196234) q[1];
sx q[1];
rz(-1.3170769) q[1];
sx q[1];
rz(-1.9016686) q[1];
rz(1.5628846) q[2];
sx q[2];
rz(-1.2894165) q[2];
sx q[2];
rz(-3.0647562) q[2];
rz(0.81372502) q[3];
sx q[3];
rz(-0.56725262) q[3];
sx q[3];
rz(-0.57813416) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
