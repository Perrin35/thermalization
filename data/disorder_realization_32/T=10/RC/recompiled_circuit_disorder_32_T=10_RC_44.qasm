OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(0.16103345) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8344603) q[0];
sx q[0];
rz(-1.6856598) q[0];
sx q[0];
rz(0.97712028) q[0];
x q[1];
rz(-0.52486323) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(-0.85688574) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7339488) q[3];
sx q[3];
rz(-1.6626433) q[3];
sx q[3];
rz(-0.44029217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(2.0387409) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(0.57759181) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.577852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.455362) q[0];
x q[1];
rz(-0.30303843) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(-2.9532202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6985804) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.3742616) q[1];
rz(-pi) q[2];
rz(2.9868449) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(1.9783431) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(-2.5015586) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.144369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(0.46023603) q[0];
x q[1];
rz(-2.6687713) q[2];
sx q[2];
rz(-2.678215) q[2];
sx q[2];
rz(0.18714999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80013393) q[1];
sx q[1];
rz(-1.1143436) q[1];
sx q[1];
rz(-1.1851428) q[1];
rz(-pi) q[2];
rz(-0.010300962) q[3];
sx q[3];
rz(-0.28763887) q[3];
sx q[3];
rz(2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-2.3864746) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0621322) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(0.53702766) q[0];
x q[1];
rz(0.56447345) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(2.0828473) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48651931) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(-1.4763114) q[1];
rz(-pi) q[2];
rz(-2.643814) q[3];
sx q[3];
rz(-1.0504424) q[3];
sx q[3];
rz(-1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.283678) q[0];
sx q[0];
rz(0.65335269) q[0];
rz(-pi) q[1];
rz(2.2064477) q[2];
sx q[2];
rz(-1.9872553) q[2];
sx q[2];
rz(-0.93174975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.000948) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(2.7468202) q[1];
rz(-0.3211) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-0.54164106) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(2.0419962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8894419) q[0];
sx q[0];
rz(-0.60177207) q[0];
sx q[0];
rz(-2.2579231) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8855368) q[2];
sx q[2];
rz(-2.5883) q[2];
sx q[2];
rz(-1.3421343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53193608) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(-2.460536) q[1];
rz(-2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-0.15730102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0387715) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(-3.0022013) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.332875) q[2];
sx q[2];
rz(-2.3908797) q[2];
sx q[2];
rz(2.6964292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0485222) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(2.1919769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7492206) q[3];
sx q[3];
rz(-2.4761768) q[3];
sx q[3];
rz(2.3371646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5787443) q[0];
sx q[0];
rz(-1.0283854) q[0];
sx q[0];
rz(-1.5638652) q[0];
rz(-pi) q[1];
rz(2.2145055) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(1.5333652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5756702) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-0.5528321) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1612438) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(2.6524515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(-0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(1.3508505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8838053) q[0];
sx q[0];
rz(-1.4865849) q[0];
sx q[0];
rz(-3.0620831) q[0];
rz(-pi) q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-2.2293408) q[2];
sx q[2];
rz(1.272162) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46637725) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(1.1364163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3488995) q[3];
sx q[3];
rz(-1.6716692) q[3];
sx q[3];
rz(1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1039625) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(-2.2274341) q[0];
rz(-1.5461966) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(1.8116902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5369397) q[1];
sx q[1];
rz(-2.5250658) q[1];
sx q[1];
rz(1.489868) q[1];
rz(-1.8923558) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-2.3615169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(2.8725273) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(0.11733304) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
