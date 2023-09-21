OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.035064) q[0];
sx q[0];
rz(-2.0523235) q[0];
sx q[0];
rz(2.9805592) q[0];
rz(-1.5313907) q[1];
sx q[1];
rz(-2.664497) q[1];
sx q[1];
rz(2.6452126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3408605) q[0];
sx q[0];
rz(-2.1600318) q[0];
sx q[0];
rz(-3.0032934) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8018087) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(-1.9356188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.032174234) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.4515619) q[1];
rz(-0.22827893) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(-1.3397863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751511) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(-0.22002797) q[0];
x q[1];
rz(-2.4096476) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(-1.0674455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4430122) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.3742616) q[1];
rz(-2.3489477) q[3];
sx q[3];
rz(-0.21867293) q[3];
sx q[3];
rz(1.45989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(-2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0796233) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(-0.18297718) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9336752) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(-1.0008706) q[0];
x q[1];
rz(0.41855721) q[2];
sx q[2];
rz(-1.7757799) q[2];
sx q[2];
rz(2.1870854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80013393) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(-1.1851428) q[1];
x q[2];
rz(-2.8539682) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(-0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1220876) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6782137) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(2.148927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216114) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(2.6556334) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2797212) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(-2.9574403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0419473) q[1];
sx q[1];
rz(-1.4775659) q[1];
sx q[1];
rz(2.9779742) q[1];
rz(0.87611115) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(2.0640404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(2.7899182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36067097) q[0];
sx q[0];
rz(-2.43649) q[0];
sx q[0];
rz(2.6893925) q[0];
rz(0.93047662) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(3.004068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.56375757) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(-1.5866336) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8204927) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(2.0146446) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(1.0995964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(0.41082541) q[0];
x q[1];
rz(-2.1018283) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(2.642717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82367831) q[1];
sx q[1];
rz(-0.70033973) q[1];
sx q[1];
rz(-2.8631696) q[1];
rz(-pi) q[2];
rz(-0.78824708) q[3];
sx q[3];
rz(-1.8133834) q[3];
sx q[3];
rz(1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(2.8586094) q[2];
rz(-0.7061559) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(-0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.776942) q[0];
sx q[0];
rz(-1.4770664) q[0];
sx q[0];
rz(-2.4073497) q[0];
rz(0.21643164) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(0.12491465) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3992608) q[1];
sx q[1];
rz(-0.6472339) q[1];
sx q[1];
rz(-1.8982235) q[1];
rz(1.2792475) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(2.8214422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1823696) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(0.29327926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54942185) q[0];
sx q[0];
rz(-0.54245078) q[0];
sx q[0];
rz(0.011499238) q[0];
rz(-pi) q[1];
rz(0.92708712) q[2];
sx q[2];
rz(-1.7037399) q[2];
sx q[2];
rz(1.5333652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5659225) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-0.5528321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(1.3508505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8838053) q[0];
sx q[0];
rz(-1.4865849) q[0];
sx q[0];
rz(-3.0620831) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66019085) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(-0.24812631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1119712) q[1];
sx q[1];
rz(-2.0051149) q[1];
sx q[1];
rz(-3.1236468) q[1];
rz(-pi) q[2];
rz(-1.1397347) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-0.44719493) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(0.42778095) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34925941) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(0.40339289) q[0];
x q[1];
rz(-1.595396) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(-1.3299024) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1093724) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(2.1857775) q[1];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(-2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-1.8722664) q[2];
sx q[2];
rz(-1.8282679) q[2];
sx q[2];
rz(-1.4399583) q[2];
rz(-0.11733304) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];