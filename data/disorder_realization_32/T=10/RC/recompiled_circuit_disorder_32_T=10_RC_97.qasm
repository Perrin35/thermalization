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
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3408605) q[0];
sx q[0];
rz(-2.1600318) q[0];
sx q[0];
rz(3.0032934) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8018087) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(1.9356188) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-pi) q[2];
rz(-1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(-2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(2.853945) q[2];
rz(-1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(-2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(-1.5637406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(-1.455362) q[0];
x q[1];
rz(1.8446484) q[2];
sx q[2];
rz(-1.2784064) q[2];
sx q[2];
rz(-1.3016303) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4430122) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(-1.3742616) q[1];
rz(-pi) q[2];
rz(-0.79264499) q[3];
sx q[3];
rz(-0.21867293) q[3];
sx q[3];
rz(-1.45989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.144369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2079175) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(1.0008706) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3470596) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(-0.70659107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.548) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(-2.6542632) q[1];
rz(-pi) q[2];
rz(-0.28762443) q[3];
sx q[3];
rz(-1.5737185) q[3];
sx q[3];
rz(-0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(0.71587193) q[0];
rz(-2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0621322) q[0];
sx q[0];
rz(-1.8347164) q[0];
sx q[0];
rz(0.53702766) q[0];
x q[1];
rz(2.5771192) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(-2.0828473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0419473) q[1];
sx q[1];
rz(-1.6640267) q[1];
sx q[1];
rz(0.16361841) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99289258) q[3];
sx q[3];
rz(-1.143647) q[3];
sx q[3];
rz(3.0689193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(1.4286208) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85580101) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(2.48824) q[0];
x q[1];
rz(0.50260966) q[2];
sx q[2];
rz(-2.1447499) q[2];
sx q[2];
rz(-0.9290907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5366718) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(0.037996304) q[1];
x q[2];
rz(2.1122123) q[3];
sx q[3];
rz(-2.5684528) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9436283) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11739843) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(1.0995964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106364) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(0.41082541) q[0];
rz(-2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(2.1649233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53193608) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(0.68105662) q[1];
rz(2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-2.9842916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10282117) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(-3.0022013) q[0];
rz(-pi) q[1];
rz(-2.925161) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(0.12491465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8017756) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(0.2384618) q[1];
rz(-pi) q[2];
rz(1.2792475) q[3];
sx q[3];
rz(-0.96372094) q[3];
sx q[3];
rz(0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56284833) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(1.5777274) q[0];
x q[1];
rz(-1.3515527) q[2];
sx q[2];
rz(-0.65537894) q[2];
sx q[2];
rz(0.13742451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5659225) q[1];
sx q[1];
rz(-2.2568963) q[1];
sx q[1];
rz(-0.5528321) q[1];
x q[2];
rz(1.9803489) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(2.6524515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.7907422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25778739) q[0];
sx q[0];
rz(-1.6550078) q[0];
sx q[0];
rz(3.0620831) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66019085) q[2];
sx q[2];
rz(-1.635951) q[2];
sx q[2];
rz(-0.24812631) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0296214) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(3.1236468) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10339046) q[3];
sx q[3];
rz(-1.7915465) q[3];
sx q[3];
rz(-3.0605928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0302089) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(1.0650728) q[0];
rz(-0.42189235) q[2];
sx q[2];
rz(-1.5483529) q[2];
sx q[2];
rz(-2.9107712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-0.95658703) q[1];
sx q[1];
rz(3.0843656) q[1];
rz(0.1428991) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(0.83574502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.8511843) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(0.84531534) q[2];
sx q[2];
rz(-2.7477063) q[2];
sx q[2];
rz(0.817) q[2];
rz(-1.061822) q[3];
sx q[3];
rz(-2.90425) q[3];
sx q[3];
rz(-1.4958285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
