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
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0954477) q[0];
sx q[0];
rz(-0.60337043) q[0];
sx q[0];
rz(-1.3674111) q[0];
x q[1];
rz(0.52486323) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(-2.2847069) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-2.2757581) q[1];
sx q[1];
rz(1.4515619) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(-1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(0.57759181) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47050414) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(1.455362) q[0];
rz(-pi) q[1];
rz(-1.8446484) q[2];
sx q[2];
rz(-1.8631862) q[2];
sx q[2];
rz(-1.3016303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6985804) q[1];
sx q[1];
rz(-0.46291446) q[1];
sx q[1];
rz(1.7673311) q[1];
rz(-0.15474774) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(-0.043116365) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.144369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(2.6813566) q[0];
x q[1];
rz(1.3470596) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(-2.4350016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5971165) q[1];
sx q[1];
rz(-0.58864486) q[1];
sx q[1];
rz(2.4878923) q[1];
rz(-pi) q[2];
rz(0.010300962) q[3];
sx q[3];
rz(-2.8539538) q[3];
sx q[3];
rz(2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(0.71980643) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(0.71587193) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199812) q[0];
sx q[0];
rz(-0.5926168) q[0];
sx q[0];
rz(0.48595925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9348346) q[2];
sx q[2];
rz(-1.0996498) q[2];
sx q[2];
rz(-0.83540321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.042295481) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(-0.52109615) q[1];
rz(2.2654815) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(-2.0640404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.7129718) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-2.7899182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-2.48824) q[0];
rz(0.93514498) q[2];
sx q[2];
rz(-1.9872553) q[2];
sx q[2];
rz(-2.2098429) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5366718) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(3.1035963) q[1];
rz(-pi) q[2];
rz(2.0760173) q[3];
sx q[3];
rz(-1.2875644) q[3];
sx q[3];
rz(0.035988228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(1.0995964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27605155) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(-1.0827351) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(-0.97666937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1814327) q[1];
sx q[1];
rz(-0.90248855) q[1];
sx q[1];
rz(-1.3431576) q[1];
rz(-0.78824708) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.9842916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0387715) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(-3.0022013) q[0];
x q[1];
rz(-1.8087177) q[2];
sx q[2];
rz(-0.750713) q[2];
sx q[2];
rz(-0.44516341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3992608) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(1.8982235) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62742426) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(-1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(2.9390745) q[3];
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
rz(-1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5787443) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(-1.5638652) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9759334) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(0.13656244) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(0.80470316) q[1];
x q[2];
rz(-1.9803489) q[3];
sx q[3];
rz(-2.000994) q[3];
sx q[3];
rz(-0.48914117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.3508505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8352855) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.6552734) q[0];
x q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-2.2293408) q[2];
sx q[2];
rz(-1.8694307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1119712) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(-3.1236468) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0302089) q[0];
sx q[0];
rz(-1.9281403) q[0];
sx q[0];
rz(-1.0650728) q[0];
x q[1];
rz(3.0868297) q[2];
sx q[2];
rz(-2.7191396) q[2];
sx q[2];
rz(1.3899318) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.032220275) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(0.95581518) q[1];
x q[2];
rz(1.8923558) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
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
rz(-1.7789755) q[3];
sx q[3];
rz(-1.4559742) q[3];
sx q[3];
rz(-2.5696587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];