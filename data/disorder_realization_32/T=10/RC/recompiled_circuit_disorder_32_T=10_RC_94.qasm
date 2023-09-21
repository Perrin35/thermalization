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
rz(-1.5313907) q[1];
sx q[1];
rz(-2.664497) q[1];
sx q[1];
rz(2.6452126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0954477) q[0];
sx q[0];
rz(-0.60337043) q[0];
sx q[0];
rz(1.3674111) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33978396) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(-1.9356188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9909518) q[1];
sx q[1];
rz(-2.4283263) q[1];
sx q[1];
rz(0.13891061) q[1];
x q[2];
rz(2.9133137) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(-1.3397863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9614812) q[0];
sx q[0];
rz(-2.8948088) q[0];
sx q[0];
rz(-0.47829511) q[0];
rz(-1.8446484) q[2];
sx q[2];
rz(-1.8631862) q[2];
sx q[2];
rz(-1.3016303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8375081) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(2.0259894) q[1];
rz(-pi) q[2];
rz(-2.3489477) q[3];
sx q[3];
rz(-0.21867293) q[3];
sx q[3];
rz(-1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.95298195) q[1];
sx q[1];
rz(-1.9972237) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9336752) q[0];
sx q[0];
rz(-1.9662004) q[0];
sx q[0];
rz(-1.0008706) q[0];
rz(-2.6687713) q[2];
sx q[2];
rz(-2.678215) q[2];
sx q[2];
rz(0.18714999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3414587) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(1.9564499) q[1];
rz(-pi) q[2];
rz(-0.28762443) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1220876) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(2.148927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794605) q[0];
sx q[0];
rz(-1.8347164) q[0];
sx q[0];
rz(-0.53702766) q[0];
x q[1];
rz(-0.9348346) q[2];
sx q[2];
rz(-1.0996498) q[2];
sx q[2];
rz(0.83540321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6550733) q[1];
sx q[1];
rz(-1.4078948) q[1];
sx q[1];
rz(-1.6652813) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99289258) q[3];
sx q[3];
rz(-1.143647) q[3];
sx q[3];
rz(-0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(2.7899182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-0.65335269) q[0];
rz(0.93514498) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(2.2098429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56375757) q[1];
sx q[1];
rz(-1.1760684) q[1];
sx q[1];
rz(1.5549591) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0655754) q[3];
sx q[3];
rz(-1.2875644) q[3];
sx q[3];
rz(-3.1056044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-1.126948) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11739843) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(2.0419962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.106364) q[0];
sx q[0];
rz(-2.0237676) q[0];
sx q[0];
rz(-2.7307672) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0397644) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(-0.49887564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1814327) q[1];
sx q[1];
rz(-2.2391041) q[1];
sx q[1];
rz(1.3431576) q[1];
x q[2];
rz(-2.3533456) q[3];
sx q[3];
rz(-1.8133834) q[3];
sx q[3];
rz(1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-2.9842916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-1.6645263) q[0];
sx q[0];
rz(-2.4073497) q[0];
rz(-pi) q[1];
rz(-1.332875) q[2];
sx q[2];
rz(-0.750713) q[2];
sx q[2];
rz(0.44516341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3992608) q[1];
sx q[1];
rz(-0.6472339) q[1];
sx q[1];
rz(1.8982235) q[1];
x q[2];
rz(-1.8623452) q[3];
sx q[3];
rz(-0.96372094) q[3];
sx q[3];
rz(-2.8214422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(2.9390745) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.3515527) q[2];
sx q[2];
rz(-0.65537894) q[2];
sx q[2];
rz(3.0041681) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5756702) q[1];
sx q[1];
rz(-2.2568963) q[1];
sx q[1];
rz(-2.5887606) q[1];
x q[2];
rz(-0.46383143) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(-2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(2.8914715) q[3];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.3508505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8838053) q[0];
sx q[0];
rz(-1.6550078) q[0];
sx q[0];
rz(-3.0620831) q[0];
rz(-pi) q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(1.8694307) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46637725) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(-1.1364163) q[1];
x q[2];
rz(-1.7926932) q[3];
sx q[3];
rz(-1.6716692) q[3];
sx q[3];
rz(-1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-0.16960493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7923332) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(-0.40339289) q[0];
rz(-pi) q[1];
rz(-0.42189235) q[2];
sx q[2];
rz(-1.5932398) q[2];
sx q[2];
rz(2.9107712) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.032220275) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(-2.1857775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(-2.3615169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(-2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8725273) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(-1.3626171) q[3];
sx q[3];
rz(-1.6856185) q[3];
sx q[3];
rz(0.57193397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];