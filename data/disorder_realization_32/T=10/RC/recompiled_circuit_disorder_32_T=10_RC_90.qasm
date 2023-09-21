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
rz(2.0523235) q[0];
sx q[0];
rz(9.2637445) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007322) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(0.13829921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33978396) q[2];
sx q[2];
rz(-1.3801563) q[2];
sx q[2];
rz(1.9356188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6160994) q[1];
sx q[1];
rz(-1.6615189) q[1];
sx q[1];
rz(0.70848042) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22827893) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(-1.3397863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.455362) q[0];
rz(-2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(2.9532202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8375081) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(1.1156032) q[1];
rz(-pi) q[2];
rz(1.7277667) q[3];
sx q[3];
rz(-1.7236712) q[3];
sx q[3];
rz(-2.4863941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-0.64003402) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.144369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17830081) q[0];
sx q[0];
rz(-0.68094567) q[0];
sx q[0];
rz(0.91239022) q[0];
x q[1];
rz(-0.47282131) q[2];
sx q[2];
rz(-0.46337767) q[2];
sx q[2];
rz(-2.9544427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5971165) q[1];
sx q[1];
rz(-0.58864486) q[1];
sx q[1];
rz(2.4878923) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5738437) q[3];
sx q[3];
rz(-1.8584195) q[3];
sx q[3];
rz(-2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-0.90467492) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7870165) q[0];
sx q[0];
rz(-1.0542608) q[0];
sx q[0];
rz(-1.2660962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86187141) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0992972) q[1];
sx q[1];
rz(-0.18810939) q[1];
sx q[1];
rz(-0.52109615) q[1];
x q[2];
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
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(2.7899182) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2131166) q[0];
sx q[0];
rz(-0.94841829) q[0];
sx q[0];
rz(-1.2147796) q[0];
x q[1];
rz(0.93514498) q[2];
sx q[2];
rz(-1.9872553) q[2];
sx q[2];
rz(-2.2098429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56375757) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(-1.5549591) q[1];
x q[2];
rz(-2.8204927) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(-1.0995964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8655411) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(1.0827351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8855368) q[2];
sx q[2];
rz(-2.5883) q[2];
sx q[2];
rz(-1.7994583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53193608) q[1];
sx q[1];
rz(-1.7488639) q[1];
sx q[1];
rz(-2.460536) q[1];
rz(1.2332637) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(-0.28901643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(-0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-2.9842916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10282117) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(-3.0022013) q[0];
x q[1];
rz(-0.21643164) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(-0.12491465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8017756) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(-2.9031309) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5141684) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5921708) q[0];
sx q[0];
rz(-0.54245078) q[0];
sx q[0];
rz(-0.011499238) q[0];
rz(-2.9759334) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(0.13656244) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3775023) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(-0.80470316) q[1];
rz(-pi) q[2];
rz(0.71473748) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(-1.8470256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(0.088767178) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.3508505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8838053) q[0];
sx q[0];
rz(-1.6550078) q[0];
sx q[0];
rz(-0.079509602) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0355989) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(1.7352599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6752154) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(2.0051763) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3488995) q[3];
sx q[3];
rz(-1.6716692) q[3];
sx q[3];
rz(-1.6290806) q[3];
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
rz(1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(2.9719877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113838) q[0];
sx q[0];
rz(-1.9281403) q[0];
sx q[0];
rz(-2.0765199) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42189235) q[2];
sx q[2];
rz(-1.5483529) q[2];
sx q[2];
rz(-0.23082146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5369397) q[1];
sx q[1];
rz(-0.61652684) q[1];
sx q[1];
rz(1.489868) q[1];
rz(-0.1428991) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(2.4492241) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(2.0326116) q[3];
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
rz(-pi/2) q[0];
sx q[2];
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