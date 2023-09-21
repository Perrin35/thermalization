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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8344603) q[0];
sx q[0];
rz(-1.6856598) q[0];
sx q[0];
rz(-2.1644724) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6167294) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(0.85688574) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.4515619) q[1];
rz(1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-0.28764763) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.6862306) q[0];
rz(0.73194506) q[2];
sx q[2];
rz(-2.7436896) q[2];
sx q[2];
rz(-2.0741472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6985804) q[1];
sx q[1];
rz(-0.46291446) q[1];
sx q[1];
rz(1.7673311) q[1];
rz(-pi) q[2];
rz(-2.9868449) q[3];
sx q[3];
rz(-1.4156716) q[3];
sx q[3];
rz(2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(0.46023603) q[0];
rz(-2.6687713) q[2];
sx q[2];
rz(-2.678215) q[2];
sx q[2];
rz(0.18714999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5444762) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(-2.4878923) q[1];
x q[2];
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
rz(-pi/2) q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(2.604565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86187141) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(-0.18415235) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6550733) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(1.4763114) q[1];
rz(-pi) q[2];
x q[2];
rz(2.643814) q[3];
sx q[3];
rz(-1.0504424) q[3];
sx q[3];
rz(-1.9073515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(2.3748421) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2131166) q[0];
sx q[0];
rz(-0.94841829) q[0];
sx q[0];
rz(1.2147796) q[0];
rz(-2.2064477) q[2];
sx q[2];
rz(-1.1543373) q[2];
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
rz(1.000948) q[1];
sx q[1];
rz(-1.556177) q[1];
sx q[1];
rz(2.7468202) q[1];
rz(-pi) q[2];
rz(1.0655754) q[3];
sx q[3];
rz(-1.8540283) q[3];
sx q[3];
rz(0.035988228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(0.41082541) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(0.97666937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1814327) q[1];
sx q[1];
rz(-0.90248855) q[1];
sx q[1];
rz(-1.798435) q[1];
rz(1.908329) q[3];
sx q[3];
rz(-2.3300155) q[3];
sx q[3];
rz(-0.28901643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(0.15730102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10282117) q[0];
sx q[0];
rz(-2.4024995) q[0];
sx q[0];
rz(-3.0022013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21643164) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(-3.016678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.339817) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(-0.2384618) q[1];
rz(-1.2792475) q[3];
sx q[3];
rz(-0.96372094) q[3];
sx q[3];
rz(2.8214422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(2.8483134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1300668) q[0];
sx q[0];
rz(-1.5767326) q[0];
sx q[0];
rz(-2.5991711) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9759334) q[2];
sx q[2];
rz(-2.207901) q[2];
sx q[2];
rz(-3.0050302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5756702) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-2.5887606) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1612438) q[3];
sx q[3];
rz(-2.000994) q[3];
sx q[3];
rz(0.48914117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(0.5789825) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042689) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.7907422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8838053) q[0];
sx q[0];
rz(-1.6550078) q[0];
sx q[0];
rz(-0.079509602) q[0];
rz(0.10599372) q[2];
sx q[2];
rz(-2.4786737) q[2];
sx q[2];
rz(1.7352599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6752154) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(1.1364163) q[1];
rz(1.3488995) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-0.16960493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34925941) q[0];
sx q[0];
rz(-1.0997286) q[0];
sx q[0];
rz(-2.7381998) q[0];
x q[1];
rz(-3.0868297) q[2];
sx q[2];
rz(-2.7191396) q[2];
sx q[2];
rz(1.7516608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5369397) q[1];
sx q[1];
rz(-0.61652684) q[1];
sx q[1];
rz(-1.489868) q[1];
x q[2];
rz(1.8923558) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.8511843) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-0.26906536) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(3.0242596) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
