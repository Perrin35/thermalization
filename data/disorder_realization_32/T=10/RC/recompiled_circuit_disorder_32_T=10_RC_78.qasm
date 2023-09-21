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
rz(-2.6452126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8007322) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(3.0032934) q[0];
rz(-pi) q[1];
rz(-0.33978396) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(-1.2059739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9909518) q[1];
sx q[1];
rz(-2.4283263) q[1];
sx q[1];
rz(3.002682) q[1];
x q[2];
rz(-2.9133137) q[3];
sx q[3];
rz(-2.7242959) q[3];
sx q[3];
rz(-1.3397863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.577852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801115) q[0];
sx q[0];
rz(-2.8948088) q[0];
sx q[0];
rz(2.6632975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(0.18837243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9175375) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(0.097150306) q[1];
rz(1.413826) q[3];
sx q[3];
rz(-1.4179215) q[3];
sx q[3];
rz(-2.4863941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-0.18297718) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5367033) q[0];
sx q[0];
rz(-2.0920144) q[0];
sx q[0];
rz(0.46023603) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6687713) q[2];
sx q[2];
rz(-0.46337767) q[2];
sx q[2];
rz(-2.9544427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5444762) q[1];
sx q[1];
rz(-0.58864486) q[1];
sx q[1];
rz(0.65370037) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5677489) q[3];
sx q[3];
rz(-1.2831732) q[3];
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
rz(-2.7147839) q[3];
sx q[3];
rz(-0.75511801) q[3];
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
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(2.604565) q[0];
x q[1];
rz(-2.2067581) q[2];
sx q[2];
rz(-1.0996498) q[2];
sx q[2];
rz(-0.83540321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0996453) q[1];
sx q[1];
rz(-1.6640267) q[1];
sx q[1];
rz(-0.16361841) q[1];
rz(2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(1.9073515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.283678) q[0];
sx q[0];
rz(0.65335269) q[0];
rz(-0.93047662) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(-3.004068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5366718) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(0.037996304) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0293803) q[3];
sx q[3];
rz(-0.57313985) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-2.0419962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106364) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(-2.7307672) q[0];
rz(-pi) q[1];
rz(-2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(2.1649233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53193608) q[1];
sx q[1];
rz(-1.7488639) q[1];
sx q[1];
rz(-0.68105662) q[1];
x q[2];
rz(0.78824708) q[3];
sx q[3];
rz(-1.8133834) q[3];
sx q[3];
rz(-1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.5299861) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(0.15730102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10282117) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(3.0022013) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8087177) q[2];
sx q[2];
rz(-2.3908797) q[2];
sx q[2];
rz(2.6964292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7423319) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(1.2433692) q[1];
rz(-pi) q[2];
rz(-1.2792475) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5787443) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(1.5777274) q[0];
rz(-pi) q[1];
rz(-2.9759334) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(-3.0050302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20626942) q[1];
sx q[1];
rz(-0.85201293) q[1];
sx q[1];
rz(2.1410336) q[1];
rz(-pi) q[2];
rz(-0.71473748) q[3];
sx q[3];
rz(-2.5566031) q[3];
sx q[3];
rz(1.294567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-3.0528255) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(-0.44257277) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.7907422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8352855) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.4863192) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10599372) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(-1.7352599) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0296214) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(0.017945826) q[1];
rz(-1.1397347) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.3859008) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302089) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(-1.0650728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.054762997) q[2];
sx q[2];
rz(-0.4224531) q[2];
sx q[2];
rz(1.3899318) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1093724) q[1];
sx q[1];
rz(-1.5240372) q[1];
sx q[1];
rz(-0.95581518) q[1];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(2.3615169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(1.2693263) q[2];
sx q[2];
rz(-1.8282679) q[2];
sx q[2];
rz(-1.4399583) q[2];
rz(1.3626171) q[3];
sx q[3];
rz(-1.4559742) q[3];
sx q[3];
rz(-2.5696587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
