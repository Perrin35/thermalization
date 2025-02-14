OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2990155) q[0];
sx q[0];
rz(3.314078) q[0];
sx q[0];
rz(9.4077851) q[0];
rz(-2.6301771) q[1];
sx q[1];
rz(-2.6452682) q[1];
sx q[1];
rz(2.6561148) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66390178) q[0];
sx q[0];
rz(-1.6833651) q[0];
sx q[0];
rz(-1.4590053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67022143) q[2];
sx q[2];
rz(-2.1677758) q[2];
sx q[2];
rz(-0.94979294) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37066165) q[1];
sx q[1];
rz(-1.071613) q[1];
sx q[1];
rz(0.73727495) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0270013) q[3];
sx q[3];
rz(-1.262731) q[3];
sx q[3];
rz(1.9730572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7444676) q[2];
sx q[2];
rz(-0.10245377) q[2];
sx q[2];
rz(2.3490119) q[2];
rz(-2.1553195) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(3.0084394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.88696402) q[0];
sx q[0];
rz(-1.6034842) q[0];
sx q[0];
rz(3.0358553) q[0];
rz(2.3836783) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(2.6656718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0237813) q[0];
sx q[0];
rz(-2.887368) q[0];
sx q[0];
rz(2.7455968) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7143102) q[2];
sx q[2];
rz(-2.0594308) q[2];
sx q[2];
rz(1.2554864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5761199) q[1];
sx q[1];
rz(-0.58153897) q[1];
sx q[1];
rz(-1.3493058) q[1];
x q[2];
rz(0.29637166) q[3];
sx q[3];
rz(-2.1934273) q[3];
sx q[3];
rz(-2.6105256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14625749) q[2];
sx q[2];
rz(-0.47440752) q[2];
sx q[2];
rz(-1.8246626) q[2];
rz(-2.3138192) q[3];
sx q[3];
rz(-1.0931949) q[3];
sx q[3];
rz(-0.83724666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4227609) q[0];
sx q[0];
rz(-1.3439002) q[0];
sx q[0];
rz(-1.2368917) q[0];
rz(-1.749136) q[1];
sx q[1];
rz(-1.4207276) q[1];
sx q[1];
rz(-0.69033355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.515432) q[0];
sx q[0];
rz(-1.0179839) q[0];
sx q[0];
rz(1.0286858) q[0];
rz(1.5788011) q[2];
sx q[2];
rz(-2.6782616) q[2];
sx q[2];
rz(-0.066502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18728367) q[1];
sx q[1];
rz(-0.25588965) q[1];
sx q[1];
rz(-2.4268716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.407904) q[3];
sx q[3];
rz(-2.0533959) q[3];
sx q[3];
rz(-0.35728282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6231125) q[2];
sx q[2];
rz(-1.6192351) q[2];
sx q[2];
rz(-0.073866455) q[2];
rz(1.9833924) q[3];
sx q[3];
rz(-1.9246212) q[3];
sx q[3];
rz(-1.4069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.912792) q[0];
sx q[0];
rz(-0.055963628) q[0];
sx q[0];
rz(-2.9486616) q[0];
rz(-1.2436766) q[1];
sx q[1];
rz(-1.396215) q[1];
sx q[1];
rz(2.9085433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4547494) q[0];
sx q[0];
rz(-0.17486869) q[0];
sx q[0];
rz(-0.57342477) q[0];
x q[1];
rz(0.86200447) q[2];
sx q[2];
rz(-1.5420015) q[2];
sx q[2];
rz(-2.6301602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.081542) q[1];
sx q[1];
rz(-2.8257408) q[1];
sx q[1];
rz(2.9595023) q[1];
rz(-pi) q[2];
rz(-0.47073029) q[3];
sx q[3];
rz(-2.1370557) q[3];
sx q[3];
rz(1.2486718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9243246) q[2];
sx q[2];
rz(-1.7622207) q[2];
sx q[2];
rz(3.0775089) q[2];
rz(-0.25121769) q[3];
sx q[3];
rz(-0.46291864) q[3];
sx q[3];
rz(0.29138756) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578882) q[0];
sx q[0];
rz(-1.2936445) q[0];
sx q[0];
rz(1.8573014) q[0];
rz(-3.1341556) q[1];
sx q[1];
rz(-1.1859272) q[1];
sx q[1];
rz(0.95710212) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30516) q[0];
sx q[0];
rz(-1.1344988) q[0];
sx q[0];
rz(2.7899105) q[0];
x q[1];
rz(-0.77218036) q[2];
sx q[2];
rz(-0.77754687) q[2];
sx q[2];
rz(2.5302056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51077183) q[1];
sx q[1];
rz(-1.242852) q[1];
sx q[1];
rz(-1.6751218) q[1];
x q[2];
rz(-0.65517218) q[3];
sx q[3];
rz(-1.4219033) q[3];
sx q[3];
rz(-0.53194351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0127504) q[2];
sx q[2];
rz(-0.70657554) q[2];
sx q[2];
rz(-1.0303222) q[2];
rz(0.71850592) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(2.4350186) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71501032) q[0];
sx q[0];
rz(-0.74786818) q[0];
sx q[0];
rz(-0.36710516) q[0];
rz(2.7857065) q[1];
sx q[1];
rz(-1.117319) q[1];
sx q[1];
rz(-1.5497367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0192896) q[0];
sx q[0];
rz(-0.82989026) q[0];
sx q[0];
rz(1.6229963) q[0];
rz(-pi) q[1];
rz(2.0075304) q[2];
sx q[2];
rz(-1.6730089) q[2];
sx q[2];
rz(-1.5222766) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15799668) q[1];
sx q[1];
rz(-1.776773) q[1];
sx q[1];
rz(-0.41974824) q[1];
x q[2];
rz(-2.995302) q[3];
sx q[3];
rz(-2.6669569) q[3];
sx q[3];
rz(2.4865347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97673577) q[2];
sx q[2];
rz(-1.2931436) q[2];
sx q[2];
rz(3.1381651) q[2];
rz(0.88611832) q[3];
sx q[3];
rz(-2.0485853) q[3];
sx q[3];
rz(-1.4440943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057328) q[0];
sx q[0];
rz(-1.6678565) q[0];
sx q[0];
rz(0.71520299) q[0];
rz(-0.12229478) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(0.14762793) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58260949) q[0];
sx q[0];
rz(-1.9096073) q[0];
sx q[0];
rz(-2.7937583) q[0];
rz(-pi) q[1];
rz(0.4345456) q[2];
sx q[2];
rz(-1.927863) q[2];
sx q[2];
rz(0.94851102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9761889) q[1];
sx q[1];
rz(-0.87748412) q[1];
sx q[1];
rz(-2.90392) q[1];
rz(3.0711898) q[3];
sx q[3];
rz(-1.284666) q[3];
sx q[3];
rz(-2.4918258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82835308) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(3.0301136) q[2];
rz(-0.76505032) q[3];
sx q[3];
rz(-1.7181516) q[3];
sx q[3];
rz(-3.1077207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400804) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(1.0028268) q[0];
rz(2.3560933) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(-1.921382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1357546) q[0];
sx q[0];
rz(-1.5501115) q[0];
sx q[0];
rz(-0.51405859) q[0];
rz(-1.6613879) q[2];
sx q[2];
rz(-1.5085417) q[2];
sx q[2];
rz(-1.412815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0682098) q[1];
sx q[1];
rz(-2.073895) q[1];
sx q[1];
rz(0.36243172) q[1];
x q[2];
rz(-1.7986913) q[3];
sx q[3];
rz(-2.493495) q[3];
sx q[3];
rz(-0.15290393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9057374) q[2];
sx q[2];
rz(-1.7451655) q[2];
sx q[2];
rz(-0.033795707) q[2];
rz(-2.3679521) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(2.0788367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3646024) q[0];
sx q[0];
rz(-0.62043014) q[0];
sx q[0];
rz(-0.34238368) q[0];
rz(1.4211593) q[1];
sx q[1];
rz(-0.8232638) q[1];
sx q[1];
rz(-1.8470496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7086805) q[0];
sx q[0];
rz(-1.4231235) q[0];
sx q[0];
rz(-2.7765034) q[0];
rz(-pi) q[1];
rz(2.2850288) q[2];
sx q[2];
rz(-0.98443778) q[2];
sx q[2];
rz(2.7991653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0659938) q[1];
sx q[1];
rz(-0.79388777) q[1];
sx q[1];
rz(0.17236472) q[1];
x q[2];
rz(-2.8043069) q[3];
sx q[3];
rz(-1.4082068) q[3];
sx q[3];
rz(0.43460571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2932059) q[2];
sx q[2];
rz(-1.3853955) q[2];
sx q[2];
rz(2.6050341) q[2];
rz(-0.41283354) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(-1.1134953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31529108) q[0];
sx q[0];
rz(-1.8719801) q[0];
sx q[0];
rz(-1.4971365) q[0];
rz(2.3313088) q[1];
sx q[1];
rz(-1.5565926) q[1];
sx q[1];
rz(2.2946045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627914) q[0];
sx q[0];
rz(-1.7399995) q[0];
sx q[0];
rz(-0.16296084) q[0];
rz(-pi) q[1];
rz(-0.92310793) q[2];
sx q[2];
rz(-2.0294242) q[2];
sx q[2];
rz(-0.066740008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0125186) q[1];
sx q[1];
rz(-1.8400987) q[1];
sx q[1];
rz(2.49519) q[1];
x q[2];
rz(-0.096135898) q[3];
sx q[3];
rz(-1.8435974) q[3];
sx q[3];
rz(0.48855272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.102313) q[2];
sx q[2];
rz(-1.8259093) q[2];
sx q[2];
rz(-0.63423356) q[2];
rz(-2.5041194) q[3];
sx q[3];
rz(-2.0135148) q[3];
sx q[3];
rz(2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98916003) q[0];
sx q[0];
rz(-1.1897054) q[0];
sx q[0];
rz(-1.1783896) q[0];
rz(2.4329026) q[1];
sx q[1];
rz(-1.7955753) q[1];
sx q[1];
rz(1.2830455) q[1];
rz(-0.75015776) q[2];
sx q[2];
rz(-2.6476319) q[2];
sx q[2];
rz(-2.8110197) q[2];
rz(-1.9509964) q[3];
sx q[3];
rz(-0.96933848) q[3];
sx q[3];
rz(-2.5179767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
