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
rz(1.5155091) q[0];
sx q[0];
rz(3.4184366) q[0];
sx q[0];
rz(9.3079216) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(-2.6503704) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913585) q[0];
sx q[0];
rz(-1.629279) q[0];
sx q[0];
rz(2.0598777) q[0];
x q[1];
rz(-0.38552706) q[2];
sx q[2];
rz(-0.51947278) q[2];
sx q[2];
rz(0.2310209) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1894816) q[1];
sx q[1];
rz(-2.5006388) q[1];
sx q[1];
rz(2.9753995) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1163636) q[3];
sx q[3];
rz(-1.2275891) q[3];
sx q[3];
rz(2.0398839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-0.22988764) q[2];
rz(-0.057279438) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(-2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(0.21695319) q[0];
rz(-0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(2.4864973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74467662) q[0];
sx q[0];
rz(-1.7707877) q[0];
sx q[0];
rz(-2.295664) q[0];
rz(-pi) q[1];
rz(-2.7883165) q[2];
sx q[2];
rz(-2.6124119) q[2];
sx q[2];
rz(-0.25111094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0634708) q[1];
sx q[1];
rz(-2.1912716) q[1];
sx q[1];
rz(2.2038151) q[1];
rz(-0.83637303) q[3];
sx q[3];
rz(-1.4640088) q[3];
sx q[3];
rz(1.2272782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(-1.2201747) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.87797457) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(-0.61400145) q[0];
rz(1.7738495) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(2.6460389) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3494095) q[0];
sx q[0];
rz(-1.111545) q[0];
sx q[0];
rz(0.71190603) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96769963) q[2];
sx q[2];
rz(-2.357238) q[2];
sx q[2];
rz(-1.3788144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12901356) q[1];
sx q[1];
rz(-1.3316175) q[1];
sx q[1];
rz(-2.9496179) q[1];
rz(-0.63276799) q[3];
sx q[3];
rz(-2.583021) q[3];
sx q[3];
rz(2.4726923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-2.5490254) q[2];
rz(-2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-0.0058280514) q[0];
rz(-2.5616772) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(-2.6204806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0590044) q[0];
sx q[0];
rz(-0.0085235611) q[0];
sx q[0];
rz(-2.1471291) q[0];
rz(2.1717511) q[2];
sx q[2];
rz(-1.5464939) q[2];
sx q[2];
rz(-1.5497335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5856984) q[1];
sx q[1];
rz(-1.1194849) q[1];
sx q[1];
rz(2.0183619) q[1];
rz(0.24346015) q[3];
sx q[3];
rz(-0.32103466) q[3];
sx q[3];
rz(2.3784172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(-2.4288154) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37079674) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(0.2821736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9045893) q[0];
sx q[0];
rz(-2.1891199) q[0];
sx q[0];
rz(-0.34070054) q[0];
rz(-pi) q[1];
rz(-2.2065978) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(-0.73551169) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4683405) q[1];
sx q[1];
rz(-0.70827228) q[1];
sx q[1];
rz(0.10797031) q[1];
rz(-2.3318841) q[3];
sx q[3];
rz(-2.1729709) q[3];
sx q[3];
rz(2.6701327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(3.0401163) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(3.0275893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(-1.9675323) q[0];
rz(-0.38408285) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-3.1016268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5863881) q[0];
sx q[0];
rz(-0.38643906) q[0];
sx q[0];
rz(-2.1512335) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92982473) q[2];
sx q[2];
rz(-1.8397959) q[2];
sx q[2];
rz(0.19942927) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0607073) q[1];
sx q[1];
rz(-2.1645912) q[1];
sx q[1];
rz(-1.1546561) q[1];
x q[2];
rz(0.82620718) q[3];
sx q[3];
rz(-1.4109269) q[3];
sx q[3];
rz(-1.2197956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(-0.8197909) q[2];
rz(2.36006) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6509197) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.3167205) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(-2.7046611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56037134) q[0];
sx q[0];
rz(-2.9476555) q[0];
sx q[0];
rz(-1.3790129) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11796911) q[2];
sx q[2];
rz(-1.6802537) q[2];
sx q[2];
rz(-0.80889672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26035151) q[1];
sx q[1];
rz(-2.2931744) q[1];
sx q[1];
rz(0.48697014) q[1];
rz(2.1070621) q[3];
sx q[3];
rz(-1.1987743) q[3];
sx q[3];
rz(2.0780502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(3.0069922) q[2];
rz(-1.2234737) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80158919) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(0.18184161) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(0.15886074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5349992) q[0];
sx q[0];
rz(-2.6887021) q[0];
sx q[0];
rz(2.9843828) q[0];
rz(-pi) q[1];
rz(1.2863289) q[2];
sx q[2];
rz(-1.9783696) q[2];
sx q[2];
rz(1.4409232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4113058) q[1];
sx q[1];
rz(-1.3540233) q[1];
sx q[1];
rz(2.8764399) q[1];
rz(-pi) q[2];
rz(-0.55399968) q[3];
sx q[3];
rz(-1.7007728) q[3];
sx q[3];
rz(-0.25773559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2813256) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-0.88480985) q[2];
rz(-1.4191215) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.849843) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(-0.9935317) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(0.55307585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38163227) q[0];
sx q[0];
rz(-2.6928718) q[0];
sx q[0];
rz(-2.5928709) q[0];
x q[1];
rz(-1.5037914) q[2];
sx q[2];
rz(-0.83344668) q[2];
sx q[2];
rz(-0.60947567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40402624) q[1];
sx q[1];
rz(-1.8282923) q[1];
sx q[1];
rz(2.5016258) q[1];
rz(1.2129799) q[3];
sx q[3];
rz(-0.67104895) q[3];
sx q[3];
rz(0.34085654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(1.1698394) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700664) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(2.7023244) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(3.0781504) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3812696) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(2.3593581) q[0];
rz(0.25524885) q[2];
sx q[2];
rz(-2.5016157) q[2];
sx q[2];
rz(0.6682804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0902151) q[1];
sx q[1];
rz(-1.8974638) q[1];
sx q[1];
rz(3.1326381) q[1];
x q[2];
rz(0.48096913) q[3];
sx q[3];
rz(-2.0837373) q[3];
sx q[3];
rz(-0.45403593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.6434742) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3907923) q[0];
sx q[0];
rz(-1.1815429) q[0];
sx q[0];
rz(-1.0829157) q[0];
rz(0.87286585) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(0.93449705) q[2];
sx q[2];
rz(-2.1515982) q[2];
sx q[2];
rz(1.2610255) q[2];
rz(-0.57742248) q[3];
sx q[3];
rz(-1.631177) q[3];
sx q[3];
rz(-2.9384818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
