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
rz(-2.8647487) q[0];
sx q[0];
rz(-0.11685637) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(2.7873971) q[1];
sx q[1];
rz(8.9335557) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899333) q[0];
sx q[0];
rz(-2.058968) q[0];
sx q[0];
rz(-3.0753646) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38552706) q[2];
sx q[2];
rz(-0.51947278) q[2];
sx q[2];
rz(-2.9105718) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75233203) q[1];
sx q[1];
rz(-1.4717143) q[1];
sx q[1];
rz(-0.63431753) q[1];
x q[2];
rz(1.6412723) q[3];
sx q[3];
rz(-2.797496) q[3];
sx q[3];
rz(1.9650353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(0.22988764) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(2.9246395) q[0];
rz(3.1087061) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(-0.65509534) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0467619) q[0];
sx q[0];
rz(-2.3944986) q[0];
sx q[0];
rz(1.8674891) q[0];
rz(1.3711602) q[2];
sx q[2];
rz(-1.0773653) q[2];
sx q[2];
rz(2.9889929) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8959143) q[1];
sx q[1];
rz(-1.0686456) q[1];
sx q[1];
rz(2.4163626) q[1];
x q[2];
rz(2.3052196) q[3];
sx q[3];
rz(-1.6775838) q[3];
sx q[3];
rz(1.9143145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32404262) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(1.2980655) q[2];
rz(1.2201747) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(0.61400145) q[0];
rz(-1.7738495) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(-2.6460389) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3494095) q[0];
sx q[0];
rz(-1.111545) q[0];
sx q[0];
rz(0.71190603) q[0];
rz(0.96769963) q[2];
sx q[2];
rz(-2.357238) q[2];
sx q[2];
rz(-1.3788144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6537955) q[1];
sx q[1];
rz(-1.7572409) q[1];
sx q[1];
rz(1.814278) q[1];
rz(-2.6747781) q[3];
sx q[3];
rz(-1.2520077) q[3];
sx q[3];
rz(-2.7960645) q[3];
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
rz(2.6872915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4950824) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(0.5799154) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(-2.6204806) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189188) q[0];
sx q[0];
rz(-1.5661514) q[0];
sx q[0];
rz(1.5779431) q[0];
x q[1];
rz(-0.96984152) q[2];
sx q[2];
rz(-1.5950987) q[2];
sx q[2];
rz(1.5497335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.191491) q[1];
sx q[1];
rz(-1.9708212) q[1];
sx q[1];
rz(0.49329503) q[1];
rz(2.8981325) q[3];
sx q[3];
rz(-2.820558) q[3];
sx q[3];
rz(-0.76317549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(2.8233675) q[2];
rz(0.71277726) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(-0.036291432) q[0];
rz(2.6151784) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(-2.859419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.605114) q[0];
sx q[0];
rz(-1.8465586) q[0];
sx q[0];
rz(0.92425297) q[0];
x q[1];
rz(-2.2065978) q[2];
sx q[2];
rz(-2.3381669) q[2];
sx q[2];
rz(-2.406081) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1212738) q[1];
sx q[1];
rz(-1.500638) q[1];
sx q[1];
rz(0.70538818) q[1];
rz(-pi) q[2];
rz(-2.3822278) q[3];
sx q[3];
rz(-0.96624422) q[3];
sx q[3];
rz(-0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7669749) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-0.96595079) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58979708) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.1740603) q[0];
rz(0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(0.039965872) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55520457) q[0];
sx q[0];
rz(-0.38643906) q[0];
sx q[0];
rz(0.9903592) q[0];
x q[1];
rz(-2.0027809) q[2];
sx q[2];
rz(-2.453865) q[2];
sx q[2];
rz(-1.0291531) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0607073) q[1];
sx q[1];
rz(-2.1645912) q[1];
sx q[1];
rz(1.9869366) q[1];
rz(-2.3153855) q[3];
sx q[3];
rz(-1.7306657) q[3];
sx q[3];
rz(-1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(0.8197909) q[2];
rz(-2.36006) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(-1.8248722) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(-2.7046611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7765771) q[0];
sx q[0];
rz(-1.7611338) q[0];
sx q[0];
rz(-0.037419407) q[0];
rz(-3.0236235) q[2];
sx q[2];
rz(-1.461339) q[2];
sx q[2];
rz(0.80889672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8812411) q[1];
sx q[1];
rz(-2.2931744) q[1];
sx q[1];
rz(-2.6546225) q[1];
x q[2];
rz(2.1070621) q[3];
sx q[3];
rz(-1.9428184) q[3];
sx q[3];
rz(1.0635424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-3.0069922) q[2];
rz(-1.918119) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(-1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(-2.8002411) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-2.9827319) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(-0.44802702) q[0];
rz(0.57639052) q[2];
sx q[2];
rz(-2.6491671) q[2];
sx q[2];
rz(-2.3362291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6312069) q[1];
sx q[1];
rz(-0.34086168) q[1];
sx q[1];
rz(2.4426961) q[1];
rz(-2.587593) q[3];
sx q[3];
rz(-1.4408198) q[3];
sx q[3];
rz(-0.25773559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86026704) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(-0.9935317) q[0];
rz(1.7811071) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(-0.55307585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7599604) q[0];
sx q[0];
rz(-2.6928718) q[0];
sx q[0];
rz(2.5928709) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5037914) q[2];
sx q[2];
rz(-0.83344668) q[2];
sx q[2];
rz(0.60947567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7874541) q[1];
sx q[1];
rz(-0.95515697) q[1];
sx q[1];
rz(-1.8880185) q[1];
rz(0.93135466) q[3];
sx q[3];
rz(-1.7903312) q[3];
sx q[3];
rz(-2.1965248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-0.23703144) q[2];
rz(-1.6009181) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(1.1698394) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700664) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(3.0226829) q[0];
rz(-2.7023244) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-0.063442245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760323) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(-2.3593581) q[0];
x q[1];
rz(-1.756606) q[2];
sx q[2];
rz(-2.1868621) q[2];
sx q[2];
rz(-2.1587929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0513775) q[1];
sx q[1];
rz(-1.8974638) q[1];
sx q[1];
rz(0.0089545687) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2584553) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(-2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1305609) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(1.4981184) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(0.87286585) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(-0.68442576) q[2];
sx q[2];
rz(-1.050907) q[2];
sx q[2];
rz(2.4466865) q[2];
rz(0.11029966) q[3];
sx q[3];
rz(-0.58021373) q[3];
sx q[3];
rz(1.8662682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
