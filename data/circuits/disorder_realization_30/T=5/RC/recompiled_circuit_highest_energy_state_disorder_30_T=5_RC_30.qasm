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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(0.80817428) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(-1.5789403) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8573541) q[0];
sx q[0];
rz(-2.1942217) q[0];
sx q[0];
rz(-2.4592072) q[0];
rz(-pi) q[1];
rz(1.8380661) q[2];
sx q[2];
rz(-1.7198945) q[2];
sx q[2];
rz(1.9097415) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6526323) q[1];
sx q[1];
rz(-2.6217786) q[1];
sx q[1];
rz(1.7955154) q[1];
rz(1.6049625) q[3];
sx q[3];
rz(-1.6164268) q[3];
sx q[3];
rz(0.51866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.943104) q[2];
sx q[2];
rz(-0.21185943) q[2];
sx q[2];
rz(-2.1073821) q[2];
rz(2.4487623) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(1.0606162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387872) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(0.57648188) q[0];
rz(-3.1209962) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(2.1131262) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8688558) q[0];
sx q[0];
rz(-1.3139561) q[0];
sx q[0];
rz(-1.7935852) q[0];
rz(-1.1972381) q[2];
sx q[2];
rz(-1.2837063) q[2];
sx q[2];
rz(0.32880983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8081144) q[1];
sx q[1];
rz(-1.2385744) q[1];
sx q[1];
rz(-1.9971041) q[1];
rz(-1.1787492) q[3];
sx q[3];
rz(-0.96091753) q[3];
sx q[3];
rz(1.4489685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2241406) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(1.2646328) q[2];
rz(-2.1680016) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(2.8784331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6384386) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(-0.85357443) q[0];
rz(1.0517993) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(-0.015017088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029808345) q[0];
sx q[0];
rz(-1.5795603) q[0];
sx q[0];
rz(-0.86017227) q[0];
x q[1];
rz(0.21785801) q[2];
sx q[2];
rz(-0.098420489) q[2];
sx q[2];
rz(-0.86505666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5530921) q[1];
sx q[1];
rz(-0.71770826) q[1];
sx q[1];
rz(-0.11228447) q[1];
rz(-pi) q[2];
rz(-1.6347359) q[3];
sx q[3];
rz(-1.8600762) q[3];
sx q[3];
rz(2.6562986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0217648) q[2];
sx q[2];
rz(-1.4889577) q[2];
sx q[2];
rz(0.28142288) q[2];
rz(-1.4540539) q[3];
sx q[3];
rz(-1.2097996) q[3];
sx q[3];
rz(0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948831) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(2.967714) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(1.4283659) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9857432) q[0];
sx q[0];
rz(-1.5850889) q[0];
sx q[0];
rz(1.3784268) q[0];
x q[1];
rz(-2.8160353) q[2];
sx q[2];
rz(-0.37806219) q[2];
sx q[2];
rz(-2.4488673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97221365) q[1];
sx q[1];
rz(-1.5323109) q[1];
sx q[1];
rz(-2.4840218) q[1];
rz(0.81528109) q[3];
sx q[3];
rz(-0.51183701) q[3];
sx q[3];
rz(2.3533604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7017158) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-2.7244549) q[2];
rz(0.16962984) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(0.72181845) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193264) q[0];
sx q[0];
rz(-2.7665783) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(-2.3720692) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(1.5187029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5734539) q[0];
sx q[0];
rz(-0.72069695) q[0];
sx q[0];
rz(-1.9573523) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74414545) q[2];
sx q[2];
rz(-2.2501906) q[2];
sx q[2];
rz(1.3199922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7243821) q[1];
sx q[1];
rz(-1.3816178) q[1];
sx q[1];
rz(-0.030102109) q[1];
rz(-pi) q[2];
rz(1.4368016) q[3];
sx q[3];
rz(-1.6686642) q[3];
sx q[3];
rz(1.1214453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7421444) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(-2.4230912) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(0.79513454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28498483) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(1.5307776) q[0];
rz(-1.4467622) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(1.4379427) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7890678) q[0];
sx q[0];
rz(-1.97856) q[0];
sx q[0];
rz(-3.1192794) q[0];
rz(-pi) q[1];
rz(-2.6126768) q[2];
sx q[2];
rz(-2.4806266) q[2];
sx q[2];
rz(1.5691552) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6561778) q[1];
sx q[1];
rz(-2.6256769) q[1];
sx q[1];
rz(-2.9713216) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75957652) q[3];
sx q[3];
rz(-1.5686706) q[3];
sx q[3];
rz(1.1689982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1630254) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(-1.9095518) q[2];
rz(1.1466522) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(-0.31057096) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6335886) q[0];
sx q[0];
rz(-2.3242943) q[0];
sx q[0];
rz(2.001413) q[0];
rz(-0.89705244) q[1];
sx q[1];
rz(-1.8042118) q[1];
sx q[1];
rz(-0.030489771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.545118) q[0];
sx q[0];
rz(-1.4545355) q[0];
sx q[0];
rz(-1.7932939) q[0];
x q[1];
rz(-0.42564904) q[2];
sx q[2];
rz(-2.4535745) q[2];
sx q[2];
rz(-1.4378215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6086238) q[1];
sx q[1];
rz(-2.2734959) q[1];
sx q[1];
rz(-0.38680248) q[1];
rz(2.3014841) q[3];
sx q[3];
rz(-1.2222154) q[3];
sx q[3];
rz(2.9510562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19411479) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-2.1066809) q[3];
sx q[3];
rz(1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3461935) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.7561703) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.4477891) q[1];
sx q[1];
rz(0.48929712) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74936825) q[0];
sx q[0];
rz(-2.4697692) q[0];
sx q[0];
rz(0.031731204) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6750181) q[2];
sx q[2];
rz(-0.33793338) q[2];
sx q[2];
rz(0.80160415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41114901) q[1];
sx q[1];
rz(-0.48161067) q[1];
sx q[1];
rz(1.7754619) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0396949) q[3];
sx q[3];
rz(-1.3227665) q[3];
sx q[3];
rz(2.0534648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74799246) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(2.6753814) q[2];
rz(-1.6097869) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(-2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521249) q[0];
sx q[0];
rz(-0.89972275) q[0];
sx q[0];
rz(0.049064431) q[0];
rz(-0.24066726) q[1];
sx q[1];
rz(-2.1235695) q[1];
sx q[1];
rz(-3.1191471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81240679) q[0];
sx q[0];
rz(-0.57984771) q[0];
sx q[0];
rz(1.5270698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9188377) q[2];
sx q[2];
rz(-1.204426) q[2];
sx q[2];
rz(2.6719246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7514403) q[1];
sx q[1];
rz(-0.94822394) q[1];
sx q[1];
rz(2.1639813) q[1];
x q[2];
rz(0.53655728) q[3];
sx q[3];
rz(-2.381455) q[3];
sx q[3];
rz(0.11790568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(0.9683041) q[2];
rz(-2.3051895) q[3];
sx q[3];
rz(-0.29763779) q[3];
sx q[3];
rz(-1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4150998) q[0];
sx q[0];
rz(-1.2905755) q[0];
sx q[0];
rz(0.37877628) q[0];
rz(3.1355766) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(-2.5078497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7597719) q[0];
sx q[0];
rz(-2.2359747) q[0];
sx q[0];
rz(1.1367255) q[0];
rz(-pi) q[1];
x q[1];
rz(1.815755) q[2];
sx q[2];
rz(-2.5542521) q[2];
sx q[2];
rz(0.4191242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56713518) q[1];
sx q[1];
rz(-2.0525041) q[1];
sx q[1];
rz(2.5888799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59468693) q[3];
sx q[3];
rz(-2.1592038) q[3];
sx q[3];
rz(0.4247409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70574957) q[2];
sx q[2];
rz(-1.6195932) q[2];
sx q[2];
rz(2.5283234) q[2];
rz(-2.4663726) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(2.3238382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9153862) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(0.040493852) q[2];
sx q[2];
rz(-2.6815985) q[2];
sx q[2];
rz(0.56618377) q[2];
rz(-3.0499115) q[3];
sx q[3];
rz(-2.1909919) q[3];
sx q[3];
rz(1.7237678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
