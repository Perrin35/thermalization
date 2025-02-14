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
rz(-2.3174602) q[0];
sx q[0];
rz(-1.1622518) q[0];
sx q[0];
rz(-0.21221575) q[0];
rz(-2.4556887) q[1];
sx q[1];
rz(-0.75690126) q[1];
sx q[1];
rz(-0.25605717) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4418418) q[0];
sx q[0];
rz(-2.3810593) q[0];
sx q[0];
rz(-0.89450128) q[0];
rz(-3.0546064) q[2];
sx q[2];
rz(-2.1300007) q[2];
sx q[2];
rz(-2.5285697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.040671445) q[1];
sx q[1];
rz(-1.36449) q[1];
sx q[1];
rz(1.7357769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1027378) q[3];
sx q[3];
rz(-1.3016455) q[3];
sx q[3];
rz(-2.0520963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2519007) q[2];
sx q[2];
rz(-2.8163741) q[2];
sx q[2];
rz(-1.088781) q[2];
rz(2.661656) q[3];
sx q[3];
rz(-1.1705385) q[3];
sx q[3];
rz(1.1133194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29555175) q[0];
sx q[0];
rz(-0.22838455) q[0];
sx q[0];
rz(-2.8725655) q[0];
rz(0.61664063) q[1];
sx q[1];
rz(-2.4863305) q[1];
sx q[1];
rz(0.90589595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038074) q[0];
sx q[0];
rz(-1.4965987) q[0];
sx q[0];
rz(-0.78671771) q[0];
x q[1];
rz(2.3712641) q[2];
sx q[2];
rz(-0.69122756) q[2];
sx q[2];
rz(0.28797418) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62123728) q[1];
sx q[1];
rz(-1.3150421) q[1];
sx q[1];
rz(1.4489555) q[1];
rz(-pi) q[2];
rz(-1.18371) q[3];
sx q[3];
rz(-1.0446232) q[3];
sx q[3];
rz(-2.5691606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1632094) q[2];
sx q[2];
rz(-2.3613112) q[2];
sx q[2];
rz(-0.015446375) q[2];
rz(-0.77203006) q[3];
sx q[3];
rz(-2.6431712) q[3];
sx q[3];
rz(1.8933403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.5323199) q[0];
sx q[0];
rz(-2.5072704) q[0];
sx q[0];
rz(2.3442205) q[0];
rz(-2.4728921) q[1];
sx q[1];
rz(-1.0182321) q[1];
sx q[1];
rz(-0.55577898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3674372) q[0];
sx q[0];
rz(-1.3688911) q[0];
sx q[0];
rz(0.76011472) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3923116) q[2];
sx q[2];
rz(-0.67214078) q[2];
sx q[2];
rz(1.3591421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2712471) q[1];
sx q[1];
rz(-1.5588798) q[1];
sx q[1];
rz(-2.4153277) q[1];
x q[2];
rz(-1.9706412) q[3];
sx q[3];
rz(-2.2388864) q[3];
sx q[3];
rz(0.081351697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7734311) q[2];
sx q[2];
rz(-1.975446) q[2];
sx q[2];
rz(-1.2097166) q[2];
rz(-0.11956735) q[3];
sx q[3];
rz(-0.60496324) q[3];
sx q[3];
rz(-2.1623478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6794353) q[0];
sx q[0];
rz(-1.4428416) q[0];
sx q[0];
rz(-0.47732842) q[0];
rz(3.0408472) q[1];
sx q[1];
rz(-1.0557749) q[1];
sx q[1];
rz(-3.005952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076178032) q[0];
sx q[0];
rz(-1.9913902) q[0];
sx q[0];
rz(2.2432144) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56265752) q[2];
sx q[2];
rz(-0.51501319) q[2];
sx q[2];
rz(-0.80660179) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6203719) q[1];
sx q[1];
rz(-1.4209827) q[1];
sx q[1];
rz(0.42370875) q[1];
rz(-pi) q[2];
rz(-2.3587386) q[3];
sx q[3];
rz(-1.4042036) q[3];
sx q[3];
rz(-1.5448017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2123083) q[2];
sx q[2];
rz(-1.6529275) q[2];
sx q[2];
rz(-0.51367122) q[2];
rz(2.9786927) q[3];
sx q[3];
rz(-2.8859911) q[3];
sx q[3];
rz(2.0388849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61360079) q[0];
sx q[0];
rz(-3.1311212) q[0];
sx q[0];
rz(0.42798671) q[0];
rz(-1.1780659) q[1];
sx q[1];
rz(-1.4617498) q[1];
sx q[1];
rz(-1.4129432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75668272) q[0];
sx q[0];
rz(-1.9506978) q[0];
sx q[0];
rz(-1.0713473) q[0];
x q[1];
rz(-1.6128849) q[2];
sx q[2];
rz(-1.5344991) q[2];
sx q[2];
rz(2.5840379) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.481491) q[1];
sx q[1];
rz(-0.30687422) q[1];
sx q[1];
rz(2.9576388) q[1];
x q[2];
rz(-1.1023476) q[3];
sx q[3];
rz(-2.0878359) q[3];
sx q[3];
rz(-0.11786961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.596375) q[2];
sx q[2];
rz(-1.7866106) q[2];
sx q[2];
rz(-0.4893111) q[2];
rz(0.69928586) q[3];
sx q[3];
rz(-1.0438865) q[3];
sx q[3];
rz(-1.6635241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3938703) q[0];
sx q[0];
rz(-1.007653) q[0];
sx q[0];
rz(1.7973768) q[0];
rz(-2.4463553) q[1];
sx q[1];
rz(-1.2691011) q[1];
sx q[1];
rz(-0.42323798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7376331) q[0];
sx q[0];
rz(-2.9610017) q[0];
sx q[0];
rz(2.7904912) q[0];
rz(2.3219522) q[2];
sx q[2];
rz(-1.4965349) q[2];
sx q[2];
rz(-0.008645388) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5343439) q[1];
sx q[1];
rz(-2.3869276) q[1];
sx q[1];
rz(-0.21790803) q[1];
rz(2.845351) q[3];
sx q[3];
rz(-2.9518173) q[3];
sx q[3];
rz(1.1452183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3678579) q[2];
sx q[2];
rz(-2.9025142) q[2];
sx q[2];
rz(-1.9284922) q[2];
rz(0.086094543) q[3];
sx q[3];
rz(-1.2505069) q[3];
sx q[3];
rz(-1.7566173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64751476) q[0];
sx q[0];
rz(-0.32951847) q[0];
sx q[0];
rz(-2.9373017) q[0];
rz(0.34945166) q[1];
sx q[1];
rz(-1.5705669) q[1];
sx q[1];
rz(1.516516) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1384927) q[0];
sx q[0];
rz(-0.47687832) q[0];
sx q[0];
rz(1.80507) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2247844) q[2];
sx q[2];
rz(-2.6497239) q[2];
sx q[2];
rz(-2.7654344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2353846) q[1];
sx q[1];
rz(-1.5245807) q[1];
sx q[1];
rz(-1.4704204) q[1];
rz(-pi) q[2];
rz(1.54957) q[3];
sx q[3];
rz(-0.57820613) q[3];
sx q[3];
rz(2.8631722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1427631) q[2];
sx q[2];
rz(-2.5122354) q[2];
sx q[2];
rz(0.91267419) q[2];
rz(3.1257889) q[3];
sx q[3];
rz(-1.088257) q[3];
sx q[3];
rz(0.44939941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.7939746) q[0];
sx q[0];
rz(-2.0517218) q[0];
sx q[0];
rz(2.3225978) q[0];
rz(-2.7836169) q[1];
sx q[1];
rz(-2.0134108) q[1];
sx q[1];
rz(-2.5111759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25787658) q[0];
sx q[0];
rz(-1.64332) q[0];
sx q[0];
rz(0.016091165) q[0];
rz(-pi) q[1];
rz(-2.0840896) q[2];
sx q[2];
rz(-2.3529511) q[2];
sx q[2];
rz(0.12257931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.158327) q[1];
sx q[1];
rz(-0.7559146) q[1];
sx q[1];
rz(2.1793038) q[1];
x q[2];
rz(0.96004769) q[3];
sx q[3];
rz(-2.4195101) q[3];
sx q[3];
rz(-0.86632767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7716052) q[2];
sx q[2];
rz(-2.6321754) q[2];
sx q[2];
rz(0.39311692) q[2];
rz(-2.3734132) q[3];
sx q[3];
rz(-2.3438175) q[3];
sx q[3];
rz(1.2099077) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2792252) q[0];
sx q[0];
rz(-0.51280642) q[0];
sx q[0];
rz(-0.38238907) q[0];
rz(-0.44876107) q[1];
sx q[1];
rz(-0.02350137) q[1];
sx q[1];
rz(0.99865595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6210839) q[0];
sx q[0];
rz(-1.3338085) q[0];
sx q[0];
rz(-0.67279664) q[0];
rz(1.8085747) q[2];
sx q[2];
rz(-0.3467803) q[2];
sx q[2];
rz(1.717076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0599715) q[1];
sx q[1];
rz(-1.8619027) q[1];
sx q[1];
rz(-1.3636834) q[1];
x q[2];
rz(1.5019997) q[3];
sx q[3];
rz(-0.88310396) q[3];
sx q[3];
rz(-2.6042134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0465595) q[2];
sx q[2];
rz(-0.79594374) q[2];
sx q[2];
rz(2.156192) q[2];
rz(0.69917786) q[3];
sx q[3];
rz(-0.85088426) q[3];
sx q[3];
rz(3.0881171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.2153274) q[0];
sx q[0];
rz(-0.43020058) q[0];
sx q[0];
rz(0.77558023) q[0];
rz(-0.27014488) q[1];
sx q[1];
rz(-1.4560478) q[1];
sx q[1];
rz(2.5712579) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4264589) q[0];
sx q[0];
rz(-2.0394148) q[0];
sx q[0];
rz(-0.56632407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.774046) q[2];
sx q[2];
rz(-2.654944) q[2];
sx q[2];
rz(-1.987285) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3667222) q[1];
sx q[1];
rz(-1.6022283) q[1];
sx q[1];
rz(2.9079939) q[1];
x q[2];
rz(3.0489286) q[3];
sx q[3];
rz(-0.69159019) q[3];
sx q[3];
rz(-3.0244477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21738416) q[2];
sx q[2];
rz(-1.1576098) q[2];
sx q[2];
rz(0.50979924) q[2];
rz(2.9057251) q[3];
sx q[3];
rz(-2.9678952) q[3];
sx q[3];
rz(-0.97714669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3167284) q[0];
sx q[0];
rz(-1.4126128) q[0];
sx q[0];
rz(1.8870507) q[0];
rz(0.53312373) q[1];
sx q[1];
rz(-1.6785379) q[1];
sx q[1];
rz(1.0482845) q[1];
rz(3.127373) q[2];
sx q[2];
rz(-1.996025) q[2];
sx q[2];
rz(1.7964515) q[2];
rz(0.20797603) q[3];
sx q[3];
rz(-2.6317876) q[3];
sx q[3];
rz(-0.77710487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
