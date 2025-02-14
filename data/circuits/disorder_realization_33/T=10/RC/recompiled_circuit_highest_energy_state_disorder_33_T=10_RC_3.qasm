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
rz(-3.1193982) q[0];
sx q[0];
rz(-0.57214195) q[0];
sx q[0];
rz(-0.194508) q[0];
rz(-1.0021078) q[1];
sx q[1];
rz(-0.78295308) q[1];
sx q[1];
rz(3.1257544) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793711) q[0];
sx q[0];
rz(-0.88484513) q[0];
sx q[0];
rz(-2.5995194) q[0];
rz(-pi) q[1];
rz(0.39762605) q[2];
sx q[2];
rz(-0.83774746) q[2];
sx q[2];
rz(-0.83329569) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3091649) q[1];
sx q[1];
rz(-2.2057071) q[1];
sx q[1];
rz(-3.0135703) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7884618) q[3];
sx q[3];
rz(-2.7870745) q[3];
sx q[3];
rz(2.5237172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.013449239) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(-2.5246998) q[2];
rz(-1.2477751) q[3];
sx q[3];
rz(-2.3168677) q[3];
sx q[3];
rz(-2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040267471) q[0];
sx q[0];
rz(-1.9693002) q[0];
sx q[0];
rz(0.030666703) q[0];
rz(-1.4847633) q[1];
sx q[1];
rz(-2.8430884) q[1];
sx q[1];
rz(-0.31918496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3094003) q[0];
sx q[0];
rz(-2.083389) q[0];
sx q[0];
rz(-0.18621791) q[0];
rz(-pi) q[1];
rz(1.1775374) q[2];
sx q[2];
rz(-1.8135836) q[2];
sx q[2];
rz(-2.3525306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0966677) q[1];
sx q[1];
rz(-2.6003692) q[1];
sx q[1];
rz(-0.9344395) q[1];
rz(2.6131661) q[3];
sx q[3];
rz(-1.0956665) q[3];
sx q[3];
rz(-1.8003866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9926051) q[2];
sx q[2];
rz(-0.4036029) q[2];
sx q[2];
rz(-0.65275025) q[2];
rz(0.8695237) q[3];
sx q[3];
rz(-1.0966938) q[3];
sx q[3];
rz(2.7147527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94153786) q[0];
sx q[0];
rz(-2.8917199) q[0];
sx q[0];
rz(-3.0043434) q[0];
rz(-0.50357729) q[1];
sx q[1];
rz(-0.63882393) q[1];
sx q[1];
rz(-2.8403179) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3225801) q[0];
sx q[0];
rz(-2.2323881) q[0];
sx q[0];
rz(-1.1741175) q[0];
x q[1];
rz(-0.059570407) q[2];
sx q[2];
rz(-1.859331) q[2];
sx q[2];
rz(2.3602006) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2259146) q[1];
sx q[1];
rz(-1.7748194) q[1];
sx q[1];
rz(1.7197439) q[1];
x q[2];
rz(1.4120164) q[3];
sx q[3];
rz(-2.232312) q[3];
sx q[3];
rz(0.56453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0990024) q[2];
sx q[2];
rz(-0.66163915) q[2];
sx q[2];
rz(0.46690565) q[2];
rz(3.1344938) q[3];
sx q[3];
rz(-0.10741281) q[3];
sx q[3];
rz(-2.1527009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0840833) q[0];
sx q[0];
rz(-0.0060265344) q[0];
sx q[0];
rz(2.4868593) q[0];
rz(-1.5701125) q[1];
sx q[1];
rz(-1.4288158) q[1];
sx q[1];
rz(-0.44081259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8705298) q[0];
sx q[0];
rz(-1.5952627) q[0];
sx q[0];
rz(1.5464252) q[0];
rz(-0.47870584) q[2];
sx q[2];
rz(-1.5619593) q[2];
sx q[2];
rz(0.98052006) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9928575) q[1];
sx q[1];
rz(-1.2760218) q[1];
sx q[1];
rz(1.8026505) q[1];
rz(-pi) q[2];
rz(-0.55542262) q[3];
sx q[3];
rz(-0.54612386) q[3];
sx q[3];
rz(2.6046942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4948027) q[2];
sx q[2];
rz(-1.2401293) q[2];
sx q[2];
rz(0.49171641) q[2];
rz(1.3382737) q[3];
sx q[3];
rz(-2.0526142) q[3];
sx q[3];
rz(-1.4664388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3717475) q[0];
sx q[0];
rz(-2.0553698) q[0];
sx q[0];
rz(-2.0628498) q[0];
rz(-0.046401333) q[1];
sx q[1];
rz(-1.5240144) q[1];
sx q[1];
rz(-1.1001512) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3691424) q[0];
sx q[0];
rz(-0.11695172) q[0];
sx q[0];
rz(-3.0625615) q[0];
x q[1];
rz(-2.4920667) q[2];
sx q[2];
rz(-2.1500271) q[2];
sx q[2];
rz(1.1729826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0667923) q[1];
sx q[1];
rz(-2.2612345) q[1];
sx q[1];
rz(1.2509173) q[1];
rz(1.4579822) q[3];
sx q[3];
rz(-1.3028909) q[3];
sx q[3];
rz(2.9873141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2940353) q[2];
sx q[2];
rz(-2.7498249) q[2];
sx q[2];
rz(2.7993287) q[2];
rz(-1.3058454) q[3];
sx q[3];
rz(-1.2021474) q[3];
sx q[3];
rz(2.0391298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24868988) q[0];
sx q[0];
rz(-1.1868287) q[0];
sx q[0];
rz(2.7500395) q[0];
rz(2.4941817) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(2.2364) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86588106) q[0];
sx q[0];
rz(-1.0834595) q[0];
sx q[0];
rz(2.5193627) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91475418) q[2];
sx q[2];
rz(-2.0080655) q[2];
sx q[2];
rz(0.8503051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1282922) q[1];
sx q[1];
rz(-0.97620539) q[1];
sx q[1];
rz(-1.250729) q[1];
rz(-0.16512434) q[3];
sx q[3];
rz(-0.29224631) q[3];
sx q[3];
rz(2.0783034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68760005) q[2];
sx q[2];
rz(-2.8314721) q[2];
sx q[2];
rz(-2.0886759) q[2];
rz(2.8478801) q[3];
sx q[3];
rz(-0.83455825) q[3];
sx q[3];
rz(-0.53148758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.3783136) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(-3.1148425) q[0];
rz(-1.9552975) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(2.9881086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492316) q[0];
sx q[0];
rz(-0.57825297) q[0];
sx q[0];
rz(-2.078767) q[0];
x q[1];
rz(-1.6233968) q[2];
sx q[2];
rz(-1.0211049) q[2];
sx q[2];
rz(2.652183) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2177362) q[1];
sx q[1];
rz(-2.1485275) q[1];
sx q[1];
rz(2.2958638) q[1];
rz(3.1128902) q[3];
sx q[3];
rz(-1.6316292) q[3];
sx q[3];
rz(0.83567747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2457876) q[2];
sx q[2];
rz(-2.9366326) q[2];
sx q[2];
rz(1.8765571) q[2];
rz(0.17299077) q[3];
sx q[3];
rz(-1.750662) q[3];
sx q[3];
rz(-1.0783819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0377401) q[0];
sx q[0];
rz(-2.9321892) q[0];
sx q[0];
rz(3.0331392) q[0];
rz(1.3718038) q[1];
sx q[1];
rz(-1.3550974) q[1];
sx q[1];
rz(3.1350737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.485584) q[0];
sx q[0];
rz(-2.1173899) q[0];
sx q[0];
rz(0.54838108) q[0];
rz(-pi) q[1];
rz(-1.677732) q[2];
sx q[2];
rz(-1.7902014) q[2];
sx q[2];
rz(-2.1768513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35850098) q[1];
sx q[1];
rz(-2.3066562) q[1];
sx q[1];
rz(-0.77642828) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1503073) q[3];
sx q[3];
rz(-2.3081927) q[3];
sx q[3];
rz(0.54617907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45695496) q[2];
sx q[2];
rz(-2.4943887) q[2];
sx q[2];
rz(1.7015438) q[2];
rz(2.7889934) q[3];
sx q[3];
rz(-0.86360258) q[3];
sx q[3];
rz(2.6889804) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088293485) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(1.9875059) q[0];
rz(1.4510918) q[1];
sx q[1];
rz(-0.68203753) q[1];
sx q[1];
rz(0.83017224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80487554) q[0];
sx q[0];
rz(-1.5815539) q[0];
sx q[0];
rz(-1.5740001) q[0];
rz(0.50135055) q[2];
sx q[2];
rz(-1.4536469) q[2];
sx q[2];
rz(-0.43767489) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3026003) q[1];
sx q[1];
rz(-0.60878372) q[1];
sx q[1];
rz(-2.3997612) q[1];
rz(1.8228047) q[3];
sx q[3];
rz(-2.4652836) q[3];
sx q[3];
rz(-1.0715493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1900968) q[2];
sx q[2];
rz(-2.2647965) q[2];
sx q[2];
rz(0.67227501) q[2];
rz(2.7049474) q[3];
sx q[3];
rz(-1.3029706) q[3];
sx q[3];
rz(0.28948998) q[3];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3960246) q[0];
sx q[0];
rz(-2.4085299) q[0];
sx q[0];
rz(0.42098862) q[0];
rz(-1.1517395) q[1];
sx q[1];
rz(-2.9321509) q[1];
sx q[1];
rz(-0.71796012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.46653) q[0];
sx q[0];
rz(-0.27846876) q[0];
sx q[0];
rz(1.4896926) q[0];
rz(0.69952632) q[2];
sx q[2];
rz(-1.2612064) q[2];
sx q[2];
rz(1.5173315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58484287) q[1];
sx q[1];
rz(-1.442601) q[1];
sx q[1];
rz(-1.2570279) q[1];
rz(-pi) q[2];
rz(1.1208543) q[3];
sx q[3];
rz(-1.8186343) q[3];
sx q[3];
rz(-0.42603394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1154321) q[2];
sx q[2];
rz(-2.5371976) q[2];
sx q[2];
rz(1.3755414) q[2];
rz(-2.9653505) q[3];
sx q[3];
rz(-2.3459489) q[3];
sx q[3];
rz(-3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3714704) q[0];
sx q[0];
rz(-1.6902516) q[0];
sx q[0];
rz(-1.1023735) q[0];
rz(0.86028987) q[1];
sx q[1];
rz(-1.5033036) q[1];
sx q[1];
rz(-1.1183429) q[1];
rz(3.0768298) q[2];
sx q[2];
rz(-2.3751866) q[2];
sx q[2];
rz(-2.3768718) q[2];
rz(1.9574439) q[3];
sx q[3];
rz(-1.9083229) q[3];
sx q[3];
rz(1.543386) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
