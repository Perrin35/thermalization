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
rz(0.022194447) q[0];
sx q[0];
rz(-2.5694507) q[0];
sx q[0];
rz(0.194508) q[0];
rz(2.1394849) q[1];
sx q[1];
rz(-2.3586396) q[1];
sx q[1];
rz(0.015838239) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7730071) q[0];
sx q[0];
rz(-1.9814648) q[0];
sx q[0];
rz(2.3334614) q[0];
rz(-pi) q[1];
rz(-0.79722793) q[2];
sx q[2];
rz(-1.2789032) q[2];
sx q[2];
rz(-0.46353093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33783198) q[1];
sx q[1];
rz(-1.4678218) q[1];
sx q[1];
rz(-0.93195685) q[1];
x q[2];
rz(-3.0618265) q[3];
sx q[3];
rz(-1.2249984) q[3];
sx q[3];
rz(2.7553204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.013449239) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(-2.5246998) q[2];
rz(1.8938176) q[3];
sx q[3];
rz(-0.82472491) q[3];
sx q[3];
rz(2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1013252) q[0];
sx q[0];
rz(-1.1722925) q[0];
sx q[0];
rz(-0.030666703) q[0];
rz(-1.4847633) q[1];
sx q[1];
rz(-0.29850423) q[1];
sx q[1];
rz(0.31918496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83219233) q[0];
sx q[0];
rz(-1.0582036) q[0];
sx q[0];
rz(0.18621791) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1445851) q[2];
sx q[2];
rz(-0.45882031) q[2];
sx q[2];
rz(-2.8850537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.038668949) q[1];
sx q[1];
rz(-1.2596411) q[1];
sx q[1];
rz(2.021108) q[1];
rz(-pi) q[2];
rz(-0.79547151) q[3];
sx q[3];
rz(-0.69517259) q[3];
sx q[3];
rz(2.2473587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9926051) q[2];
sx q[2];
rz(-2.7379898) q[2];
sx q[2];
rz(2.4888424) q[2];
rz(-2.2720689) q[3];
sx q[3];
rz(-2.0448989) q[3];
sx q[3];
rz(0.42683998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94153786) q[0];
sx q[0];
rz(-2.8917199) q[0];
sx q[0];
rz(0.13724929) q[0];
rz(-2.6380154) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(0.30127475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190126) q[0];
sx q[0];
rz(-2.2323881) q[0];
sx q[0];
rz(1.9674752) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8598156) q[2];
sx q[2];
rz(-1.6279015) q[2];
sx q[2];
rz(-0.77243519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62472407) q[1];
sx q[1];
rz(-1.4249603) q[1];
sx q[1];
rz(-0.20624344) q[1];
x q[2];
rz(2.4739315) q[3];
sx q[3];
rz(-1.4457089) q[3];
sx q[3];
rz(-1.1043105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0425903) q[2];
sx q[2];
rz(-2.4799535) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057509363) q[0];
sx q[0];
rz(-3.1355661) q[0];
sx q[0];
rz(0.65473336) q[0];
rz(1.5701125) q[1];
sx q[1];
rz(-1.4288158) q[1];
sx q[1];
rz(0.44081259) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0869319) q[0];
sx q[0];
rz(-3.107061) q[0];
sx q[0];
rz(-2.3582929) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.12241) q[2];
sx q[2];
rz(-2.6628116) q[2];
sx q[2];
rz(-0.60730241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.310439) q[1];
sx q[1];
rz(-0.37292994) q[1];
sx q[1];
rz(-0.64784494) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55542262) q[3];
sx q[3];
rz(-2.5954688) q[3];
sx q[3];
rz(0.53689843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4948027) q[2];
sx q[2];
rz(-1.2401293) q[2];
sx q[2];
rz(0.49171641) q[2];
rz(1.803319) q[3];
sx q[3];
rz(-1.0889784) q[3];
sx q[3];
rz(1.6751539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3717475) q[0];
sx q[0];
rz(-2.0553698) q[0];
sx q[0];
rz(1.0787429) q[0];
rz(3.0951913) q[1];
sx q[1];
rz(-1.5240144) q[1];
sx q[1];
rz(2.0414415) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71985258) q[0];
sx q[0];
rz(-1.5800086) q[0];
sx q[0];
rz(-0.11658998) q[0];
rz(-2.2583986) q[2];
sx q[2];
rz(-1.0400758) q[2];
sx q[2];
rz(0.0038366537) q[2];
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
x q[2];
rz(1.6836105) q[3];
sx q[3];
rz(-1.3028909) q[3];
sx q[3];
rz(0.15427854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2940353) q[2];
sx q[2];
rz(-2.7498249) q[2];
sx q[2];
rz(-2.7993287) q[2];
rz(1.3058454) q[3];
sx q[3];
rz(-1.2021474) q[3];
sx q[3];
rz(-2.0391298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24868988) q[0];
sx q[0];
rz(-1.1868287) q[0];
sx q[0];
rz(-0.3915531) q[0];
rz(0.64741099) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(-2.2364) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1126328) q[0];
sx q[0];
rz(-2.1117983) q[0];
sx q[0];
rz(-0.99287005) q[0];
x q[1];
rz(0.53296419) q[2];
sx q[2];
rz(-0.98532444) q[2];
sx q[2];
rz(1.0356071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47899095) q[1];
sx q[1];
rz(-2.4756065) q[1];
sx q[1];
rz(-2.7061092) q[1];
rz(-1.5213826) q[3];
sx q[3];
rz(-1.8589528) q[3];
sx q[3];
rz(0.8909944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4539926) q[2];
sx q[2];
rz(-2.8314721) q[2];
sx q[2];
rz(1.0529168) q[2];
rz(0.29371253) q[3];
sx q[3];
rz(-2.3070344) q[3];
sx q[3];
rz(2.6101051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783136) q[0];
sx q[0];
rz(-1.8623619) q[0];
sx q[0];
rz(3.1148425) q[0];
rz(1.1862952) q[1];
sx q[1];
rz(-2.9588638) q[1];
sx q[1];
rz(0.15348405) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5492316) q[0];
sx q[0];
rz(-0.57825297) q[0];
sx q[0];
rz(-2.078767) q[0];
x q[1];
rz(-2.5912845) q[2];
sx q[2];
rz(-1.5259503) q[2];
sx q[2];
rz(1.1088849) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3378422) q[1];
sx q[1];
rz(-0.98176793) q[1];
sx q[1];
rz(2.4250125) q[1];
rz(-pi) q[2];
rz(0.028702486) q[3];
sx q[3];
rz(-1.5099635) q[3];
sx q[3];
rz(-2.3059152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.89580506) q[2];
sx q[2];
rz(-2.9366326) q[2];
sx q[2];
rz(1.8765571) q[2];
rz(-2.9686019) q[3];
sx q[3];
rz(-1.750662) q[3];
sx q[3];
rz(2.0632108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10385253) q[0];
sx q[0];
rz(-0.2094035) q[0];
sx q[0];
rz(3.0331392) q[0];
rz(1.7697889) q[1];
sx q[1];
rz(-1.7864952) q[1];
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
rz(-1.0242027) q[0];
sx q[0];
rz(-0.54838108) q[0];
rz(-pi) q[1];
rz(0.22062515) q[2];
sx q[2];
rz(-1.4664337) q[2];
sx q[2];
rz(-2.5121784) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5287787) q[1];
sx q[1];
rz(-1.0135796) q[1];
sx q[1];
rz(2.2293836) q[1];
x q[2];
rz(0.78289276) q[3];
sx q[3];
rz(-1.2638348) q[3];
sx q[3];
rz(-1.8249194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45695496) q[2];
sx q[2];
rz(-0.64720398) q[2];
sx q[2];
rz(1.4400488) q[2];
rz(0.35259926) q[3];
sx q[3];
rz(-0.86360258) q[3];
sx q[3];
rz(-2.6889804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088293485) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(-1.9875059) q[0];
rz(1.6905009) q[1];
sx q[1];
rz(-0.68203753) q[1];
sx q[1];
rz(2.3114204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3756374) q[0];
sx q[0];
rz(-1.5739999) q[0];
sx q[0];
rz(0.010757627) q[0];
rz(-pi) q[1];
rz(1.7042034) q[2];
sx q[2];
rz(-2.0683943) q[2];
sx q[2];
rz(1.069151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7651319) q[1];
sx q[1];
rz(-1.9674976) q[1];
sx q[1];
rz(-2.6668616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8228047) q[3];
sx q[3];
rz(-0.67630905) q[3];
sx q[3];
rz(-1.0715493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95149583) q[2];
sx q[2];
rz(-0.87679619) q[2];
sx q[2];
rz(-2.4693176) q[2];
rz(2.7049474) q[3];
sx q[3];
rz(-1.3029706) q[3];
sx q[3];
rz(0.28948998) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.745568) q[0];
sx q[0];
rz(-0.73306274) q[0];
sx q[0];
rz(0.42098862) q[0];
rz(-1.9898532) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(-0.71796012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072453) q[0];
sx q[0];
rz(-1.8483254) q[0];
sx q[0];
rz(-0.023157774) q[0];
rz(-1.9667718) q[2];
sx q[2];
rz(-2.2309592) q[2];
sx q[2];
rz(0.3044314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1970994) q[1];
sx q[1];
rz(-1.8819032) q[1];
sx q[1];
rz(3.0068946) q[1];
x q[2];
rz(-2.8676492) q[3];
sx q[3];
rz(-1.1355601) q[3];
sx q[3];
rz(2.1147605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.026160508) q[2];
sx q[2];
rz(-0.60439503) q[2];
sx q[2];
rz(-1.7660512) q[2];
rz(2.9653505) q[3];
sx q[3];
rz(-2.3459489) q[3];
sx q[3];
rz(3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3714704) q[0];
sx q[0];
rz(-1.6902516) q[0];
sx q[0];
rz(-1.1023735) q[0];
rz(-2.2813028) q[1];
sx q[1];
rz(-1.5033036) q[1];
sx q[1];
rz(-1.1183429) q[1];
rz(-0.064762887) q[2];
sx q[2];
rz(-2.3751866) q[2];
sx q[2];
rz(-2.3768718) q[2];
rz(0.82127251) q[3];
sx q[3];
rz(-0.50758624) q[3];
sx q[3];
rz(2.4313455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
