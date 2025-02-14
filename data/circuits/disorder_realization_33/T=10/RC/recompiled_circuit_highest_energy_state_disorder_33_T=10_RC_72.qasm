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
rz(2.1394849) q[1];
sx q[1];
rz(-2.3586396) q[1];
sx q[1];
rz(0.015838239) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3685856) q[0];
sx q[0];
rz(-1.9814648) q[0];
sx q[0];
rz(0.80813129) q[0];
rz(0.39762605) q[2];
sx q[2];
rz(-2.3038452) q[2];
sx q[2];
rz(0.83329569) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8324278) q[1];
sx q[1];
rz(-0.93588557) q[1];
sx q[1];
rz(-0.12802235) q[1];
rz(-pi) q[2];
rz(1.7884618) q[3];
sx q[3];
rz(-0.35451815) q[3];
sx q[3];
rz(0.61787546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.013449239) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(0.61689287) q[2];
rz(-1.8938176) q[3];
sx q[3];
rz(-0.82472491) q[3];
sx q[3];
rz(-2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1013252) q[0];
sx q[0];
rz(-1.1722925) q[0];
sx q[0];
rz(3.110926) q[0];
rz(1.6568294) q[1];
sx q[1];
rz(-2.8430884) q[1];
sx q[1];
rz(-0.31918496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3108513) q[0];
sx q[0];
rz(-1.7328528) q[0];
sx q[0];
rz(1.0507163) q[0];
x q[1];
rz(1.1775374) q[2];
sx q[2];
rz(-1.3280091) q[2];
sx q[2];
rz(-0.78906203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.044925) q[1];
sx q[1];
rz(-0.54122347) q[1];
sx q[1];
rz(0.9344395) q[1];
x q[2];
rz(2.3461211) q[3];
sx q[3];
rz(-0.69517259) q[3];
sx q[3];
rz(-0.89423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14898758) q[2];
sx q[2];
rz(-0.4036029) q[2];
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
rz(-2.2000548) q[0];
sx q[0];
rz(-0.24987276) q[0];
sx q[0];
rz(-0.13724929) q[0];
rz(-2.6380154) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(-2.8403179) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8190126) q[0];
sx q[0];
rz(-0.9092046) q[0];
sx q[0];
rz(-1.1741175) q[0];
rz(1.281777) q[2];
sx q[2];
rz(-1.5136912) q[2];
sx q[2];
rz(0.77243519) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91567809) q[1];
sx q[1];
rz(-1.7748194) q[1];
sx q[1];
rz(1.4218487) q[1];
rz(-pi) q[2];
rz(-1.4120164) q[3];
sx q[3];
rz(-0.90928066) q[3];
sx q[3];
rz(0.56453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0990024) q[2];
sx q[2];
rz(-2.4799535) q[2];
sx q[2];
rz(-2.674687) q[2];
rz(-0.0070988797) q[3];
sx q[3];
rz(-0.10741281) q[3];
sx q[3];
rz(-2.1527009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-3.0840833) q[0];
sx q[0];
rz(-0.0060265344) q[0];
sx q[0];
rz(0.65473336) q[0];
rz(-1.5714802) q[1];
sx q[1];
rz(-1.7127769) q[1];
sx q[1];
rz(-0.44081259) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0869319) q[0];
sx q[0];
rz(-3.107061) q[0];
sx q[0];
rz(2.3582929) q[0];
x q[1];
rz(0.47870584) q[2];
sx q[2];
rz(-1.5796333) q[2];
sx q[2];
rz(-2.1610726) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9928575) q[1];
sx q[1];
rz(-1.2760218) q[1];
sx q[1];
rz(-1.8026505) q[1];
rz(-pi) q[2];
rz(2.58617) q[3];
sx q[3];
rz(-0.54612386) q[3];
sx q[3];
rz(-0.53689843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64678994) q[2];
sx q[2];
rz(-1.2401293) q[2];
sx q[2];
rz(-2.6498762) q[2];
rz(-1.803319) q[3];
sx q[3];
rz(-1.0889784) q[3];
sx q[3];
rz(-1.6751539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717475) q[0];
sx q[0];
rz(-2.0553698) q[0];
sx q[0];
rz(2.0628498) q[0];
rz(3.0951913) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(1.1001512) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4217401) q[0];
sx q[0];
rz(-1.561584) q[0];
sx q[0];
rz(0.11658998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3170839) q[2];
sx q[2];
rz(-0.84133278) q[2];
sx q[2];
rz(2.1192604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8535159) q[1];
sx q[1];
rz(-1.8156681) q[1];
sx q[1];
rz(2.4254786) q[1];
x q[2];
rz(-1.6836105) q[3];
sx q[3];
rz(-1.8387018) q[3];
sx q[3];
rz(-2.9873141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2940353) q[2];
sx q[2];
rz(-2.7498249) q[2];
sx q[2];
rz(2.7993287) q[2];
rz(-1.8357473) q[3];
sx q[3];
rz(-1.2021474) q[3];
sx q[3];
rz(-2.0391298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8929028) q[0];
sx q[0];
rz(-1.9547639) q[0];
sx q[0];
rz(0.3915531) q[0];
rz(-0.64741099) q[1];
sx q[1];
rz(-0.4386816) q[1];
sx q[1];
rz(-2.2364) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0150512) q[0];
sx q[0];
rz(-2.371696) q[0];
sx q[0];
rz(-2.4036744) q[0];
rz(-0.91692704) q[2];
sx q[2];
rz(-2.3715137) q[2];
sx q[2];
rz(-2.9240312) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6626017) q[1];
sx q[1];
rz(-0.66598611) q[1];
sx q[1];
rz(-2.7061092) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5213826) q[3];
sx q[3];
rz(-1.8589528) q[3];
sx q[3];
rz(-2.2505983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68760005) q[2];
sx q[2];
rz(-0.31012055) q[2];
sx q[2];
rz(-2.0886759) q[2];
rz(0.29371253) q[3];
sx q[3];
rz(-2.3070344) q[3];
sx q[3];
rz(2.6101051) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783136) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(0.026750201) q[0];
rz(1.9552975) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(-2.9881086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9625378) q[0];
sx q[0];
rz(-2.0686596) q[0];
sx q[0];
rz(0.30740096) q[0];
rz(-pi) q[1];
rz(0.085603733) q[2];
sx q[2];
rz(-0.55194469) q[2];
sx q[2];
rz(2.7526223) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19928923) q[1];
sx q[1];
rz(-2.2484712) q[1];
sx q[1];
rz(2.3476362) q[1];
x q[2];
rz(-1.5099385) q[3];
sx q[3];
rz(-1.5994457) q[3];
sx q[3];
rz(2.4047283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89580506) q[2];
sx q[2];
rz(-2.9366326) q[2];
sx q[2];
rz(1.2650355) q[2];
rz(0.17299077) q[3];
sx q[3];
rz(-1.3909307) q[3];
sx q[3];
rz(1.0783819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.10385253) q[0];
sx q[0];
rz(-0.2094035) q[0];
sx q[0];
rz(0.10845342) q[0];
rz(-1.3718038) q[1];
sx q[1];
rz(-1.3550974) q[1];
sx q[1];
rz(0.0065189204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9193503) q[0];
sx q[0];
rz(-1.1092289) q[0];
sx q[0];
rz(2.190175) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4638607) q[2];
sx q[2];
rz(-1.3513913) q[2];
sx q[2];
rz(0.96474136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5287787) q[1];
sx q[1];
rz(-1.0135796) q[1];
sx q[1];
rz(0.91220906) q[1];
rz(-1.9912854) q[3];
sx q[3];
rz(-0.83339993) q[3];
sx q[3];
rz(-0.54617907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45695496) q[2];
sx q[2];
rz(-2.4943887) q[2];
sx q[2];
rz(1.7015438) q[2];
rz(0.35259926) q[3];
sx q[3];
rz(-0.86360258) q[3];
sx q[3];
rz(-2.6889804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0532992) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(-1.1540867) q[0];
rz(-1.6905009) q[1];
sx q[1];
rz(-0.68203753) q[1];
sx q[1];
rz(-2.3114204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80487554) q[0];
sx q[0];
rz(-1.5600388) q[0];
sx q[0];
rz(1.5740001) q[0];
x q[1];
rz(-2.9014456) q[2];
sx q[2];
rz(-2.6278751) q[2];
sx q[2];
rz(1.7983537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37646078) q[1];
sx q[1];
rz(-1.9674976) q[1];
sx q[1];
rz(-0.47473102) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91009801) q[3];
sx q[3];
rz(-1.7275095) q[3];
sx q[3];
rz(-2.4441738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95149583) q[2];
sx q[2];
rz(-0.87679619) q[2];
sx q[2];
rz(2.4693176) q[2];
rz(-0.43664524) q[3];
sx q[3];
rz(-1.3029706) q[3];
sx q[3];
rz(0.28948998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.3960246) q[0];
sx q[0];
rz(-0.73306274) q[0];
sx q[0];
rz(-2.720604) q[0];
rz(1.9898532) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(0.71796012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97372598) q[0];
sx q[0];
rz(-1.5485248) q[0];
sx q[0];
rz(1.8483961) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46109445) q[2];
sx q[2];
rz(-0.75427079) q[2];
sx q[2];
rz(-0.2939156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7803256) q[1];
sx q[1];
rz(-2.8034489) q[1];
sx q[1];
rz(1.9664155) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0438536) q[3];
sx q[3];
rz(-0.5095616) q[3];
sx q[3];
rz(-1.614712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.026160508) q[2];
sx q[2];
rz(-0.60439503) q[2];
sx q[2];
rz(-1.7660512) q[2];
rz(-0.17624217) q[3];
sx q[3];
rz(-0.79564375) q[3];
sx q[3];
rz(-3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(1.7701223) q[0];
sx q[0];
rz(-1.4513411) q[0];
sx q[0];
rz(2.0392192) q[0];
rz(0.86028987) q[1];
sx q[1];
rz(-1.5033036) q[1];
sx q[1];
rz(-1.1183429) q[1];
rz(-1.5085718) q[2];
sx q[2];
rz(-2.3351861) q[2];
sx q[2];
rz(0.67493949) q[2];
rz(-1.1841487) q[3];
sx q[3];
rz(-1.9083229) q[3];
sx q[3];
rz(1.543386) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
