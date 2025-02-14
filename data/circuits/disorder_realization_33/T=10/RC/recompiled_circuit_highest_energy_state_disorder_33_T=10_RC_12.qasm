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
rz(3.7137346) q[0];
sx q[0];
rz(9.619286) q[0];
rz(2.1394849) q[1];
sx q[1];
rz(3.9245457) q[1];
sx q[1];
rz(9.4406162) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7730071) q[0];
sx q[0];
rz(-1.1601279) q[0];
sx q[0];
rz(-2.3334614) q[0];
rz(-pi) q[1];
rz(-0.39762605) q[2];
sx q[2];
rz(-2.3038452) q[2];
sx q[2];
rz(-0.83329569) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8037607) q[1];
sx q[1];
rz(-1.4678218) q[1];
sx q[1];
rz(-2.2096358) q[1];
rz(-pi) q[2];
rz(1.7884618) q[3];
sx q[3];
rz(-2.7870745) q[3];
sx q[3];
rz(-0.61787546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1281434) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(2.5246998) q[2];
rz(1.8938176) q[3];
sx q[3];
rz(-2.3168677) q[3];
sx q[3];
rz(-2.8491546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040267471) q[0];
sx q[0];
rz(-1.1722925) q[0];
sx q[0];
rz(3.110926) q[0];
rz(-1.6568294) q[1];
sx q[1];
rz(-0.29850423) q[1];
sx q[1];
rz(-0.31918496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46542612) q[0];
sx q[0];
rz(-2.5990708) q[0];
sx q[0];
rz(1.2529621) q[0];
rz(-1.1775374) q[2];
sx q[2];
rz(-1.3280091) q[2];
sx q[2];
rz(0.78906203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7564076) q[1];
sx q[1];
rz(-1.1435724) q[1];
sx q[1];
rz(0.34308497) q[1];
rz(0.52842657) q[3];
sx q[3];
rz(-2.0459262) q[3];
sx q[3];
rz(1.3412061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9926051) q[2];
sx q[2];
rz(-2.7379898) q[2];
sx q[2];
rz(0.65275025) q[2];
rz(0.8695237) q[3];
sx q[3];
rz(-2.0448989) q[3];
sx q[3];
rz(-2.7147527) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2000548) q[0];
sx q[0];
rz(-2.8917199) q[0];
sx q[0];
rz(0.13724929) q[0];
rz(-2.6380154) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(0.30127475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190126) q[0];
sx q[0];
rz(-2.2323881) q[0];
sx q[0];
rz(1.9674752) q[0];
rz(1.3728447) q[2];
sx q[2];
rz(-2.8471409) q[2];
sx q[2];
rz(-2.1535923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91567809) q[1];
sx q[1];
rz(-1.7748194) q[1];
sx q[1];
rz(1.7197439) q[1];
rz(-1.4120164) q[3];
sx q[3];
rz(-2.232312) q[3];
sx q[3];
rz(-0.56453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0425903) q[2];
sx q[2];
rz(-2.4799535) q[2];
sx q[2];
rz(0.46690565) q[2];
rz(0.0070988797) q[3];
sx q[3];
rz(-0.10741281) q[3];
sx q[3];
rz(-0.98889178) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0840833) q[0];
sx q[0];
rz(-0.0060265344) q[0];
sx q[0];
rz(0.65473336) q[0];
rz(1.5714802) q[1];
sx q[1];
rz(-1.4288158) q[1];
sx q[1];
rz(2.7007801) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2991371) q[0];
sx q[0];
rz(-1.5464325) q[0];
sx q[0];
rz(-3.117119) q[0];
rz(3.12241) q[2];
sx q[2];
rz(-2.6628116) q[2];
sx q[2];
rz(0.60730241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.310439) q[1];
sx q[1];
rz(-2.7686627) q[1];
sx q[1];
rz(2.4937477) q[1];
rz(-pi) q[2];
rz(1.8809404) q[3];
sx q[3];
rz(-2.0278483) q[3];
sx q[3];
rz(-0.09121212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4948027) q[2];
sx q[2];
rz(-1.9014634) q[2];
sx q[2];
rz(-0.49171641) q[2];
rz(-1.3382737) q[3];
sx q[3];
rz(-2.0526142) q[3];
sx q[3];
rz(1.4664388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717475) q[0];
sx q[0];
rz(-1.0862229) q[0];
sx q[0];
rz(1.0787429) q[0];
rz(3.0951913) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(-2.0414415) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85202269) q[0];
sx q[0];
rz(-1.6873813) q[0];
sx q[0];
rz(-1.5800716) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2583986) q[2];
sx q[2];
rz(-2.1015169) q[2];
sx q[2];
rz(-3.137756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55441862) q[1];
sx q[1];
rz(-2.3918416) q[1];
sx q[1];
rz(-2.7778704) q[1];
x q[2];
rz(1.6836105) q[3];
sx q[3];
rz(-1.8387018) q[3];
sx q[3];
rz(-0.15427854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2940353) q[2];
sx q[2];
rz(-0.39176771) q[2];
sx q[2];
rz(-2.7993287) q[2];
rz(-1.3058454) q[3];
sx q[3];
rz(-1.9394453) q[3];
sx q[3];
rz(-2.0391298) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24868988) q[0];
sx q[0];
rz(-1.9547639) q[0];
sx q[0];
rz(2.7500395) q[0];
rz(0.64741099) q[1];
sx q[1];
rz(-0.4386816) q[1];
sx q[1];
rz(2.2364) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86588106) q[0];
sx q[0];
rz(-2.0581332) q[0];
sx q[0];
rz(0.62222991) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2268385) q[2];
sx q[2];
rz(-1.1335271) q[2];
sx q[2];
rz(0.8503051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4005112) q[1];
sx q[1];
rz(-1.3071187) q[1];
sx q[1];
rz(2.5225894) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8531032) q[3];
sx q[3];
rz(-1.6181711) q[3];
sx q[3];
rz(2.4758439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68760005) q[2];
sx q[2];
rz(-2.8314721) q[2];
sx q[2];
rz(-2.0886759) q[2];
rz(0.29371253) q[3];
sx q[3];
rz(-0.83455825) q[3];
sx q[3];
rz(-2.6101051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76327908) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(0.026750201) q[0];
rz(1.9552975) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(0.15348405) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1790548) q[0];
sx q[0];
rz(-1.0729331) q[0];
sx q[0];
rz(2.8341917) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6233968) q[2];
sx q[2];
rz(-1.0211049) q[2];
sx q[2];
rz(-0.48940966) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3378422) q[1];
sx q[1];
rz(-2.1598247) q[1];
sx q[1];
rz(-2.4250125) q[1];
rz(-pi) q[2];
rz(1.1304705) q[3];
sx q[3];
rz(-0.067256602) q[3];
sx q[3];
rz(-2.7471144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89580506) q[2];
sx q[2];
rz(-0.20496002) q[2];
sx q[2];
rz(-1.8765571) q[2];
rz(0.17299077) q[3];
sx q[3];
rz(-1.3909307) q[3];
sx q[3];
rz(-2.0632108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0377401) q[0];
sx q[0];
rz(-2.9321892) q[0];
sx q[0];
rz(0.10845342) q[0];
rz(1.7697889) q[1];
sx q[1];
rz(-1.7864952) q[1];
sx q[1];
rz(3.1350737) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20978808) q[0];
sx q[0];
rz(-0.75388718) q[0];
sx q[0];
rz(0.86236282) q[0];
rz(-pi) q[1];
rz(1.4638607) q[2];
sx q[2];
rz(-1.3513913) q[2];
sx q[2];
rz(2.1768513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7830917) q[1];
sx q[1];
rz(-0.83493646) q[1];
sx q[1];
rz(2.3651644) q[1];
x q[2];
rz(-2.7192332) q[3];
sx q[3];
rz(-0.82882753) q[3];
sx q[3];
rz(-3.1008848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6846377) q[2];
sx q[2];
rz(-0.64720398) q[2];
sx q[2];
rz(-1.7015438) q[2];
rz(-2.7889934) q[3];
sx q[3];
rz(-2.2779901) q[3];
sx q[3];
rz(-0.45261228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088293485) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(1.1540867) q[0];
rz(1.6905009) q[1];
sx q[1];
rz(-2.4595551) q[1];
sx q[1];
rz(-2.3114204) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367171) q[0];
sx q[0];
rz(-1.5815539) q[0];
sx q[0];
rz(-1.5675926) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24014704) q[2];
sx q[2];
rz(-0.51371759) q[2];
sx q[2];
rz(1.3432389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3026003) q[1];
sx q[1];
rz(-0.60878372) q[1];
sx q[1];
rz(-0.7418315) q[1];
rz(1.8228047) q[3];
sx q[3];
rz(-0.67630905) q[3];
sx q[3];
rz(-2.0700434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95149583) q[2];
sx q[2];
rz(-0.87679619) q[2];
sx q[2];
rz(-2.4693176) q[2];
rz(0.43664524) q[3];
sx q[3];
rz(-1.3029706) q[3];
sx q[3];
rz(-0.28948998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.745568) q[0];
sx q[0];
rz(-0.73306274) q[0];
sx q[0];
rz(0.42098862) q[0];
rz(1.1517395) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(2.4236325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508681) q[0];
sx q[0];
rz(-1.2932672) q[0];
sx q[0];
rz(0.023157774) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46109445) q[2];
sx q[2];
rz(-2.3873219) q[2];
sx q[2];
rz(-2.8476771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58484287) q[1];
sx q[1];
rz(-1.6989917) q[1];
sx q[1];
rz(-1.2570279) q[1];
rz(-0.27394343) q[3];
sx q[3];
rz(-1.1355601) q[3];
sx q[3];
rz(1.0268322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.026160508) q[2];
sx q[2];
rz(-0.60439503) q[2];
sx q[2];
rz(-1.3755414) q[2];
rz(-2.9653505) q[3];
sx q[3];
rz(-2.3459489) q[3];
sx q[3];
rz(-3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
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
rz(-0.86028987) q[1];
sx q[1];
rz(-1.6382891) q[1];
sx q[1];
rz(2.0232497) q[1];
rz(-0.064762887) q[2];
sx q[2];
rz(-2.3751866) q[2];
sx q[2];
rz(-2.3768718) q[2];
rz(-0.82127251) q[3];
sx q[3];
rz(-2.6340064) q[3];
sx q[3];
rz(-0.71024714) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
