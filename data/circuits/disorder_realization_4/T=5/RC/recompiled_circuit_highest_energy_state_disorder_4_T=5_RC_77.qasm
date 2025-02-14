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
rz(0.29851222) q[0];
sx q[0];
rz(-1.8520344) q[0];
sx q[0];
rz(0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(7.0206849) q[1];
sx q[1];
rz(7.7590396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97853547) q[0];
sx q[0];
rz(-3.0594303) q[0];
sx q[0];
rz(-2.1294247) q[0];
rz(-pi) q[1];
rz(-1.9010004) q[2];
sx q[2];
rz(-1.6088952) q[2];
sx q[2];
rz(2.7862501) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45881328) q[1];
sx q[1];
rz(-0.87535697) q[1];
sx q[1];
rz(1.0473482) q[1];
rz(-pi) q[2];
rz(-2.6232417) q[3];
sx q[3];
rz(-2.1795594) q[3];
sx q[3];
rz(-1.4639554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6797356) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(-2.8409345) q[2];
rz(-0.65827185) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(-2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41475007) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-2.9657189) q[0];
rz(-2.580592) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(0.19634136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6628274) q[0];
sx q[0];
rz(-1.2365554) q[0];
sx q[0];
rz(-2.0408551) q[0];
x q[1];
rz(-0.11949338) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(1.3136065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2159816) q[1];
sx q[1];
rz(-0.51541939) q[1];
sx q[1];
rz(-1.5748596) q[1];
rz(-pi) q[2];
rz(-0.34336932) q[3];
sx q[3];
rz(-1.4447304) q[3];
sx q[3];
rz(0.32293636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.089513) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(-1.4698131) q[2];
rz(2.6164264) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(0.66450459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0199652) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(2.5110733) q[0];
rz(1.1446674) q[1];
sx q[1];
rz(-1.2455218) q[1];
sx q[1];
rz(-3.0847881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58724126) q[0];
sx q[0];
rz(-1.2040753) q[0];
sx q[0];
rz(1.1318867) q[0];
rz(-pi) q[1];
rz(-1.0089097) q[2];
sx q[2];
rz(-0.76912921) q[2];
sx q[2];
rz(3.1388856) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2395798) q[1];
sx q[1];
rz(-2.6551592) q[1];
sx q[1];
rz(-0.3063267) q[1];
rz(1.9980691) q[3];
sx q[3];
rz(-1.3766854) q[3];
sx q[3];
rz(2.8097092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28222617) q[2];
sx q[2];
rz(-1.6702009) q[2];
sx q[2];
rz(2.7024506) q[2];
rz(1.7450843) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(0.5784353) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(2.1266345) q[0];
rz(0.062648423) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(-0.28422022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361523) q[0];
sx q[0];
rz(-1.9198737) q[0];
sx q[0];
rz(0.45850584) q[0];
rz(-pi) q[1];
rz(-0.49782984) q[2];
sx q[2];
rz(-1.4477028) q[2];
sx q[2];
rz(0.13409889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7220347) q[1];
sx q[1];
rz(-2.7830208) q[1];
sx q[1];
rz(0.034073985) q[1];
x q[2];
rz(1.3174414) q[3];
sx q[3];
rz(-1.5572963) q[3];
sx q[3];
rz(2.6576192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3011498) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(-0.31989756) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(2.1130051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-2.5909162) q[0];
sx q[0];
rz(-0.40529761) q[0];
rz(-0.48794508) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(-0.36283666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3317446) q[0];
sx q[0];
rz(-0.54143006) q[0];
sx q[0];
rz(0.62744139) q[0];
x q[1];
rz(-1.1929197) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(-3.0966126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.664714) q[1];
sx q[1];
rz(-1.5221609) q[1];
sx q[1];
rz(-2.9205079) q[1];
rz(2.4461306) q[3];
sx q[3];
rz(-1.3331279) q[3];
sx q[3];
rz(2.9836751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.9826865) q[2];
sx q[2];
rz(-0.72251594) q[2];
rz(0.51269382) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(-1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(0.27477086) q[0];
rz(2.784506) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(1.1579317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5362893) q[0];
sx q[0];
rz(-0.9035631) q[0];
sx q[0];
rz(-0.31179223) q[0];
x q[1];
rz(2.544417) q[2];
sx q[2];
rz(-0.9998601) q[2];
sx q[2];
rz(1.2025637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5722583) q[1];
sx q[1];
rz(-2.5037994) q[1];
sx q[1];
rz(0.34697726) q[1];
rz(-pi) q[2];
rz(0.88648667) q[3];
sx q[3];
rz(-0.89596701) q[3];
sx q[3];
rz(1.0533028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-0.58970279) q[2];
rz(1.3127182) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(-0.5635128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24484672) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(0.27031159) q[0];
rz(2.7145794) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(2.7332773) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5184831) q[0];
sx q[0];
rz(-2.7528931) q[0];
sx q[0];
rz(2.8220909) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0403942) q[2];
sx q[2];
rz(-2.2658341) q[2];
sx q[2];
rz(-2.0027431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9198614) q[1];
sx q[1];
rz(-1.9834922) q[1];
sx q[1];
rz(0.92642118) q[1];
rz(-pi) q[2];
rz(-1.6736843) q[3];
sx q[3];
rz(-0.85082084) q[3];
sx q[3];
rz(2.306769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(2.2173524) q[2];
rz(-0.83485323) q[3];
sx q[3];
rz(-2.3122841) q[3];
sx q[3];
rz(-2.1448962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95558178) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(2.7741449) q[0];
rz(-0.5683178) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(-0.41499358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5552705) q[0];
sx q[0];
rz(-1.8316226) q[0];
sx q[0];
rz(-1.1010948) q[0];
x q[1];
rz(-0.52936036) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(2.8279378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69232443) q[1];
sx q[1];
rz(-1.0639166) q[1];
sx q[1];
rz(-2.436422) q[1];
rz(1.5563732) q[3];
sx q[3];
rz(-0.34581772) q[3];
sx q[3];
rz(-2.7399969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(0.83317327) q[2];
rz(2.7821275) q[3];
sx q[3];
rz(-1.0901674) q[3];
sx q[3];
rz(1.6370714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-3.1267082) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(0.10717779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2994605) q[0];
sx q[0];
rz(-2.3981574) q[0];
sx q[0];
rz(1.7920262) q[0];
rz(-pi) q[1];
x q[1];
rz(2.31865) q[2];
sx q[2];
rz(-1.9230033) q[2];
sx q[2];
rz(2.4885881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5276423) q[1];
sx q[1];
rz(-2.6384934) q[1];
sx q[1];
rz(-0.31219074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5507372) q[3];
sx q[3];
rz(-2.8300432) q[3];
sx q[3];
rz(0.24190608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(-2.2819819) q[2];
rz(-2.882242) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(-2.7249961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6421826) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.1580178) q[0];
rz(2.6462789) q[1];
sx q[1];
rz(-1.9751578) q[1];
sx q[1];
rz(0.29702979) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3619694) q[0];
sx q[0];
rz(-1.8898148) q[0];
sx q[0];
rz(-0.0042798288) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13003243) q[2];
sx q[2];
rz(-0.30270019) q[2];
sx q[2];
rz(-1.5720194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37849879) q[1];
sx q[1];
rz(-1.2796254) q[1];
sx q[1];
rz(-1.7353472) q[1];
x q[2];
rz(0.016344109) q[3];
sx q[3];
rz(-2.2333003) q[3];
sx q[3];
rz(1.7733044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(-0.86088172) q[2];
rz(-0.46368972) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(-1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1562445) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(1.8078177) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(2.4598224) q[2];
sx q[2];
rz(-0.78777704) q[2];
sx q[2];
rz(-2.4252979) q[2];
rz(1.9577338) q[3];
sx q[3];
rz(-2.0328184) q[3];
sx q[3];
rz(0.55749374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
