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
rz(0.98195568) q[1];
sx q[1];
rz(-0.73749956) q[1];
sx q[1];
rz(1.4758543) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1630572) q[0];
sx q[0];
rz(-0.082162372) q[0];
sx q[0];
rz(2.1294247) q[0];
x q[1];
rz(-1.453773) q[2];
sx q[2];
rz(-0.3323148) q[2];
sx q[2];
rz(-1.8154643) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6827794) q[1];
sx q[1];
rz(-0.87535697) q[1];
sx q[1];
rz(2.0942445) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51835097) q[3];
sx q[3];
rz(-2.1795594) q[3];
sx q[3];
rz(-1.6776372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6797356) q[2];
sx q[2];
rz(-1.4718141) q[2];
sx q[2];
rz(2.8409345) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-0.39977795) q[3];
sx q[3];
rz(0.58832735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41475007) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(2.9657189) q[0];
rz(-0.56100065) q[1];
sx q[1];
rz(-0.70204061) q[1];
sx q[1];
rz(2.9452513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6628274) q[0];
sx q[0];
rz(-1.9050373) q[0];
sx q[0];
rz(-1.1007376) q[0];
x q[1];
rz(-3.0220993) q[2];
sx q[2];
rz(-1.6867078) q[2];
sx q[2];
rz(-1.3136065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92561103) q[1];
sx q[1];
rz(-0.51541939) q[1];
sx q[1];
rz(1.566733) q[1];
x q[2];
rz(-1.4370055) q[3];
sx q[3];
rz(-1.9113298) q[3];
sx q[3];
rz(1.2029369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0520797) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(-2.6164264) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(0.66450459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0199652) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(2.5110733) q[0];
rz(-1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81679427) q[0];
sx q[0];
rz(-1.9787119) q[0];
sx q[0];
rz(2.7403031) q[0];
x q[1];
rz(2.6654451) q[2];
sx q[2];
rz(-2.2000929) q[2];
sx q[2];
rz(0.72222939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60384974) q[1];
sx q[1];
rz(-1.7122388) q[1];
sx q[1];
rz(2.6745922) q[1];
rz(-2.9288559) q[3];
sx q[3];
rz(-1.1520583) q[3];
sx q[3];
rz(1.9902844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28222617) q[2];
sx q[2];
rz(-1.6702009) q[2];
sx q[2];
rz(-0.43914208) q[2];
rz(-1.7450843) q[3];
sx q[3];
rz(-1.2625932) q[3];
sx q[3];
rz(0.5784353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948931) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(1.0149581) q[0];
rz(-3.0789442) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(-0.28422022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0089909) q[0];
sx q[0];
rz(-1.1418482) q[0];
sx q[0];
rz(-1.9563849) q[0];
rz(-pi) q[1];
rz(2.6437628) q[2];
sx q[2];
rz(-1.6938899) q[2];
sx q[2];
rz(-0.13409889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7220347) q[1];
sx q[1];
rz(-0.35857182) q[1];
sx q[1];
rz(-3.1075187) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.013945097) q[3];
sx q[3];
rz(-1.317465) q[3];
sx q[3];
rz(-1.0903181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3011498) q[2];
sx q[2];
rz(-1.0089077) q[2];
sx q[2];
rz(0.76048771) q[2];
rz(0.31989756) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29192057) q[0];
sx q[0];
rz(-2.5909162) q[0];
sx q[0];
rz(0.40529761) q[0];
rz(-2.6536476) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(-2.778756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3170119) q[0];
sx q[0];
rz(-1.2634227) q[0];
sx q[0];
rz(0.45305832) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5561781) q[2];
sx q[2];
rz(-1.8903132) q[2];
sx q[2];
rz(-1.8210757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47687864) q[1];
sx q[1];
rz(-1.6194317) q[1];
sx q[1];
rz(2.9205079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7801301) q[3];
sx q[3];
rz(-2.4130954) q[3];
sx q[3];
rz(2.0036774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.001658) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(-0.72251594) q[2];
rz(2.6288988) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672511) q[0];
sx q[0];
rz(-0.63684547) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(-0.35708669) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(-1.1579317) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9102218) q[0];
sx q[0];
rz(-1.3274258) q[0];
sx q[0];
rz(-0.87941186) q[0];
rz(-2.544417) q[2];
sx q[2];
rz(-0.9998601) q[2];
sx q[2];
rz(-1.2025637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5722583) q[1];
sx q[1];
rz(-2.5037994) q[1];
sx q[1];
rz(-0.34697726) q[1];
x q[2];
rz(-0.66863184) q[3];
sx q[3];
rz(-0.9210081) q[3];
sx q[3];
rz(1.1710407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-2.5518899) q[2];
rz(1.8288745) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(-2.5780799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967459) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(0.27031159) q[0];
rz(0.42701328) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-2.7332773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9664551) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(-1.4428663) q[0];
rz(-pi) q[1];
rz(-0.87323453) q[2];
sx q[2];
rz(-1.6484652) q[2];
sx q[2];
rz(2.6447062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85874365) q[1];
sx q[1];
rz(-0.74902422) q[1];
sx q[1];
rz(2.2006459) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.116577) q[3];
sx q[3];
rz(-0.72598476) q[3];
sx q[3];
rz(0.67949142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2030877) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(0.92424029) q[2];
rz(0.83485323) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95558178) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(0.36744776) q[0];
rz(-2.5732749) q[1];
sx q[1];
rz(-1.3325997) q[1];
sx q[1];
rz(-0.41499358) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85431722) q[0];
sx q[0];
rz(-1.1181896) q[0];
sx q[0];
rz(-0.29083473) q[0];
rz(-0.52936036) q[2];
sx q[2];
rz(-1.0870013) q[2];
sx q[2];
rz(2.8279378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7447978) q[1];
sx q[1];
rz(-2.2992981) q[1];
sx q[1];
rz(0.70835967) q[1];
rz(-pi) q[2];
rz(-1.9165809) q[3];
sx q[3];
rz(-1.5659075) q[3];
sx q[3];
rz(-1.1556311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(0.83317327) q[2];
rz(-0.35946515) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-0.014884431) q[0];
rz(2.0170276) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(-0.10717779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2994605) q[0];
sx q[0];
rz(-2.3981574) q[0];
sx q[0];
rz(-1.3495665) q[0];
rz(-pi) q[1];
rz(2.6769019) q[2];
sx q[2];
rz(-2.26311) q[2];
sx q[2];
rz(-1.9141045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9669144) q[1];
sx q[1];
rz(-1.0941097) q[1];
sx q[1];
rz(-1.7382453) q[1];
rz(-pi) q[2];
rz(0.26132432) q[3];
sx q[3];
rz(-1.7423986) q[3];
sx q[3];
rz(1.8971407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.3078835) q[2];
sx q[2];
rz(-0.85961071) q[2];
rz(2.882242) q[3];
sx q[3];
rz(-1.7813464) q[3];
sx q[3];
rz(0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49941007) q[0];
sx q[0];
rz(-0.12633093) q[0];
sx q[0];
rz(1.1580178) q[0];
rz(-2.6462789) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(0.29702979) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3619694) q[0];
sx q[0];
rz(-1.2517778) q[0];
sx q[0];
rz(3.1373128) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0115602) q[2];
sx q[2];
rz(-0.30270019) q[2];
sx q[2];
rz(-1.5695733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37849879) q[1];
sx q[1];
rz(-1.2796254) q[1];
sx q[1];
rz(-1.4062455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.591743) q[3];
sx q[3];
rz(-0.6626752) q[3];
sx q[3];
rz(-1.7998723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.872252) q[2];
sx q[2];
rz(-2.7597235) q[2];
sx q[2];
rz(2.2807109) q[2];
rz(0.46368972) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(-2.004682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9853482) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(1.3337749) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(-2.1352519) q[2];
sx q[2];
rz(-2.153572) q[2];
sx q[2];
rz(-1.5700271) q[2];
rz(-2.648223) q[3];
sx q[3];
rz(-1.2262288) q[3];
sx q[3];
rz(-0.83362383) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
