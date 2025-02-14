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
rz(-2.3772216) q[0];
sx q[0];
rz(-1.8005014) q[0];
sx q[0];
rz(-2.2361225) q[0];
rz(-1.634693) q[1];
sx q[1];
rz(-0.96087471) q[1];
sx q[1];
rz(-1.5009872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849821) q[0];
sx q[0];
rz(-2.2652103) q[0];
sx q[0];
rz(-0.25882369) q[0];
rz(-pi) q[1];
rz(-1.2086966) q[2];
sx q[2];
rz(-0.46572567) q[2];
sx q[2];
rz(-3.1052239) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0762877) q[1];
sx q[1];
rz(-1.5283547) q[1];
sx q[1];
rz(-1.5941335) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0910283) q[3];
sx q[3];
rz(-2.5135698) q[3];
sx q[3];
rz(-1.764738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(-2.8196715) q[2];
rz(2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(2.6296997) q[0];
rz(1.5882675) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(1.1221251) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.666709) q[0];
sx q[0];
rz(-2.8132952) q[0];
sx q[0];
rz(-0.31995456) q[0];
rz(-3.1013158) q[2];
sx q[2];
rz(-2.0000519) q[2];
sx q[2];
rz(-1.5778936) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.304804) q[1];
sx q[1];
rz(-1.2366017) q[1];
sx q[1];
rz(1.1172953) q[1];
rz(-pi) q[2];
rz(-2.0418704) q[3];
sx q[3];
rz(-2.2212914) q[3];
sx q[3];
rz(-0.6944523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3731709) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(2.6785417) q[2];
rz(2.566346) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082044) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(-2.7114765) q[0];
rz(2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(0.63708416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642536) q[0];
sx q[0];
rz(-1.6116983) q[0];
sx q[0];
rz(-3.1275355) q[0];
x q[1];
rz(3.1153684) q[2];
sx q[2];
rz(-2.3560239) q[2];
sx q[2];
rz(-2.3303967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3780669) q[1];
sx q[1];
rz(-2.3664306) q[1];
sx q[1];
rz(-0.35825348) q[1];
rz(0.51882842) q[3];
sx q[3];
rz(-1.9137376) q[3];
sx q[3];
rz(-2.1800621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-0.72506881) q[2];
rz(1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(-2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722039) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(1.4440906) q[0];
rz(-0.048642453) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(-0.28894249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83544448) q[0];
sx q[0];
rz(-1.939538) q[0];
sx q[0];
rz(0.25935092) q[0];
x q[1];
rz(-3.1132675) q[2];
sx q[2];
rz(-1.420317) q[2];
sx q[2];
rz(-0.17720824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1370586) q[1];
sx q[1];
rz(-1.4618131) q[1];
sx q[1];
rz(-1.3957028) q[1];
rz(3.0722202) q[3];
sx q[3];
rz(-0.66646229) q[3];
sx q[3];
rz(1.762076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.284953) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(-0.3802158) q[2];
rz(-2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83288348) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(-1.0943476) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9834866) q[0];
sx q[0];
rz(-1.1855584) q[0];
sx q[0];
rz(-1.4225106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78754707) q[2];
sx q[2];
rz(-2.6733477) q[2];
sx q[2];
rz(0.59216532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.588394) q[1];
sx q[1];
rz(-1.0244245) q[1];
sx q[1];
rz(2.9562034) q[1];
x q[2];
rz(-0.15954475) q[3];
sx q[3];
rz(-2.2189848) q[3];
sx q[3];
rz(1.2153347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(2.1649427) q[2];
rz(-0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831182) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.2716768) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(-1.6082825) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5137702) q[0];
sx q[0];
rz(-1.5830399) q[0];
sx q[0];
rz(-1.6374541) q[0];
x q[1];
rz(1.6663867) q[2];
sx q[2];
rz(-0.098473452) q[2];
sx q[2];
rz(2.5661039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6110797) q[1];
sx q[1];
rz(-0.30834282) q[1];
sx q[1];
rz(-0.18402305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.815237) q[3];
sx q[3];
rz(-0.77271739) q[3];
sx q[3];
rz(-2.6561007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(2.7563654) q[2];
rz(0.084608229) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-2.718524) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794681) q[0];
sx q[0];
rz(-1.6982268) q[0];
sx q[0];
rz(-2.3412933) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0881617) q[2];
sx q[2];
rz(-2.3950999) q[2];
sx q[2];
rz(1.5160402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7518378) q[1];
sx q[1];
rz(-1.0026649) q[1];
sx q[1];
rz(1.1144494) q[1];
x q[2];
rz(2.2479731) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(-2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6323382) q[2];
sx q[2];
rz(-1.961144) q[2];
sx q[2];
rz(0.21305591) q[2];
rz(-2.3325855) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(-1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(-0.2991547) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(2.1655653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19071391) q[0];
sx q[0];
rz(-1.6813206) q[0];
sx q[0];
rz(-2.8948363) q[0];
rz(-1.3930182) q[2];
sx q[2];
rz(-1.7209098) q[2];
sx q[2];
rz(0.77265384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9240771) q[1];
sx q[1];
rz(-1.922632) q[1];
sx q[1];
rz(-2.887243) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3221413) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(-1.6909042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30190793) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(3.0111266) q[2];
rz(-1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94006222) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(1.6424204) q[0];
rz(-2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(-1.3444208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3657065) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(1.1613174) q[0];
rz(-pi) q[1];
rz(-0.70602472) q[2];
sx q[2];
rz(-2.386552) q[2];
sx q[2];
rz(1.4817099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.498913) q[1];
sx q[1];
rz(-2.1331926) q[1];
sx q[1];
rz(-2.845473) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5193941) q[3];
sx q[3];
rz(-2.1902764) q[3];
sx q[3];
rz(-3.0752237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-0.33668533) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(-0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(-1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880599) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.6499299) q[0];
rz(-0.55117575) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(0.59250441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1051286) q[0];
sx q[0];
rz(-1.8499287) q[0];
sx q[0];
rz(1.817784) q[0];
x q[1];
rz(-2.5985322) q[2];
sx q[2];
rz(-2.0351699) q[2];
sx q[2];
rz(-2.4891702) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7043276) q[1];
sx q[1];
rz(-2.6710837) q[1];
sx q[1];
rz(0.60948845) q[1];
rz(-pi) q[2];
rz(-1.5207401) q[3];
sx q[3];
rz(-1.0721237) q[3];
sx q[3];
rz(3.0704481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42980117) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(3.0832624) q[2];
rz(2.77099) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534828) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.7755605) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(1.8283394) q[2];
sx q[2];
rz(-1.0124442) q[2];
sx q[2];
rz(-0.53895216) q[2];
rz(2.1128863) q[3];
sx q[3];
rz(-2.0301314) q[3];
sx q[3];
rz(0.83740656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
