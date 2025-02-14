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
rz(-2.0877617) q[0];
sx q[0];
rz(-0.55019903) q[0];
sx q[0];
rz(-1.3681816) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36461386) q[0];
sx q[0];
rz(-1.783051) q[0];
sx q[0];
rz(1.7745166) q[0];
x q[1];
rz(-2.3362581) q[2];
sx q[2];
rz(-0.67395681) q[2];
sx q[2];
rz(2.0827302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3319407) q[1];
sx q[1];
rz(-0.9729079) q[1];
sx q[1];
rz(-0.34418587) q[1];
x q[2];
rz(-3.1038001) q[3];
sx q[3];
rz(-2.0366014) q[3];
sx q[3];
rz(-0.58096262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54174417) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(2.8372676) q[3];
sx q[3];
rz(-2.2797238) q[3];
sx q[3];
rz(-0.26722515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32352725) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(-0.7793119) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(1.13387) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4138105) q[0];
sx q[0];
rz(-1.4725794) q[0];
sx q[0];
rz(-1.1350687) q[0];
rz(-pi) q[1];
rz(-1.5120087) q[2];
sx q[2];
rz(-1.4530621) q[2];
sx q[2];
rz(1.9699772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8633921) q[1];
sx q[1];
rz(-2.8542633) q[1];
sx q[1];
rz(-1.1817688) q[1];
rz(-pi) q[2];
rz(-1.5690345) q[3];
sx q[3];
rz(-2.632318) q[3];
sx q[3];
rz(0.60520303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4506932) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(-2.9449985) q[2];
rz(0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(-1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75152385) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(0.73177904) q[0];
rz(2.7732908) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(2.7751353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55933773) q[0];
sx q[0];
rz(-2.8213503) q[0];
sx q[0];
rz(-0.59424627) q[0];
rz(-pi) q[1];
rz(0.24213893) q[2];
sx q[2];
rz(-1.79091) q[2];
sx q[2];
rz(1.5907703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7190105) q[1];
sx q[1];
rz(-1.9937464) q[1];
sx q[1];
rz(-0.99653901) q[1];
rz(-pi) q[2];
rz(-2.5440574) q[3];
sx q[3];
rz(-2.1794277) q[3];
sx q[3];
rz(-0.017692117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9590108) q[2];
sx q[2];
rz(-0.77114791) q[2];
sx q[2];
rz(-2.2491573) q[2];
rz(-0.45448947) q[3];
sx q[3];
rz(-0.82388866) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3733805) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(0.40661231) q[0];
rz(-0.82798249) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(2.3060395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53518772) q[0];
sx q[0];
rz(-1.390103) q[0];
sx q[0];
rz(3.101493) q[0];
x q[1];
rz(1.2242975) q[2];
sx q[2];
rz(-1.3801127) q[2];
sx q[2];
rz(2.4257223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20488508) q[1];
sx q[1];
rz(-2.6815273) q[1];
sx q[1];
rz(0.33554828) q[1];
rz(-pi) q[2];
rz(1.9560247) q[3];
sx q[3];
rz(-1.6611757) q[3];
sx q[3];
rz(-0.65672311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39997175) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(2.2201904) q[2];
rz(-1.7737927) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(3.0000946) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(-1.5171492) q[1];
sx q[1];
rz(-2.8327063) q[1];
sx q[1];
rz(1.3295757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3388728) q[0];
sx q[0];
rz(-0.021139806) q[0];
sx q[0];
rz(-2.6129524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9022568) q[2];
sx q[2];
rz(-2.1194601) q[2];
sx q[2];
rz(-1.885883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5160421) q[1];
sx q[1];
rz(-0.4333638) q[1];
sx q[1];
rz(-2.0651555) q[1];
x q[2];
rz(2.1383189) q[3];
sx q[3];
rz(-1.7470932) q[3];
sx q[3];
rz(0.76667537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(-0.55207437) q[2];
rz(2.3479346) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8296705) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(2.435834) q[0];
rz(-0.13748473) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(0.71298832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891083) q[0];
sx q[0];
rz(-2.6875067) q[0];
sx q[0];
rz(-2.5986791) q[0];
rz(-pi) q[1];
rz(2.001686) q[2];
sx q[2];
rz(-2.0451114) q[2];
sx q[2];
rz(-2.9805019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56604993) q[1];
sx q[1];
rz(-0.41047305) q[1];
sx q[1];
rz(1.3259352) q[1];
rz(-pi) q[2];
rz(1.2778132) q[3];
sx q[3];
rz(-1.4594541) q[3];
sx q[3];
rz(-0.35929832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9299499) q[2];
sx q[2];
rz(-2.2723618) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(-0.12187135) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(-1.6596644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.876494) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(0.70736831) q[0];
rz(2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.7792938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14868098) q[0];
sx q[0];
rz(-1.8472965) q[0];
sx q[0];
rz(-0.25603237) q[0];
rz(-pi) q[1];
rz(2.7168305) q[2];
sx q[2];
rz(-2.2519037) q[2];
sx q[2];
rz(2.8515138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8352812) q[1];
sx q[1];
rz(-1.0544485) q[1];
sx q[1];
rz(0.55037127) q[1];
x q[2];
rz(1.6860289) q[3];
sx q[3];
rz(-2.0243905) q[3];
sx q[3];
rz(0.99059243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1661487) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55050945) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(1.9367223) q[0];
rz(-0.51271802) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(0.18096322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2406684) q[0];
sx q[0];
rz(-2.5818733) q[0];
sx q[0];
rz(-0.75738917) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7082735) q[2];
sx q[2];
rz(-1.9172693) q[2];
sx q[2];
rz(0.32127831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4875723) q[1];
sx q[1];
rz(-2.1250238) q[1];
sx q[1];
rz(-0.2481064) q[1];
x q[2];
rz(-3.0309903) q[3];
sx q[3];
rz(-1.4975274) q[3];
sx q[3];
rz(1.4721118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60244954) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(0.096972018) q[2];
rz(0.55475956) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612514) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(3.0338147) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(-1.5239747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1861098) q[0];
sx q[0];
rz(-0.26676565) q[0];
sx q[0];
rz(-1.0590932) q[0];
rz(-pi) q[1];
rz(-1.0650159) q[2];
sx q[2];
rz(-0.98978251) q[2];
sx q[2];
rz(0.53321028) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7858582) q[1];
sx q[1];
rz(-1.5275035) q[1];
sx q[1];
rz(2.0329352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49160853) q[3];
sx q[3];
rz(-1.5244686) q[3];
sx q[3];
rz(2.2674048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9109351) q[2];
sx q[2];
rz(-2.6378938) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(0.23760992) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(-2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199353) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(-2.7811116) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(-1.9182659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66401615) q[0];
sx q[0];
rz(-1.7682791) q[0];
sx q[0];
rz(-2.5504179) q[0];
rz(2.3476841) q[2];
sx q[2];
rz(-2.0092699) q[2];
sx q[2];
rz(-2.7011342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66050038) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(-0.41892799) q[1];
x q[2];
rz(1.0343219) q[3];
sx q[3];
rz(-1.4015028) q[3];
sx q[3];
rz(-3.0790429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90784043) q[2];
sx q[2];
rz(-1.7641822) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(-0.20445538) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(-2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37888708) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(-1.7145722) q[1];
sx q[1];
rz(-0.49976977) q[1];
sx q[1];
rz(1.5907092) q[1];
rz(-0.84080055) q[2];
sx q[2];
rz(-1.9304095) q[2];
sx q[2];
rz(-0.87113397) q[2];
rz(-1.5132202) q[3];
sx q[3];
rz(-1.7840458) q[3];
sx q[3];
rz(1.9174867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
