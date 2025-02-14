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
rz(0.82062757) q[0];
sx q[0];
rz(-0.66250116) q[0];
sx q[0];
rz(2.878046) q[0];
rz(1.1072371) q[1];
sx q[1];
rz(4.4237408) q[1];
sx q[1];
rz(10.10034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1206664) q[0];
sx q[0];
rz(-1.5772737) q[0];
sx q[0];
rz(-1.2351888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98904977) q[2];
sx q[2];
rz(-2.0768696) q[2];
sx q[2];
rz(2.842234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61900422) q[1];
sx q[1];
rz(-1.9591041) q[1];
sx q[1];
rz(2.9791645) q[1];
rz(-pi) q[2];
rz(1.0682929) q[3];
sx q[3];
rz(-1.2811105) q[3];
sx q[3];
rz(-2.009249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88966173) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(2.7737854) q[2];
rz(-0.31428549) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(-0.17816003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80938584) q[0];
sx q[0];
rz(-3.0520913) q[0];
sx q[0];
rz(-1.7820763) q[0];
rz(1.2844194) q[1];
sx q[1];
rz(-1.1651243) q[1];
sx q[1];
rz(-0.051609106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29362677) q[0];
sx q[0];
rz(-2.81797) q[0];
sx q[0];
rz(-2.1725562) q[0];
rz(2.935595) q[2];
sx q[2];
rz(-1.2075181) q[2];
sx q[2];
rz(1.8415123) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9549841) q[1];
sx q[1];
rz(-0.37395769) q[1];
sx q[1];
rz(1.5692406) q[1];
x q[2];
rz(-1.9767999) q[3];
sx q[3];
rz(-2.0116802) q[3];
sx q[3];
rz(-1.0025546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.08903683) q[2];
sx q[2];
rz(-0.26036981) q[2];
sx q[2];
rz(-0.54711771) q[2];
rz(1.268092) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(0.33856302) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38076213) q[0];
sx q[0];
rz(-1.0365423) q[0];
sx q[0];
rz(1.3408252) q[0];
rz(0.85397589) q[1];
sx q[1];
rz(-1.2869336) q[1];
sx q[1];
rz(2.8391848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.069347) q[0];
sx q[0];
rz(-2.4503142) q[0];
sx q[0];
rz(1.1290347) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45091429) q[2];
sx q[2];
rz(-1.779777) q[2];
sx q[2];
rz(1.9318888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9309792) q[1];
sx q[1];
rz(-2.180778) q[1];
sx q[1];
rz(-2.239834) q[1];
x q[2];
rz(-2.8676978) q[3];
sx q[3];
rz(-1.0948101) q[3];
sx q[3];
rz(-2.7416503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93620187) q[2];
sx q[2];
rz(-1.8855636) q[2];
sx q[2];
rz(-0.54226792) q[2];
rz(0.99644709) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(-3.0136133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19313136) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(-1.8782072) q[0];
rz(-0.43426934) q[1];
sx q[1];
rz(-0.77249211) q[1];
sx q[1];
rz(-1.451452) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6497191) q[0];
sx q[0];
rz(-1.8611835) q[0];
sx q[0];
rz(0.47652737) q[0];
x q[1];
rz(-0.66777467) q[2];
sx q[2];
rz(-2.479907) q[2];
sx q[2];
rz(-1.1061321) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8122708) q[1];
sx q[1];
rz(-1.656161) q[1];
sx q[1];
rz(-1.4567767) q[1];
x q[2];
rz(-2.1918815) q[3];
sx q[3];
rz(-1.1009645) q[3];
sx q[3];
rz(2.7787152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0749391) q[2];
sx q[2];
rz(-2.1982919) q[2];
sx q[2];
rz(2.2554876) q[2];
rz(-0.0063627176) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(-1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56115443) q[0];
sx q[0];
rz(-0.4564603) q[0];
sx q[0];
rz(-1.891267) q[0];
rz(2.9087032) q[1];
sx q[1];
rz(-0.75585514) q[1];
sx q[1];
rz(-3.0511391) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88921493) q[0];
sx q[0];
rz(-1.4486635) q[0];
sx q[0];
rz(-1.8218763) q[0];
rz(2.6457067) q[2];
sx q[2];
rz(-2.6156938) q[2];
sx q[2];
rz(2.8484707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4793746) q[1];
sx q[1];
rz(-1.5081128) q[1];
sx q[1];
rz(-2.4880243) q[1];
x q[2];
rz(3.0556942) q[3];
sx q[3];
rz(-1.2441564) q[3];
sx q[3];
rz(-0.24490164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6470486) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(1.9336644) q[2];
rz(-2.2203994) q[3];
sx q[3];
rz(-0.84417206) q[3];
sx q[3];
rz(2.6827961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10814609) q[0];
sx q[0];
rz(-1.4220081) q[0];
sx q[0];
rz(-0.49149996) q[0];
rz(-0.72832251) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(2.947015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174343) q[0];
sx q[0];
rz(-1.3493269) q[0];
sx q[0];
rz(0.99507777) q[0];
rz(0.79030444) q[2];
sx q[2];
rz(-2.2111597) q[2];
sx q[2];
rz(-1.1932288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7214365) q[1];
sx q[1];
rz(-2.2682574) q[1];
sx q[1];
rz(-2.2583794) q[1];
x q[2];
rz(0.65690144) q[3];
sx q[3];
rz(-1.7563151) q[3];
sx q[3];
rz(2.3168062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3538094) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(0.0018399012) q[2];
rz(-1.6483824) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(0.28436896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(2.0765685) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-0.5717259) q[1];
sx q[1];
rz(-2.6782742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119448) q[0];
sx q[0];
rz(-1.72763) q[0];
sx q[0];
rz(0.98008396) q[0];
rz(-pi) q[1];
rz(0.58775522) q[2];
sx q[2];
rz(-0.81771933) q[2];
sx q[2];
rz(0.14434926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66760705) q[1];
sx q[1];
rz(-1.7069478) q[1];
sx q[1];
rz(0.48696792) q[1];
rz(-pi) q[2];
rz(-1.4557119) q[3];
sx q[3];
rz(-1.5949378) q[3];
sx q[3];
rz(0.80366443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.00039438417) q[2];
sx q[2];
rz(-1.980915) q[2];
sx q[2];
rz(-0.82994962) q[2];
rz(-1.0478033) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-2.9629663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(-2.3867699) q[0];
rz(2.7524475) q[1];
sx q[1];
rz(-1.4578338) q[1];
sx q[1];
rz(2.7438502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11618092) q[0];
sx q[0];
rz(-1.3137903) q[0];
sx q[0];
rz(-1.2776796) q[0];
x q[1];
rz(2.7806578) q[2];
sx q[2];
rz(-0.58141469) q[2];
sx q[2];
rz(-1.6045398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.039835416) q[1];
sx q[1];
rz(-1.3164742) q[1];
sx q[1];
rz(1.8524342) q[1];
x q[2];
rz(-0.92614545) q[3];
sx q[3];
rz(-0.93346262) q[3];
sx q[3];
rz(-1.6660575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8336746) q[2];
sx q[2];
rz(-0.29611823) q[2];
sx q[2];
rz(0.31416848) q[2];
rz(1.8999772) q[3];
sx q[3];
rz(-1.5887458) q[3];
sx q[3];
rz(-1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.817713) q[0];
sx q[0];
rz(-2.4535024) q[0];
sx q[0];
rz(2.896198) q[0];
rz(-1.2303526) q[1];
sx q[1];
rz(-0.10086682) q[1];
sx q[1];
rz(2.6255677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2986261) q[0];
sx q[0];
rz(-0.43382513) q[0];
sx q[0];
rz(0.4900269) q[0];
x q[1];
rz(0.74036171) q[2];
sx q[2];
rz(-2.2380793) q[2];
sx q[2];
rz(1.1297392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5616331) q[1];
sx q[1];
rz(-1.206534) q[1];
sx q[1];
rz(1.2666525) q[1];
rz(-1.6327758) q[3];
sx q[3];
rz(-1.4345285) q[3];
sx q[3];
rz(-2.800022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55769354) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(2.058513) q[2];
rz(-2.9205868) q[3];
sx q[3];
rz(-0.14328863) q[3];
sx q[3];
rz(-0.7964645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021127064) q[0];
sx q[0];
rz(-0.080862232) q[0];
sx q[0];
rz(0.13588151) q[0];
rz(-1.4295626) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(-0.28724614) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82086876) q[0];
sx q[0];
rz(-2.831376) q[0];
sx q[0];
rz(-2.2287772) q[0];
x q[1];
rz(1.3332149) q[2];
sx q[2];
rz(-1.6529978) q[2];
sx q[2];
rz(2.0978218) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96480709) q[1];
sx q[1];
rz(-0.63399708) q[1];
sx q[1];
rz(0.37668677) q[1];
rz(0.29218896) q[3];
sx q[3];
rz(-2.1351519) q[3];
sx q[3];
rz(0.60114229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91934943) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(0.83402056) q[2];
rz(-2.833994) q[3];
sx q[3];
rz(-2.7643876) q[3];
sx q[3];
rz(-2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.96497) q[0];
sx q[0];
rz(-1.8360092) q[0];
sx q[0];
rz(-1.0259884) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(0.42456498) q[2];
sx q[2];
rz(-1.6909611) q[2];
sx q[2];
rz(1.2517902) q[2];
rz(2.9959804) q[3];
sx q[3];
rz(-1.3619193) q[3];
sx q[3];
rz(-3.0633434) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
