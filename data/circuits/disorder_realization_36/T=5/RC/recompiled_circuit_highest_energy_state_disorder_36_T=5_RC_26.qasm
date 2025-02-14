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
rz(-0.90717301) q[0];
sx q[0];
rz(-2.2450759) q[0];
sx q[0];
rz(-3.1007015) q[0];
rz(1.8812802) q[1];
sx q[1];
rz(-0.44668302) q[1];
sx q[1];
rz(3.0466411) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10863607) q[0];
sx q[0];
rz(-2.6755736) q[0];
sx q[0];
rz(-2.6172383) q[0];
x q[1];
rz(-0.55741252) q[2];
sx q[2];
rz(-1.7669618) q[2];
sx q[2];
rz(-0.74166751) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8548838) q[1];
sx q[1];
rz(-2.2818533) q[1];
sx q[1];
rz(-0.49763007) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0980174) q[3];
sx q[3];
rz(-1.9736787) q[3];
sx q[3];
rz(1.948818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8591259) q[2];
sx q[2];
rz(-0.60255113) q[2];
sx q[2];
rz(1.1653384) q[2];
rz(-0.91853842) q[3];
sx q[3];
rz(-3.0375752) q[3];
sx q[3];
rz(1.2134086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391881) q[0];
sx q[0];
rz(-2.5027051) q[0];
sx q[0];
rz(-0.30493394) q[0];
rz(-2.1091499) q[1];
sx q[1];
rz(-2.86125) q[1];
sx q[1];
rz(-2.2974744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979008) q[0];
sx q[0];
rz(-1.4183174) q[0];
sx q[0];
rz(-1.3293847) q[0];
x q[1];
rz(0.12298583) q[2];
sx q[2];
rz(-1.7779576) q[2];
sx q[2];
rz(2.5442051) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55791486) q[1];
sx q[1];
rz(-2.6632383) q[1];
sx q[1];
rz(-2.5972523) q[1];
rz(-pi) q[2];
x q[2];
rz(0.02133298) q[3];
sx q[3];
rz(-2.5144684) q[3];
sx q[3];
rz(0.54652484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9742763) q[2];
sx q[2];
rz(-0.68334371) q[2];
sx q[2];
rz(0.58504504) q[2];
rz(-0.85917464) q[3];
sx q[3];
rz(-1.9935358) q[3];
sx q[3];
rz(2.008321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0228731) q[0];
sx q[0];
rz(-2.3206503) q[0];
sx q[0];
rz(3.0430479) q[0];
rz(0.56602829) q[1];
sx q[1];
rz(-2.2955344) q[1];
sx q[1];
rz(0.50618323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4249742) q[0];
sx q[0];
rz(-2.0831265) q[0];
sx q[0];
rz(1.1879735) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7012791) q[2];
sx q[2];
rz(-0.8699421) q[2];
sx q[2];
rz(0.93849692) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58193026) q[1];
sx q[1];
rz(-0.25296695) q[1];
sx q[1];
rz(2.0848666) q[1];
rz(0.99876253) q[3];
sx q[3];
rz(-1.2353829) q[3];
sx q[3];
rz(-2.4728925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11518662) q[2];
sx q[2];
rz(-1.9280484) q[2];
sx q[2];
rz(3.0565267) q[2];
rz(-2.3885942) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(-1.5299214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6104777) q[0];
sx q[0];
rz(-1.6401289) q[0];
sx q[0];
rz(2.6016972) q[0];
rz(0.029622948) q[1];
sx q[1];
rz(-0.74631515) q[1];
sx q[1];
rz(1.7866887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82096174) q[0];
sx q[0];
rz(-2.3665049) q[0];
sx q[0];
rz(-1.8863999) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71075534) q[2];
sx q[2];
rz(-1.8216933) q[2];
sx q[2];
rz(2.9352842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78614178) q[1];
sx q[1];
rz(-2.1340573) q[1];
sx q[1];
rz(-1.8763297) q[1];
rz(-0.79014312) q[3];
sx q[3];
rz(-1.6340268) q[3];
sx q[3];
rz(-2.9823398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4386091) q[2];
sx q[2];
rz(-0.97038022) q[2];
sx q[2];
rz(2.2461829) q[2];
rz(-1.6926951) q[3];
sx q[3];
rz(-1.9796895) q[3];
sx q[3];
rz(-2.7321775) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92622906) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(-3.1410425) q[0];
rz(2.6161361) q[1];
sx q[1];
rz(-2.2976687) q[1];
sx q[1];
rz(0.97132436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1570515) q[0];
sx q[0];
rz(-1.2573842) q[0];
sx q[0];
rz(2.393852) q[0];
x q[1];
rz(2.7883456) q[2];
sx q[2];
rz(-1.9611036) q[2];
sx q[2];
rz(-0.82510222) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7441874) q[1];
sx q[1];
rz(-2.1402855) q[1];
sx q[1];
rz(2.2389212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2574785) q[3];
sx q[3];
rz(-1.7755812) q[3];
sx q[3];
rz(1.1623032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690669) q[2];
sx q[2];
rz(-1.0774287) q[2];
sx q[2];
rz(1.7871008) q[2];
rz(0.6012249) q[3];
sx q[3];
rz(-1.0433334) q[3];
sx q[3];
rz(-1.575298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71156597) q[0];
sx q[0];
rz(-2.8373748) q[0];
sx q[0];
rz(2.8416204) q[0];
rz(2.1638339) q[1];
sx q[1];
rz(-0.58039665) q[1];
sx q[1];
rz(-0.93394867) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47563206) q[0];
sx q[0];
rz(-2.7802659) q[0];
sx q[0];
rz(0.40478171) q[0];
rz(0.61638083) q[2];
sx q[2];
rz(-1.167815) q[2];
sx q[2];
rz(2.4015534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63864795) q[1];
sx q[1];
rz(-2.3343122) q[1];
sx q[1];
rz(-0.32420968) q[1];
x q[2];
rz(-1.9225695) q[3];
sx q[3];
rz(-1.1922424) q[3];
sx q[3];
rz(1.5768676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0710435) q[2];
sx q[2];
rz(-1.0087548) q[2];
sx q[2];
rz(2.0096931) q[2];
rz(1.2567358) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(1.7251714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32135949) q[0];
sx q[0];
rz(-1.8978523) q[0];
sx q[0];
rz(-2.5309122) q[0];
rz(-2.3848379) q[1];
sx q[1];
rz(-0.38833955) q[1];
sx q[1];
rz(1.6441708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066787) q[0];
sx q[0];
rz(-1.4467708) q[0];
sx q[0];
rz(-1.5962315) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36368215) q[2];
sx q[2];
rz(-1.1712345) q[2];
sx q[2];
rz(-1.3043208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6886474) q[1];
sx q[1];
rz(-1.4645836) q[1];
sx q[1];
rz(2.4912686) q[1];
rz(0.11388643) q[3];
sx q[3];
rz(-2.4166757) q[3];
sx q[3];
rz(-2.8020086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9553767) q[2];
sx q[2];
rz(-2.4775041) q[2];
sx q[2];
rz(2.5353298) q[2];
rz(2.9210505) q[3];
sx q[3];
rz(-1.3371779) q[3];
sx q[3];
rz(-2.7748599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99899387) q[0];
sx q[0];
rz(-1.6479011) q[0];
sx q[0];
rz(2.2338474) q[0];
rz(-1.6084464) q[1];
sx q[1];
rz(-0.95450675) q[1];
sx q[1];
rz(-0.018891637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2174199) q[0];
sx q[0];
rz(-1.1013563) q[0];
sx q[0];
rz(-1.7007692) q[0];
x q[1];
rz(0.18757815) q[2];
sx q[2];
rz(-2.8370428) q[2];
sx q[2];
rz(1.9001324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50190137) q[1];
sx q[1];
rz(-2.0328201) q[1];
sx q[1];
rz(-1.5359182) q[1];
rz(-pi) q[2];
rz(0.83062474) q[3];
sx q[3];
rz(-2.0181832) q[3];
sx q[3];
rz(0.4448286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7682401) q[2];
sx q[2];
rz(-1.5529996) q[2];
sx q[2];
rz(0.24660435) q[2];
rz(-1.2982347) q[3];
sx q[3];
rz(-1.4539098) q[3];
sx q[3];
rz(-1.8728144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3835555) q[0];
sx q[0];
rz(-2.0691431) q[0];
sx q[0];
rz(-2.2776336) q[0];
rz(1.9268688) q[1];
sx q[1];
rz(-1.1178958) q[1];
sx q[1];
rz(0.58089677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6075077) q[0];
sx q[0];
rz(-2.0568741) q[0];
sx q[0];
rz(2.3838504) q[0];
rz(-0.038453416) q[2];
sx q[2];
rz(-2.3068301) q[2];
sx q[2];
rz(-2.9308386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8859007) q[1];
sx q[1];
rz(-1.2322353) q[1];
sx q[1];
rz(2.5060966) q[1];
rz(2.3623051) q[3];
sx q[3];
rz(-1.3715991) q[3];
sx q[3];
rz(1.3785386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2555799) q[2];
sx q[2];
rz(-2.6801127) q[2];
sx q[2];
rz(0.18730051) q[2];
rz(-0.97697941) q[3];
sx q[3];
rz(-1.5005485) q[3];
sx q[3];
rz(1.2788844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5115857) q[0];
sx q[0];
rz(-0.66660175) q[0];
sx q[0];
rz(-1.2641719) q[0];
rz(0.97201792) q[1];
sx q[1];
rz(-0.7518026) q[1];
sx q[1];
rz(1.972563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326871) q[0];
sx q[0];
rz(-1.4356614) q[0];
sx q[0];
rz(-1.9072471) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9355785) q[2];
sx q[2];
rz(-1.9309461) q[2];
sx q[2];
rz(-2.3003038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5981751) q[1];
sx q[1];
rz(-1.4663883) q[1];
sx q[1];
rz(2.2648952) q[1];
rz(-pi) q[2];
rz(0.31619483) q[3];
sx q[3];
rz(-1.4619816) q[3];
sx q[3];
rz(-1.3294045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4740037) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(-0.37707314) q[2];
rz(-0.22370473) q[3];
sx q[3];
rz(-2.5535899) q[3];
sx q[3];
rz(1.912775) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2558462) q[0];
sx q[0];
rz(-1.3930014) q[0];
sx q[0];
rz(-0.2572671) q[0];
rz(0.59355758) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(2.42057) q[2];
sx q[2];
rz(-0.99698721) q[2];
sx q[2];
rz(0.1884603) q[2];
rz(-2.138924) q[3];
sx q[3];
rz(-0.46350652) q[3];
sx q[3];
rz(1.9476611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
