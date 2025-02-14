OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(1.6530871) q[0];
rz(1.1324963) q[1];
sx q[1];
rz(-2.7096665) q[1];
sx q[1];
rz(-0.63634029) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778438) q[0];
sx q[0];
rz(-1.1882458) q[0];
sx q[0];
rz(-0.93884672) q[0];
rz(-pi) q[1];
rz(1.1875528) q[2];
sx q[2];
rz(-1.2072717) q[2];
sx q[2];
rz(-2.5853047) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8752567) q[1];
sx q[1];
rz(-0.41558973) q[1];
sx q[1];
rz(-1.2591836) q[1];
rz(-pi) q[2];
rz(0.79145517) q[3];
sx q[3];
rz(-0.88578445) q[3];
sx q[3];
rz(-2.5572957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1119614) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-2.4911528) q[2];
rz(1.316831) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(3.0397547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0874852) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(-0.45463872) q[0];
rz(-3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(2.815411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9232193) q[0];
sx q[0];
rz(-2.120864) q[0];
sx q[0];
rz(-2.0327507) q[0];
rz(-0.97450476) q[2];
sx q[2];
rz(-0.49879227) q[2];
sx q[2];
rz(-1.5128262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61039576) q[1];
sx q[1];
rz(-0.67761868) q[1];
sx q[1];
rz(-2.4536831) q[1];
rz(0.91008998) q[3];
sx q[3];
rz(-1.6636563) q[3];
sx q[3];
rz(1.4610491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0252016) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(2.1873059) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-2.818675) q[3];
sx q[3];
rz(3.1402816) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8702451) q[0];
sx q[0];
rz(-1.1529237) q[0];
sx q[0];
rz(1.9648319) q[0];
rz(-2.7765043) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(-1.6129859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919892) q[0];
sx q[0];
rz(-0.024734298) q[0];
sx q[0];
rz(-1.0036841) q[0];
rz(-pi) q[1];
rz(-1.4055785) q[2];
sx q[2];
rz(-2.6927136) q[2];
sx q[2];
rz(0.53476483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4903725) q[1];
sx q[1];
rz(-1.7224604) q[1];
sx q[1];
rz(-0.030812736) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4432676) q[3];
sx q[3];
rz(-0.089524448) q[3];
sx q[3];
rz(-1.2379299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8927346) q[2];
sx q[2];
rz(-0.55344075) q[2];
sx q[2];
rz(1.152732) q[2];
rz(1.8966127) q[3];
sx q[3];
rz(-2.1608519) q[3];
sx q[3];
rz(-0.28321701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1588441) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(0.68956476) q[0];
rz(3.116563) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(-1.4208581) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11224225) q[0];
sx q[0];
rz(-0.94932244) q[0];
sx q[0];
rz(-1.6981324) q[0];
rz(-pi) q[1];
rz(1.1941853) q[2];
sx q[2];
rz(-2.6832697) q[2];
sx q[2];
rz(-1.5397746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.082877872) q[1];
sx q[1];
rz(-2.0738462) q[1];
sx q[1];
rz(1.5392667) q[1];
rz(-pi) q[2];
rz(-2.8042364) q[3];
sx q[3];
rz(-1.8280262) q[3];
sx q[3];
rz(1.4816062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.932852) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(-1.9746732) q[2];
rz(-0.44421998) q[3];
sx q[3];
rz(-1.9053713) q[3];
sx q[3];
rz(-0.3096295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8532448) q[0];
sx q[0];
rz(-1.7720368) q[0];
sx q[0];
rz(0.41611588) q[0];
rz(-0.5101282) q[1];
sx q[1];
rz(-2.3704539) q[1];
sx q[1];
rz(0.77521926) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0508779) q[0];
sx q[0];
rz(-2.4950301) q[0];
sx q[0];
rz(-0.334686) q[0];
x q[1];
rz(-0.40521605) q[2];
sx q[2];
rz(-1.4034668) q[2];
sx q[2];
rz(2.2115603) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3415592) q[1];
sx q[1];
rz(-1.4455035) q[1];
sx q[1];
rz(0.11334669) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6121077) q[3];
sx q[3];
rz(-2.4765402) q[3];
sx q[3];
rz(2.3163026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9267209) q[2];
sx q[2];
rz(-0.99433172) q[2];
sx q[2];
rz(1.035824) q[2];
rz(2.9521247) q[3];
sx q[3];
rz(-0.47179705) q[3];
sx q[3];
rz(2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254461) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(-2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(2.0004418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095207) q[0];
sx q[0];
rz(-0.70200015) q[0];
sx q[0];
rz(-1.899748) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1823523) q[2];
sx q[2];
rz(-1.5634166) q[2];
sx q[2];
rz(0.5317229) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49887744) q[1];
sx q[1];
rz(-1.6756454) q[1];
sx q[1];
rz(-0.14123209) q[1];
x q[2];
rz(1.136888) q[3];
sx q[3];
rz(-1.5818137) q[3];
sx q[3];
rz(2.9316255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5100539) q[2];
sx q[2];
rz(-2.4271991) q[2];
sx q[2];
rz(-1.7047403) q[2];
rz(2.0884183) q[3];
sx q[3];
rz(-1.6097693) q[3];
sx q[3];
rz(-2.6385782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18008867) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(0.12251138) q[0];
rz(-0.65613121) q[1];
sx q[1];
rz(-1.2686814) q[1];
sx q[1];
rz(1.3414541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91514523) q[0];
sx q[0];
rz(-1.6065803) q[0];
sx q[0];
rz(0.045374845) q[0];
x q[1];
rz(1.4627005) q[2];
sx q[2];
rz(-1.8524395) q[2];
sx q[2];
rz(3.0710655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0602051) q[1];
sx q[1];
rz(-1.3045746) q[1];
sx q[1];
rz(2.2039444) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9060494) q[3];
sx q[3];
rz(-0.45484124) q[3];
sx q[3];
rz(-1.7419635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5682257) q[2];
sx q[2];
rz(-1.919701) q[2];
sx q[2];
rz(1.4637671) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(-1.8370139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463335) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(2.1268225) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(-2.0789304) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39358562) q[0];
sx q[0];
rz(-2.9655122) q[0];
sx q[0];
rz(-0.3349456) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6133461) q[2];
sx q[2];
rz(-0.88111231) q[2];
sx q[2];
rz(1.4572342) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6148551) q[1];
sx q[1];
rz(-1.4397845) q[1];
sx q[1];
rz(-3.0796771) q[1];
x q[2];
rz(2.1001535) q[3];
sx q[3];
rz(-0.58515195) q[3];
sx q[3];
rz(-3.0781259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.3883611) q[2];
sx q[2];
rz(1.3512705) q[2];
rz(-1.9225559) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(-2.3082699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(-2.2341527) q[0];
rz(-2.9145248) q[1];
sx q[1];
rz(-2.0957969) q[1];
sx q[1];
rz(-2.2423832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9010021) q[0];
sx q[0];
rz(-1.0623029) q[0];
sx q[0];
rz(-2.2700538) q[0];
rz(-pi) q[1];
rz(1.4883409) q[2];
sx q[2];
rz(-1.4926148) q[2];
sx q[2];
rz(-2.9417335) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6688706) q[1];
sx q[1];
rz(-0.54520118) q[1];
sx q[1];
rz(1.9342058) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3189018) q[3];
sx q[3];
rz(-1.5334276) q[3];
sx q[3];
rz(-2.5189375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1850618) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(2.5409307) q[2];
rz(1.0968084) q[3];
sx q[3];
rz(-2.5513702) q[3];
sx q[3];
rz(0.87305951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.3908865) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(-0.98439687) q[0];
rz(-0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(-1.9459928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0988783) q[0];
sx q[0];
rz(-2.1568568) q[0];
sx q[0];
rz(1.1436966) q[0];
x q[1];
rz(-0.97862794) q[2];
sx q[2];
rz(-0.76972967) q[2];
sx q[2];
rz(-0.20537381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3477563) q[1];
sx q[1];
rz(-1.4897921) q[1];
sx q[1];
rz(-1.9543314) q[1];
x q[2];
rz(-0.77287425) q[3];
sx q[3];
rz(-1.2307271) q[3];
sx q[3];
rz(-0.085207663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.055858) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(2.1585507) q[2];
rz(-0.25162697) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(2.1380641) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49878237) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(0.21435195) q[2];
sx q[2];
rz(-0.8690693) q[2];
sx q[2];
rz(-2.2590841) q[2];
rz(-0.33369251) q[3];
sx q[3];
rz(-1.8751388) q[3];
sx q[3];
rz(1.4061385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
