OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9031154) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(1.205501) q[0];
x q[1];
rz(-1.5654972) q[2];
sx q[2];
rz(-1.0056579) q[2];
sx q[2];
rz(2.0084755) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.002418) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-0.26267085) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7847071) q[3];
sx q[3];
rz(-0.89324739) q[3];
sx q[3];
rz(2.4524636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-3.1343592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7318856) q[0];
sx q[0];
rz(-1.9623555) q[0];
sx q[0];
rz(-3.0407228) q[0];
rz(-pi) q[1];
rz(-2.8480808) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(-2.7345865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0011598) q[1];
sx q[1];
rz(-1.7908485) q[1];
sx q[1];
rz(-2.7697255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(-0.12156045) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.7720222) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(2.8799768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3900671) q[0];
sx q[0];
rz(-1.5528423) q[0];
sx q[0];
rz(3.1319588) q[0];
rz(-pi) q[1];
rz(-3.0171266) q[2];
sx q[2];
rz(-0.83041588) q[2];
sx q[2];
rz(2.894573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1640472) q[1];
sx q[1];
rz(-2.3512212) q[1];
sx q[1];
rz(2.9708944) q[1];
rz(-pi) q[2];
rz(1.5393799) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(-1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.2483695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.473339) q[0];
sx q[0];
rz(-0.11675294) q[0];
sx q[0];
rz(-0.97572414) q[0];
x q[1];
rz(0.20935697) q[2];
sx q[2];
rz(-1.7484089) q[2];
sx q[2];
rz(-2.3043485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.215938) q[1];
sx q[1];
rz(-2.2207894) q[1];
sx q[1];
rz(1.710612) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8321886) q[3];
sx q[3];
rz(-0.86463118) q[3];
sx q[3];
rz(-0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(-1.2083496) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(-0.77520448) q[0];
rz(-2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(2.2391589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8309801) q[0];
sx q[0];
rz(-1.5251625) q[0];
sx q[0];
rz(-2.5020585) q[0];
rz(-0.4544223) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(-2.1511252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.21197) q[1];
sx q[1];
rz(-2.3096482) q[1];
sx q[1];
rz(-0.79939876) q[1];
rz(-0.0028354672) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(1.0032652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(-1.6852089) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4393282) q[0];
sx q[0];
rz(-1.4685255) q[0];
sx q[0];
rz(-3.0754473) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1226419) q[2];
sx q[2];
rz(-1.1345703) q[2];
sx q[2];
rz(0.052554616) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6946053) q[1];
sx q[1];
rz(-1.9896549) q[1];
sx q[1];
rz(2.6161731) q[1];
rz(0.73268907) q[3];
sx q[3];
rz(-1.9273888) q[3];
sx q[3];
rz(-2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(-1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-2.1988791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276886) q[0];
sx q[0];
rz(-1.9644992) q[0];
sx q[0];
rz(-0.64481553) q[0];
rz(-pi) q[1];
rz(-1.7985293) q[2];
sx q[2];
rz(-1.7274389) q[2];
sx q[2];
rz(1.1802955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.645694) q[1];
sx q[1];
rz(-1.2518479) q[1];
sx q[1];
rz(-1.1369399) q[1];
x q[2];
rz(2.0122635) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(2.7260775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-0.38267246) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(0.27012816) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89284183) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(2.6177004) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21247272) q[2];
sx q[2];
rz(-2.5310235) q[2];
sx q[2];
rz(0.072710466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9406416) q[1];
sx q[1];
rz(-2.0307699) q[1];
sx q[1];
rz(0.56401395) q[1];
x q[2];
rz(0.038589434) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(-0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(1.930687) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(1.680826) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54366771) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(-3.0790867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3912348) q[2];
sx q[2];
rz(-1.1968687) q[2];
sx q[2];
rz(2.0342846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2673805) q[1];
sx q[1];
rz(-1.0460209) q[1];
sx q[1];
rz(-1.5449779) q[1];
rz(-2.3567696) q[3];
sx q[3];
rz(-1.2900969) q[3];
sx q[3];
rz(0.37898889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(0.6955859) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(2.4263884) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(0.38063231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779553) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(1.3689343) q[0];
rz(-pi) q[1];
rz(-2.5194174) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(3.0058793) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14643529) q[1];
sx q[1];
rz(-2.6019115) q[1];
sx q[1];
rz(1.1848918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57772824) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(-1.1364737) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-1.0271172) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-0.2480416) q[2];
sx q[2];
rz(-1.7938062) q[2];
sx q[2];
rz(1.2563406) q[2];
rz(1.1313664) q[3];
sx q[3];
rz(-1.7527179) q[3];
sx q[3];
rz(-2.6408165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];