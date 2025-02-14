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
rz(1.6765321) q[0];
sx q[0];
rz(-2.9002011) q[0];
sx q[0];
rz(-0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(3.1225966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26532224) q[0];
sx q[0];
rz(-0.92060584) q[0];
sx q[0];
rz(2.8702535) q[0];
rz(-0.094474205) q[2];
sx q[2];
rz(-1.7565389) q[2];
sx q[2];
rz(-2.8217725) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7380812) q[1];
sx q[1];
rz(-0.36952239) q[1];
sx q[1];
rz(2.6415537) q[1];
rz(-pi) q[2];
rz(-0.98145841) q[3];
sx q[3];
rz(-2.8638726) q[3];
sx q[3];
rz(0.026175682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12048177) q[2];
sx q[2];
rz(-0.80664539) q[2];
sx q[2];
rz(-0.79180229) q[2];
rz(2.0558489) q[3];
sx q[3];
rz(-1.2191685) q[3];
sx q[3];
rz(-2.048548) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2708112) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(-2.812401) q[0];
rz(1.6393433) q[1];
sx q[1];
rz(-2.2396125) q[1];
sx q[1];
rz(2.7194068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6227111) q[0];
sx q[0];
rz(-0.71463138) q[0];
sx q[0];
rz(-2.463752) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29624002) q[2];
sx q[2];
rz(-2.5736817) q[2];
sx q[2];
rz(-2.8366249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93041673) q[1];
sx q[1];
rz(-0.32725829) q[1];
sx q[1];
rz(-0.14463592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6738191) q[3];
sx q[3];
rz(-2.3405082) q[3];
sx q[3];
rz(0.57263206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72055703) q[2];
sx q[2];
rz(-1.9562419) q[2];
sx q[2];
rz(-2.0015008) q[2];
rz(3.1371878) q[3];
sx q[3];
rz(-1.5604228) q[3];
sx q[3];
rz(0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(2.7797101) q[0];
sx q[0];
rz(-3.0411868) q[0];
sx q[0];
rz(-2.7959339) q[0];
rz(-2.0480305) q[1];
sx q[1];
rz(-2.3193017) q[1];
sx q[1];
rz(1.0543157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39798007) q[0];
sx q[0];
rz(-0.52310399) q[0];
sx q[0];
rz(2.4260957) q[0];
x q[1];
rz(-1.4619751) q[2];
sx q[2];
rz(-0.76078712) q[2];
sx q[2];
rz(1.6016122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24040996) q[1];
sx q[1];
rz(-2.6844271) q[1];
sx q[1];
rz(2.9967876) q[1];
rz(-pi) q[2];
x q[2];
rz(0.02324795) q[3];
sx q[3];
rz(-0.90475268) q[3];
sx q[3];
rz(-0.2153399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4652319) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(-2.5993627) q[2];
rz(-0.17255653) q[3];
sx q[3];
rz(-1.6514643) q[3];
sx q[3];
rz(1.1057314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.3365823) q[0];
sx q[0];
rz(-1.35291) q[0];
sx q[0];
rz(-1.2611058) q[0];
rz(2.005596) q[1];
sx q[1];
rz(-1.1120064) q[1];
sx q[1];
rz(-3.0109308) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9016584) q[0];
sx q[0];
rz(-1.8404418) q[0];
sx q[0];
rz(-1.1872227) q[0];
rz(-pi) q[1];
rz(-1.7950141) q[2];
sx q[2];
rz(-1.0750689) q[2];
sx q[2];
rz(2.5350646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9798292) q[1];
sx q[1];
rz(-1.5391304) q[1];
sx q[1];
rz(-2.0791441) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9404802) q[3];
sx q[3];
rz(-0.94025984) q[3];
sx q[3];
rz(2.4760404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38264349) q[2];
sx q[2];
rz(-0.073315695) q[2];
sx q[2];
rz(-0.22155133) q[2];
rz(2.4394636) q[3];
sx q[3];
rz(-1.5789072) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5704983) q[0];
sx q[0];
rz(-1.2449188) q[0];
sx q[0];
rz(-0.13394295) q[0];
rz(2.7722923) q[1];
sx q[1];
rz(-1.9422653) q[1];
sx q[1];
rz(-1.7514924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1655859) q[0];
sx q[0];
rz(-1.0635785) q[0];
sx q[0];
rz(2.8602703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8758135) q[2];
sx q[2];
rz(-1.3172564) q[2];
sx q[2];
rz(0.60088742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.578957) q[1];
sx q[1];
rz(-1.0308415) q[1];
sx q[1];
rz(-1.8946429) q[1];
rz(-1.8769299) q[3];
sx q[3];
rz(-1.9189034) q[3];
sx q[3];
rz(-0.26012421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8451346) q[2];
sx q[2];
rz(-2.239581) q[2];
sx q[2];
rz(1.585539) q[2];
rz(-0.29160094) q[3];
sx q[3];
rz(-2.6652938) q[3];
sx q[3];
rz(2.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.6911102) q[0];
sx q[0];
rz(-0.99523681) q[0];
sx q[0];
rz(0.52571785) q[0];
rz(-2.763343) q[1];
sx q[1];
rz(-1.5317761) q[1];
sx q[1];
rz(0.27413109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7129242) q[0];
sx q[0];
rz(-0.98085058) q[0];
sx q[0];
rz(-1.601786) q[0];
rz(-0.81118213) q[2];
sx q[2];
rz(-2.5497782) q[2];
sx q[2];
rz(2.9331911) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25933274) q[1];
sx q[1];
rz(-1.8253528) q[1];
sx q[1];
rz(-1.987129) q[1];
rz(-2.0977705) q[3];
sx q[3];
rz(-0.29114215) q[3];
sx q[3];
rz(0.038719161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5963001) q[2];
sx q[2];
rz(-1.9175074) q[2];
sx q[2];
rz(-2.8833585) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(-2.7365007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79893583) q[0];
sx q[0];
rz(-2.6234143) q[0];
sx q[0];
rz(3.0361191) q[0];
rz(-0.27490973) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(-2.8555433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29198439) q[0];
sx q[0];
rz(-1.7546856) q[0];
sx q[0];
rz(-1.6894132) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83467612) q[2];
sx q[2];
rz(-1.3501985) q[2];
sx q[2];
rz(-2.5479864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5040999) q[1];
sx q[1];
rz(-1.5063725) q[1];
sx q[1];
rz(-1.9570051) q[1];
rz(-pi) q[2];
rz(-2.0155573) q[3];
sx q[3];
rz(-2.0789903) q[3];
sx q[3];
rz(-1.944384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1470571) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(2.6915754) q[2];
rz(-0.58722812) q[3];
sx q[3];
rz(-1.8000032) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.29379544) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(1.1867123) q[0];
rz(-0.31934357) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(-0.26225463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569993) q[0];
sx q[0];
rz(-2.5525064) q[0];
sx q[0];
rz(2.7688945) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9840711) q[2];
sx q[2];
rz(-1.4461293) q[2];
sx q[2];
rz(-1.1050129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1407518) q[1];
sx q[1];
rz(-2.5927564) q[1];
sx q[1];
rz(1.728265) q[1];
rz(-pi) q[2];
rz(-0.75388925) q[3];
sx q[3];
rz(-2.5544832) q[3];
sx q[3];
rz(1.0894437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.78476) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-2.812815) q[2];
rz(1.5152991) q[3];
sx q[3];
rz(-1.7833775) q[3];
sx q[3];
rz(-2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560155) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(-0.27895862) q[0];
rz(0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(1.4252211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12869975) q[0];
sx q[0];
rz(-2.1919247) q[0];
sx q[0];
rz(-1.5170044) q[0];
rz(-0.89377706) q[2];
sx q[2];
rz(-2.2502459) q[2];
sx q[2];
rz(-2.627934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60630262) q[1];
sx q[1];
rz(-2.6043476) q[1];
sx q[1];
rz(0.37588889) q[1];
x q[2];
rz(0.29028671) q[3];
sx q[3];
rz(-0.93433524) q[3];
sx q[3];
rz(-0.94062128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7384537) q[2];
sx q[2];
rz(-0.29610115) q[2];
sx q[2];
rz(-2.9816755) q[2];
rz(-0.81965172) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(-0.73614365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6920456) q[0];
sx q[0];
rz(-1.8501546) q[0];
sx q[0];
rz(2.8458169) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(-1.5089418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49182941) q[0];
sx q[0];
rz(-1.2163645) q[0];
sx q[0];
rz(1.9948122) q[0];
x q[1];
rz(1.9231173) q[2];
sx q[2];
rz(-1.5218456) q[2];
sx q[2];
rz(-1.5201598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9450284) q[1];
sx q[1];
rz(-1.5032282) q[1];
sx q[1];
rz(3.0516207) q[1];
x q[2];
rz(-2.4911777) q[3];
sx q[3];
rz(-2.5456508) q[3];
sx q[3];
rz(-0.015344674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8136924) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(-2.5076765) q[2];
rz(3.0564803) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(0.77272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2278628) q[0];
sx q[0];
rz(-1.5999404) q[0];
sx q[0];
rz(3.0318442) q[0];
rz(1.2211424) q[1];
sx q[1];
rz(-2.2962062) q[1];
sx q[1];
rz(-2.5795945) q[1];
rz(0.22460266) q[2];
sx q[2];
rz(-1.1418268) q[2];
sx q[2];
rz(1.8662966) q[2];
rz(1.4681592) q[3];
sx q[3];
rz(-2.9087421) q[3];
sx q[3];
rz(-1.4431492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
