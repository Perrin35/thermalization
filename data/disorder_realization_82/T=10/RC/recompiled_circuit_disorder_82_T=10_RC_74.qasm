OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(-1.27966) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0227538) q[0];
sx q[0];
rz(-1.7099329) q[0];
sx q[0];
rz(-1.820178) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42135294) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(-2.8533964) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0496088) q[1];
sx q[1];
rz(-0.99124747) q[1];
sx q[1];
rz(2.0642573) q[1];
rz(0.23366191) q[3];
sx q[3];
rz(-1.3705001) q[3];
sx q[3];
rz(1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0063643) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(2.546229) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90855234) q[0];
sx q[0];
rz(-1.6516343) q[0];
sx q[0];
rz(-0.064706133) q[0];
x q[1];
rz(2.3357453) q[2];
sx q[2];
rz(-2.5666551) q[2];
sx q[2];
rz(2.2573059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65709719) q[1];
sx q[1];
rz(-2.8946981) q[1];
sx q[1];
rz(-2.2730278) q[1];
x q[2];
rz(1.3120193) q[3];
sx q[3];
rz(-0.32734713) q[3];
sx q[3];
rz(-0.60206383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(0.22274676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(-1.1654277) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9420625) q[2];
sx q[2];
rz(-1.7811333) q[2];
sx q[2];
rz(2.631275) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45914868) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(-3.0241806) q[1];
rz(-pi) q[2];
rz(-1.907903) q[3];
sx q[3];
rz(-1.3385696) q[3];
sx q[3];
rz(0.043134886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5266787) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(-2.148669) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6572861) q[2];
sx q[2];
rz(-1.2109586) q[2];
sx q[2];
rz(-0.53777018) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11776609) q[1];
sx q[1];
rz(-0.95925602) q[1];
sx q[1];
rz(-1.0002506) q[1];
rz(-pi) q[2];
rz(2.782015) q[3];
sx q[3];
rz(-1.3782856) q[3];
sx q[3];
rz(-2.0457091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.2801923) q[2];
rz(-0.81930339) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.357648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4015409) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(0.037365035) q[0];
rz(1.3327417) q[2];
sx q[2];
rz(-0.14698725) q[2];
sx q[2];
rz(0.57839314) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9822426) q[1];
sx q[1];
rz(-2.9061926) q[1];
sx q[1];
rz(1.2118641) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1702607) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0673163) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(-2.5308454) q[0];
rz(-3.0797708) q[2];
sx q[2];
rz(-2.0197868) q[2];
sx q[2];
rz(-0.87473727) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7102393) q[1];
sx q[1];
rz(-2.2696886) q[1];
sx q[1];
rz(-0.61183521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38576491) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(1.3145087) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(-1.1869173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884739) q[0];
sx q[0];
rz(-2.0679727) q[0];
sx q[0];
rz(0.32062809) q[0];
x q[1];
rz(0.059139472) q[2];
sx q[2];
rz(-1.5416317) q[2];
sx q[2];
rz(-1.189122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4759051) q[1];
sx q[1];
rz(-0.74531065) q[1];
sx q[1];
rz(1.3728632) q[1];
rz(0.39077057) q[3];
sx q[3];
rz(-0.6676538) q[3];
sx q[3];
rz(0.30927502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26178965) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(-1.6342182) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(-2.142895) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-2.2858455) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(-2.1112679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949769) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(1.0975518) q[0];
rz(-pi) q[1];
rz(-1.4872929) q[2];
sx q[2];
rz(-0.52831542) q[2];
sx q[2];
rz(-0.63937843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.2703018) q[1];
rz(-2.2780667) q[3];
sx q[3];
rz(-0.72609767) q[3];
sx q[3];
rz(1.8973779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(3.0714152) q[2];
rz(-2.3146546) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(1.0519741) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-2.0064328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90056706) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(-2.7450949) q[0];
rz(-pi) q[1];
x q[1];
rz(0.029149292) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(-1.0708035) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4371722) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(-pi) q[2];
rz(-1.6039861) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(-2.8896796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(-2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(2.8519582) q[3];
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
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(-0.19616729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553186) q[0];
sx q[0];
rz(-2.1325169) q[0];
sx q[0];
rz(-0.67027153) q[0];
rz(-pi) q[1];
rz(1.4354544) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(0.70336715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13388053) q[1];
sx q[1];
rz(-0.73426437) q[1];
sx q[1];
rz(1.5937362) q[1];
rz(-pi) q[2];
rz(-2.4862643) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(0.67656803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(-0.23399227) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-3.0951981) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
