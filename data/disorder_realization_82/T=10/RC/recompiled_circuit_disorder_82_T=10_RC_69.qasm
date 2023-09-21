OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3368971) q[0];
sx q[0];
rz(-2.1043632) q[0];
sx q[0];
rz(-0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5582433) q[0];
sx q[0];
rz(-1.3238751) q[0];
sx q[0];
rz(-2.9980744) q[0];
x q[1];
rz(-0.42135294) q[2];
sx q[2];
rz(-1.7755277) q[2];
sx q[2];
rz(-2.8533964) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0496088) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(1.0773354) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7765331) q[3];
sx q[3];
rz(-1.3418901) q[3];
sx q[3];
rz(-2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(2.485086) q[2];
rz(0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-2.6541236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90855234) q[0];
sx q[0];
rz(-1.4899583) q[0];
sx q[0];
rz(3.0768865) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7198118) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(0.031907206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0812379) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(-2.9802122) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2536212) q[3];
sx q[3];
rz(-1.4884236) q[3];
sx q[3];
rz(-1.2143283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(2.0987434) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449126) q[0];
sx q[0];
rz(-1.1059462) q[0];
sx q[0];
rz(0.19038598) q[0];
rz(-0.12292513) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(-2.9188459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6762786) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(-2.7200384) q[0];
x q[1];
rz(2.1026956) q[2];
sx q[2];
rz(-2.7173018) q[2];
sx q[2];
rz(1.588856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45914868) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(-0.1174121) q[1];
x q[2];
rz(1.907903) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(0.043134886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(-2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(1.4933043) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677778) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(-0.52315229) q[0];
rz(-pi) q[1];
rz(0.36107365) q[2];
sx q[2];
rz(-1.6517342) q[2];
sx q[2];
rz(-1.0635478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3355545) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(-0.6946509) q[1];
rz(-0.50587378) q[3];
sx q[3];
rz(-2.7357091) q[3];
sx q[3];
rz(-2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(-1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(-2.3732896) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4015409) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(3.1042276) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(-0.81894433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0799775) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(1.3498989) q[1];
x q[2];
rz(-2.1987678) q[3];
sx q[3];
rz(-2.4403265) q[3];
sx q[3];
rz(3.0758465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0129464) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7662402) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(1.2067632) q[0];
rz(-1.4432625) q[2];
sx q[2];
rz(-0.45293929) q[2];
sx q[2];
rz(2.4085101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39735079) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(2.1703297) q[1];
x q[2];
rz(-2.2699039) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(0.20760078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.05904077) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(-3.011076) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(-0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(0.39757279) q[0];
rz(-1.827084) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.9546753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163527) q[0];
sx q[0];
rz(-1.851474) q[0];
sx q[0];
rz(-2.0902082) q[0];
rz(-1.5415807) q[2];
sx q[2];
rz(-1.511682) q[2];
sx q[2];
rz(2.7616449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4759051) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(-1.3728632) q[1];
x q[2];
rz(-2.5116634) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26178965) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-1.0236434) q[1];
sx q[1];
rz(1.0303248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4521342) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(2.336691) q[0];
x q[1];
rz(2.0975935) q[2];
sx q[2];
rz(-1.5287405) q[2];
sx q[2];
rz(0.85925697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47146637) q[1];
sx q[1];
rz(-1.2738436) q[1];
sx q[1];
rz(-0.15836497) q[1];
x q[2];
rz(-0.97708948) q[3];
sx q[3];
rz(-1.1247375) q[3];
sx q[3];
rz(2.2462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6475911) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-1.1351599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90056706) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(0.39649773) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1124434) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(-1.0708035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2444297) q[1];
sx q[1];
rz(-2.0532554) q[1];
sx q[1];
rz(-0.06628118) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5376066) q[3];
sx q[3];
rz(-2.5182708) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90821663) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88403945) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(0.89417017) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7061383) q[2];
sx q[2];
rz(-1.6207098) q[2];
sx q[2];
rz(0.70336715) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10298221) q[1];
sx q[1];
rz(-0.83676941) q[1];
sx q[1];
rz(0.02070133) q[1];
x q[2];
rz(-0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-0.47714864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(0.23399227) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(-1.3410939) q[3];
sx q[3];
rz(-1.5256186) q[3];
sx q[3];
rz(-1.516173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];